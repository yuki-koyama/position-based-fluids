#include "kernels.hpp"
#include "neighbor-search-engine.hpp"
#include "particle.hpp"
#include "particles-alembic-manager.hpp"
#include "types.hpp"
#include <iostream>
#include <parallel-util.hpp>
#include <timer.hpp>
#include <vector>

constexpr auto calcKernel     = calcPoly6Kernel;
constexpr auto calcGradKernel = calcGradSpikyKernel;

using NeighborSearchEngine = UniformGridNeighborSearchEngine;

Scalar calcDensity(const int                    target_index,
                   const std::vector<Particle>& particles,
                   const NeighborSearchEngine&  neighbor_search_engine,
                   const Scalar                 radius)
{
    const auto& p_target = particles[target_index];

    Scalar density = 0.0;
    for (int neighbor_index : neighbor_search_engine.retrieveNeighbors(target_index)) {
        const auto& p = particles[neighbor_index];

        density += p.m * calcKernel(p_target.p - p.p, radius);
    }

    return density;
}

Scalar calcConstraint(const int                    target_index,
                      const std::vector<Particle>& particles,
                      const NeighborSearchEngine&  neighbor_search_engine,
                      const Scalar                 rest_density,
                      const Scalar                 radius)
{
    const Scalar density = calcDensity(target_index, particles, neighbor_search_engine, radius);

    return (density / rest_density) - 1.0;
}

Vec3 calcGradConstraint(const int                    target_index,
                        const int                    var_index,
                        const std::vector<Particle>& particles,
                        const NeighborSearchEngine&  neighbor_search_engine,
                        const Scalar                 rest_density,
                        const Scalar                 radius)
{
    const auto& p_target = particles[target_index];

    if (target_index == var_index) {
        Vec3 sum = Vec3::Zero();
        for (int neighbor_index : neighbor_search_engine.retrieveNeighbors(target_index)) {
            const auto& p = particles[neighbor_index];

            sum += p.m * calcGradKernel(p_target.p - p.p, radius);
        }

        return sum / rest_density;
    } else {
        const auto& p = particles[var_index];

        return -p.m * calcGradKernel(p_target.p - p.p, radius) / rest_density;
    }
}

void printAverageNumNeighbors(const NeighborSearchEngine& neighbor_search_engine)
{
    const int num_particles = neighbor_search_engine.getNumParticles();

    VecX nums(num_particles);
    for (int i = 0; i < num_particles; ++i) {
        nums[i] = neighbor_search_engine.retrieveNeighbors(i).size();
    }
    std::cout << "Average(#neighbors): " << nums.mean() << std::endl;
}

void printAverageDensity(const std::vector<Particle>& particles,
                         const NeighborSearchEngine&  neighbor_search_engine,
                         const Scalar                 radius)
{
    VecX buffer(particles.size());
    for (int i = 0; i < particles.size(); ++i) {
        buffer[i] = calcDensity(i, particles, neighbor_search_engine, radius);
    }
    std::cout << "Average(density): " << buffer.mean() << std::endl;
}

void step(const Scalar dt, std::vector<Particle>& particles)
{
    constexpr int    num_iters       = 2;
    constexpr Scalar radius          = 0.10;
    constexpr Scalar rest_density    = 1000.0;
    constexpr Scalar epsilon_cfm     = 1e+05;
    constexpr Scalar damping         = 0.999;
    constexpr Scalar viscosity_coeff = 0.050;
    constexpr bool   verbose         = false;

    const int num_particles = particles.size();

    for (int i = 0; i < num_particles; ++i) {
        particles[i].v = particles[i].v + dt * Vec3(0.0, -9.8, 0.0);
        particles[i].p = particles[i].x + dt * particles[i].v;
    }

    // Prepare a neighborhood search engine
    NeighborSearchEngine neighbor_search_engine(radius, particles);

    // Find neighborhoods of every particle
    neighbor_search_engine.searchNeighbors();

    if constexpr (verbose) {
        printAverageNumNeighbors(neighbor_search_engine);
        printAverageDensity(particles, neighbor_search_engine, radius);
    }

    for (int k = 0; k < num_iters; ++k) {
        // Calculate lambda
        VecX lambda(num_particles);

        const auto calc_lambda = [&](const int i) {
            const auto&  p         = particles[i];
            const Scalar numerator = calcConstraint(i, particles, neighbor_search_engine, rest_density, radius);

            Scalar denominator = 0.0;
            for (int neighbor_index : neighbor_search_engine.retrieveNeighbors(i)) {
                const Vec3 grad =
                    calcGradConstraint(i, neighbor_index, particles, neighbor_search_engine, rest_density, radius);

                // Note: In Eq.12, the inverse mass is dropped for simplicity
                denominator += (1.0 / particles[neighbor_index].m) * grad.squaredNorm();
            }
            // Note: Add an epsilon value for relaxation (see Eq.11)
            // TODO: Confirm this equation
            denominator += epsilon_cfm;

            lambda[i] = -numerator / denominator;
        };
        parallelutil::parallel_for(num_particles, calc_lambda);

        // Calculate delta p (note: Jacobi style)
        MatX delta_p(3, num_particles);

        const auto calc_delta_p = [&](const int i) {
            const auto& p             = particles[i];
            const auto& neighbors     = neighbor_search_engine.retrieveNeighbors(i);
            const int   num_neighbors = neighbors.size();

            // Calculate the artificial tensile pressure correction constant
            constexpr Scalar corr_n = 4.0;
            constexpr Scalar corr_h = 0.30;
            const Scalar     corr_k = p.m * 1.0e-04; // Note: This equation has no ground and may not work well
            const Scalar     corr_w = calcKernel(corr_h * radius * Vec3::UnitX(), radius);

            // Calculate the sum of pressure effect (Eq.12)
            MatX buffer(3, num_neighbors);
            for (int j = 0; j < num_neighbors; ++j) {
                const int neighbor_index = neighbors[j];

                // Calculate the artificial tensile pressure correction
                const Scalar kernel_val = calcKernel(p.p - particles[neighbor_index].p, radius);
                const Scalar ratio      = kernel_val / corr_w;
                const Scalar corr_coeff = -corr_k * std::pow(ratio, corr_n);

                const Scalar coeff = particles[neighbor_index].m * (lambda[i] + lambda[neighbor_index] + corr_coeff);

                buffer.col(j) = coeff * calcGradKernel(p.p - particles[neighbor_index].p, radius);
            }
            const Vec3 sum = buffer.rowwise().sum();

            // Calculate delta p of this particle
            delta_p.col(i) = (1.0 / p.m) * (1.0 / rest_density) * sum;
        };
        parallelutil::parallel_for(num_particles, calc_delta_p);

        // Apply delta p in the Jacobi style
        const auto apply_delta_p = [&](const int i) {
            particles[i].p += delta_p.col(i);
        };
        parallelutil::parallel_for(num_particles, apply_delta_p);

        for (int i = 0; i < num_particles; ++i) {
            auto& p = particles[i];

            // Detect and handle environmental collisions (in a very naive way)
            p.p = p.p.cwiseMax(Vec3(-1.0, 0.0, -1.0));
            p.p = p.p.cwiseMin(Vec3(+1.0, 8.0, +1.0));
        }
    }

    // Update positions and velocities
    for (int i = 0; i < num_particles; ++i) {
        particles[i].v = damping * (particles[i].p - particles[i].x) / dt;
        particles[i].x = particles[i].p;
    }

    // Apply the XSPH viscosity effect
    VecX densities(num_particles);
    MatX delta_v(3, num_particles);

    parallelutil::parallel_for(num_particles, [&](const int i) {
        densities[i] = calcDensity(i, particles, neighbor_search_engine, radius);
    });

    parallelutil::parallel_for(num_particles, [&](const int i) {
        const auto& p             = particles[i];
        const auto& neighbors     = neighbor_search_engine.retrieveNeighbors(i);
        const int   num_neighbors = neighbors.size();

        MatX buffer(3, num_neighbors);
        for (int j = 0; j < num_neighbors; ++j) {
            const int    neighbor_index = neighbors[j];
            const Scalar kernel_val     = calcKernel(p.x - particles[neighbor_index].x, radius);
            const auto   rel_velocity   = particles[neighbor_index].v - p.v;

            buffer.col(j) = (p.m / densities[neighbor_index]) * kernel_val * rel_velocity;
        }
        const auto sum = buffer.rowwise().sum();

        delta_v.col(i) = viscosity_coeff * sum;
    });

    parallelutil::parallel_for(num_particles, [&](const int i) {
        particles[i].v += delta_v.col(i);
    });
}

int main()
{
    std::vector<Particle> particles;

    constexpr Scalar dt            = 1.0 / 60.0;
    constexpr int    num_frames    = 240;
    constexpr int    num_particles = 108000;

    std::cout << "#particles: " << num_particles << std::endl;

    // Generate and initialize particles
    particles.resize(num_particles);
    for (int i = 0; i < num_particles; ++i) {
        particles[i].i = i;
        particles[i].m = 3000.0 / static_cast<Scalar>(num_particles);
        particles[i].x = Vec3(0.5, 1.5, 0.5).cwiseProduct(Vec3::Random()) + Vec3(-0.5, 2.5, 0.0);
        particles[i].v = Vec3::Zero();
    }

    // Instantiate an alembic manager and submit the initial status
    ParticlesAlembicManager alembic_manager("./test.abc", dt, "Fluid", &particles);

    // Relax initial particle positions (a dirty solution for resolving bad initial states)
    constexpr int num_relax_steps = 20;
    for (int k = 0; k < num_relax_steps; ++k) {
        const auto   timer   = timer::Timer("Init " + std::to_string(k) + " / " + std::to_string(num_relax_steps));
        const Scalar damping = std::min(1.0, (static_cast<Scalar>(k) * 2.0 / static_cast<Scalar>(num_relax_steps)));

        // Step the simulation time forward using a very small time step
        step(1e-04 * dt, particles);

        // Damp velocities for stability
        for (auto& p : particles) {
            p.v = damping * p.v;
        }
    }
    for (auto& p : particles) {
        p.v = Vec3::Zero();
    }

    // Write the current status
    alembic_manager.submitCurrentStatus();

    // Simulate particles
    constexpr int    num_substeps = 5;
    constexpr Scalar sub_dt       = dt / static_cast<Scalar>(num_substeps);
    for (int t = 0; t < num_frames; ++t) {
        // Instantiate timer
        const auto timer = timer::Timer("Frame #" + std::to_string(t));

        // Step the simulation time
        for (int k = 0; k < num_substeps; ++k) {
            step(sub_dt, particles);
        }

        // Write the current status
        alembic_manager.submitCurrentStatus();
    }

    return 0;
}
