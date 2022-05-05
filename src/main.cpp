#include "kernels.hpp"
#include "particle.hpp"
#include "particles-alembic-manager.hpp"
#include "types.hpp"
#include <Eigen/Core>
#include <cstdint>
#include <iostream>
#include <vector>

constexpr auto calcKernel     = calcPoly6Kernel;
constexpr auto calcGradKernel = calcGradPoly6Kernel;

std::vector<std::vector<int>> findNeighborParticles(const Scalar radius, const std::vector<Particle>& particles)
{
    const int    num_particles  = particles.size();
    const Scalar radius_squared = radius * radius;

    std::vector<std::vector<int>> neighbors_list(num_particles);

    for (int i = 0; i < num_particles; ++i)
    {
        for (int j = 0; j < num_particles; ++j)
        {
            const Scalar squared_dist = (particles[i].p - particles[j].p).squaredNorm();

            if (squared_dist < radius_squared)
            {
                neighbors_list[i].push_back(j);
            }
        }
    }

    return neighbors_list;
}

Scalar calcDensity(const int                            target_index,
                   const std::vector<Particle>&         particles,
                   const std::vector<std::vector<int>>& neighbor_list,
                   const Scalar                         radius)
{
    const auto& p_target = particles[target_index];

    Scalar density = 0.0;
    for (int neighbor_index : neighbor_list[target_index])
    {
        const auto& p = particles[neighbor_index];

        density += p.m * calcKernel(p_target.p - p.p, radius);
    }

    return density;
}

Scalar calcConstraint(const int                            target_index,
                      const std::vector<Particle>&         particles,
                      const std::vector<std::vector<int>>& neighbor_list,
                      const Scalar                         rest_density,
                      const Scalar                         radius)
{
    const Scalar density = calcDensity(target_index, particles, neighbor_list, radius);

    return (density / rest_density) - 1.0;
}

Vec3 calcGradConstraint(const int                            target_index,
                        const int                            var_index,
                        const std::vector<Particle>&         particles,
                        const std::vector<std::vector<int>>& neighbor_list,
                        const Scalar                         rest_density,
                        const Scalar                         radius)
{
    const auto& p_target = particles[target_index];

    if (target_index == var_index)
    {
        Vec3 sum = Vec3::Zero();
        for (int neighbor_index : neighbor_list[target_index])
        {
            const auto& p = particles[neighbor_index];

            sum += p.m * calcGradKernel(p_target.p - p.p, radius);
        }

        return sum / rest_density;
    }
    else
    {
        const auto& p = particles[var_index];

        return -p.m * calcGradKernel(p_target.p - p.p, radius) / rest_density;
    }
}

void printAverageNumNeighbors(const std::vector<std::vector<int>>& neighbors_list)
{
    VecX nums(neighbors_list.size());
    for (int i = 0; i < neighbors_list.size(); ++i)
    {
        nums[i] = neighbors_list[i].size();
    }
    std::cout << "Average(#neighbors): " << nums.mean() << std::endl;
}

void printAverageDensity(const std::vector<Particle>&         particles,
                         const std::vector<std::vector<int>>& neighbor_list,
                         const Scalar                         radius)
{
    VecX buffer(particles.size());
    for (int i = 0; i < particles.size(); ++ i)
    {
        buffer[i] = calcDensity(i, particles, neighbor_list, radius);
    }
    std::cout << "Average(density): " << buffer.mean() << std::endl;
}

void step(const Scalar dt, std::vector<Particle>& particles)
{
    constexpr int    num_iters    = 4;
    constexpr Scalar radius       = 0.2;
    constexpr Scalar rest_density = 1000.0;
    constexpr Scalar epsilon      = 5e+01;

    const int num_particles = particles.size();

    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].v = particles[i].v + dt * Vec3(0.0, -9.8, 0.0);
        particles[i].p = particles[i].x + dt * particles[i].v;
    }

    const auto neighbors_list = findNeighborParticles(radius, particles);

    printAverageNumNeighbors(neighbors_list);
    printAverageDensity(particles, neighbors_list, radius);

    for (int k = 0; k < num_iters; ++k)
    {
        // Calculate lambda
        VecX lambda(num_particles);
        for (int i = 0; i < num_particles; ++i)
        {
            auto& p = particles[i];

            const Scalar numerator = calcConstraint(i, particles, neighbors_list, rest_density, radius);

            Scalar denominator = 0.0;
            for (int neighbor_index : neighbors_list[i])
            {
                const Vec3 grad =
                    calcGradConstraint(i, neighbor_index, particles, neighbors_list, rest_density, radius);

                // Note: In Eq.12, the inverse mass is dropped for simplicity
                denominator += (1.0 / particles[neighbor_index].m) * grad.squaredNorm();
            }
            // Note: Add an epsilon value for relaxation (see Eq.11)
            denominator += epsilon;

            lambda[i] = -numerator / denominator;
        }

        // Calculate delta p (note: Jacobi style)
        MatX delta_p(3, num_particles);
        for (int i = 0; i < num_particles; ++i)
        {
            auto& p = particles[i];

            const int num_neighbors = neighbors_list[i].size();

            assert(num_neighbors > 0);

            // Calculate the sum of pressure effect (Eq.12)
            MatX buffer(3, num_neighbors);
            for (int j = 0; j < num_neighbors; ++j)
            {
                const int neighbor_index = neighbors_list[i][j];

                // TODO: Confirm this equation
                const Scalar coeff = particles[neighbor_index].m * lambda[i] + p.m * lambda[neighbor_index];

                buffer.col(j) = coeff * calcGradKernel(p.p - particles[neighbor_index].p, radius);
            }
            const Vec3 sum = buffer.rowwise().sum();

            // Calculate delta p of this particle
            delta_p.col(i) = (1.0 / rest_density) * sum;
        }

        // Apply delta p (note: Jacobi style)
        for (int i = 0; i < num_particles; ++i)
        {
            particles[i].p += delta_p.col(i);
        }

        for (int i = 0; i < num_particles; ++i)
        {
            auto& p = particles[i];

            // Detect and handle environmental collisions (in a very naive way)
            p.p = p.p.cwiseMax(Vec3(-1.0, 0.0, -1.0));
            p.p = p.p.cwiseMin(Vec3(+1.0, 4.0, +1.0));
        }
    }

    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].v = (particles[i].p - particles[i].x) / dt;
        particles[i].x = particles[i].p;
    }
}

void test()
{
    constexpr auto   calcFunc     = calcPoly6Kernel;
    constexpr auto   calcGradFunc = calcGradPoly6Kernel;
    constexpr Scalar epsilon      = 1e-08;
    constexpr Scalar radius       = 0.2;

    const Vec3 x = radius * Vec3::Random();

    const Scalar val = calcFunc(x, radius);

    Vec3 numerical_grad;
    for (int i = 0; i < 3; ++i)
    {
        Vec3 eps = Vec3::Zero();
        eps[i]   = epsilon;

        numerical_grad[i] = (calcFunc(x + eps, radius) - val) / epsilon;
    }

    const Vec3 grad = calcGradFunc(x, radius);

    assert((grad - numerical_grad).norm() < 1e-04);
}

int main()
{
    std::vector<Particle> particles;

    constexpr Scalar dt         = 1.0 / 60.0;
    constexpr int    num_frames = 60;

    constexpr int x_size = 20;
    constexpr int y_size = 20;
    constexpr int z_size = 20;

    // Generate and initialize particles
    constexpr int num_particles = x_size * y_size * z_size;
    particles.resize(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].i = i;
        particles[i].m = 4000.0 / static_cast<Scalar>(num_particles);
        particles[i].x = Vec3(0.5, 2.0, 0.5).cwiseProduct(Vec3::Random()) + Vec3(0.0, 2.0, 0.0);
        particles[i].v = Vec3::Zero();
    }

    // Instantiate an alembic manager and submit the initial status
    ParticlesAlembicManager alembic_manager("./test.abc", dt, "Fluid", &particles);
    alembic_manager.submitCurrentStatus();

    // Simulate particles
    for (int t = 0; t < num_frames; ++t)
    {
        constexpr int    num_substeps = 2;
        constexpr Scalar sub_dt       = dt / static_cast<Scalar>(num_substeps);

        // Step the simulation time
        for (int k = 0; k < num_substeps; ++k)
        {
            step(sub_dt, particles);
        }

        // Write the current status
        alembic_manager.submitCurrentStatus();
    }

    // Print particle status
    for (auto& particle : particles)
    {
        std::cout << particle.x.transpose() << std::endl;
    }

    return 0;
}
