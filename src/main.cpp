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
            // Do not include the target particle itself as its neighborhoods
            if (i == j)
            {
                continue;
            }

            const Scalar squared_dist = (particles[i].p - particles[j].p).squaredNorm();

            if (squared_dist < radius_squared)
            {
                neighbors_list[i].push_back(j);
            }
        }
    }

    return neighbors_list;
}

Scalar calcConstraint(const int                            target_index,
                      const std::vector<Particle>&         particles,
                      const std::vector<std::vector<int>>& neighbor_list,
                      const Scalar                         rest_density,
                      const Scalar                         radius)
{
    const auto& p_target = particles[target_index];

    Scalar density = 0.0;
    for (int neighbor_index : neighbor_list[target_index])
    {
        const auto& p = particles[neighbor_index];

        density += p.m * calcKernel(p_target.p - p.p, radius);
    }

    return (density / rest_density) - 1.0;
}

Vec3 calcGradConstraint(const int                            target_index,
                        const int                            var_index,
                        const std::vector<Particle>&         particles,
                        const std::vector<std::vector<int>>& neighbor_list,
                        const Scalar                         rest_density,
                        const Scalar                         radius)
{
    // TODO: Calculate gradient
    return Vec3::Zero();
}

void step(const Scalar dt, std::vector<Particle>& particles)
{
    constexpr int    num_iters    = 10;
    constexpr Scalar radius       = 0.2;
    constexpr Scalar rest_density = 1000.0;

    const int num_particles = particles.size();

    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].v = particles[i].v + dt * Vec3(0.0, -9.8, 0.0);
        particles[i].p = particles[i].x + dt * particles[i].v;
    }

    const auto neighbors_list = findNeighborParticles(radius, particles);

    for (int k = 0; k < num_iters; ++k)
    {
        // Calculate lambda
        VecX lambda(num_particles);
        for (int i = 0; i < num_particles; ++i)
        {
            auto& p = particles[i];

            // TODO: Calculate lambda here
            lambda[i] = 0.0;
        }

        // Calculate delta p (note: Jacobi style)
        MatX delta_p(3, num_particles);
        for (int i = 0; i < num_particles; ++i)
        {
            auto& p = particles[i];

            const int num_neighbors = neighbors_list[i].size();

            // If there is no neighborhoods, just ignore the pressure
            if (num_neighbors == 0)
            {
                continue;
            }

            // Calculate the sum of pressure effect
            MatX buffer(3, num_neighbors);
            for (int j = 0; j < num_neighbors; ++j)
            {
                const int    neighbor_index = neighbors_list[i][j];
                const Scalar coeff          = lambda[i] + lambda[neighbor_index];

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

int main()
{
    std::vector<Particle> particles;

    constexpr Scalar dt = 1.0 / 60.0;

    constexpr int x_size = 8;
    constexpr int y_size = 8;
    constexpr int z_size = 8;

    // Generate and initialize particles
    constexpr int num_particles = x_size * y_size * z_size;
    particles.resize(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].i = i;
        particles[i].m = 1.0 / static_cast<Scalar>(num_particles);
        particles[i].x = Vec3::Random() + Vec3(0.0, 1.0, 0.0);
        particles[i].v = Vec3::Zero();
    }

    // Instantiate an alembic manager and submit the initial status
    ParticlesAlembicManager alembic_manager("./test.abc", dt, "Fluid", &particles);
    alembic_manager.submitCurrentStatus();

    // Simulate particles
    for (int t = 0; t < 60; ++t)
    {
        // Step the simulation time
        step(dt, particles);

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
