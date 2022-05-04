#include "particle.hpp"
#include "particles-alembic-manager.hpp"
#include <Eigen/Core>
#include <cstdint>
#include <iostream>
#include <vector>

std::vector<std::vector<int>> findNeighborParticles(const Scalar radius, const std::vector<Particle>& particles)
{
    const int num_particles = particles.size();

    std::vector<std::vector<int>> neighbors_list(num_particles);

    for (int i = 0; i < num_particles; ++i)
    {
        for (int j = 0; j < num_particles; ++j)
        {
            const Scalar squared_dist = (particles[i].p - particles[j].p).squaredNorm();

            if (squared_dist < radius * radius)
            {
                neighbors_list[i].push_back(j);
            }
        }
    }

    return neighbors_list;
}

Scalar calcPoly6Kernel(const Vec3 r, const Scalar h)
{
    constexpr Scalar pi    = 3.14159265358979323;
    constexpr Scalar coeff = 315.0 / (64.0 * pi);

    const Scalar h_squared = h * h;
    const Scalar r_squared = r.squaredNorm();

    if (r.squaredNorm() > h_squared)
    {
        return 0.0;
    }

    const Scalar h_4th_power = h_squared * h_squared;
    const Scalar h_9th_power = h_4th_power * h_4th_power * h;
    const Scalar diff        = h_squared - r_squared;
    const Scalar diff_cubed  = diff * diff * diff;

    return (coeff / h_9th_power) * diff_cubed;
}

void step(const Scalar dt, std::vector<Particle>& particles)
{
    constexpr int    num_iters    = 10;
    constexpr Scalar radius       = 0.2;
    constexpr Scalar rest_density = 1.0;

    const int num_particles = particles.size();

    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].v = particles[i].v + dt * Vec3(0.0, -9.8, 0.0);
        particles[i].p = particles[i].x + dt * particles[i].v;
    }

    const auto neighbors_list = findNeighborParticles(radius, particles);

    for (int k = 0; k < num_iters; ++k)
    {
        for (int i = 0; i < num_particles; ++i)
        {
            auto& p = particles[i];

            // Detect and handle environmental collisions (in a very naive way)
            p.p = p.p.cwiseMax(Vec3(-1.0, 0.0, -1.0));
            p.p = p.p.cwiseMin(Vec3(+1.0, 2.0, +1.0));
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
    ParticlesAlembicManager alembic_manager("./test.abc", dt, "fluid", &particles);
    alembic_manager.submitCurrentStatus();

    // Simulate particles
    for (int t = 0; t < 120; ++t)
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
