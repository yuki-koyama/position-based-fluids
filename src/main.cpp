#include "particle.hpp"
#include "particles-alembic-manager.hpp"
#include <Eigen/Core>
#include <cstdint>
#include <iostream>
#include <vector>

void step(const Scalar dt, std::vector<Particle>& particles)
{
    constexpr int num_iters = 10;

    const int num_particles = particles.size();

    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].v = particles[i].v + dt * Vec3(0.0, -9.8, 0.0);
        particles[i].p = particles[i].x + dt * particles[i].v;
    }

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
