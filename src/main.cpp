#include "particle.hpp"
#include "particles-alembic-manager.hpp"
#include <Eigen/Core>
#include <cstdint>
#include <iostream>
#include <vector>

void step(std::vector<Particle>& particles)
{
    const int num_particles = particles.size();

    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].x = Vec3::Random() + Vec3(0.0, 0.5, 0.0);
    }
}

int main()
{
    std::vector<Particle> particles;

    constexpr Scalar dt = 1.0 / 60.0;

    constexpr int x_size = 6;
    constexpr int y_size = 6;
    constexpr int z_size = 6;

    // Generate and initialize particles
    constexpr int num_particles = x_size * y_size * z_size;
    particles.resize(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].i = i;
        particles[i].m = 1.0 / static_cast<Scalar>(num_particles);
        particles[i].x = Vec3::Random() + Vec3(0.0, 0.5, 0.0);
        particles[i].v = Vec3::Zero();
    }

    // Instantiate an alembic manager and submit the initial status
    ParticlesAlembicManager alembic_manager("./test.abc", dt, "fluid", &particles);
    alembic_manager.submitCurrentStatus();

    // Simulate particles
    for (int t = 0; t < 120; ++t)
    {
        // Step the simulation time
        step(particles);

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
