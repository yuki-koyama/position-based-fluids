#include <Eigen/Core>
#include <cstdint>
#include <iostream>
#include <vector>

using Scalar = double;
using Vec3   = Eigen::Matrix<Scalar, 3, 1>;

struct Particle
{
    int    i;
    Scalar m;
    Vec3   x;
    Vec3   v;

    Vec3             p;
    std::vector<int> neighbors;
};

int main()
{
    std::vector<Particle> particles;

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
        particles[i].x = Vec3::Random();
        particles[i].v = Vec3::Zero();
    }

    // Print particle status
    for (auto& particle : particles)
    {
        std::cout << particle.x.transpose() << std::endl;
    }

    return 0;
}
