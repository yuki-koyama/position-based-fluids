#include <Eigen/Core>
#include <iostream>
#include <vector>

using Vec3 = Eigen::Vector3d;

struct Particle
{
    uint32_t i;
    double   m;
    Vec3     x;
    Vec3     v;

    Vec3                  p;
    std::vector<uint32_t> neighbors;
};

int main()
{
    std::vector<Particle> particles{20};

    for (auto& particle : particles)
    {
        std::cout << particle.x.transpose() << std::endl;
    }

    return 0;
}
