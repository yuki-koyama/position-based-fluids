#include <Eigen/Core>
#include <iostream>
#include <vector>

using Vec3 = Eigen::Vector3d;

struct Particle
{
    double m;
    Vec3   x;
    Vec3   p;
    Vec3   v;
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
