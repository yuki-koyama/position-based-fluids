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
    std::vector<Particle> particles{20};

    for (auto& particle : particles)
    {
        std::cout << particle.x.transpose() << std::endl;
    }

    return 0;
}
