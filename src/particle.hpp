#ifndef PBF_PARTICLE_HPP
#define PBF_PARTICLE_HPP

#include <Eigen/Core>
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

#endif
