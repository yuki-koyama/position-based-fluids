#ifndef PBF_PARTICLE_HPP
#define PBF_PARTICLE_HPP

#include "types.hpp"

struct Particle
{
    int    i;
    Scalar m;
    Vec3   x;
    Vec3   v;
    Vec3   p;
};

#endif
