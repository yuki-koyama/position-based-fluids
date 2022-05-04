#ifndef PBF_KERNELS_HPP
#define PBF_KERNELS_HPP

#include "types.hpp"

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

#endif
