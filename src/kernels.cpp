#include "kernels.hpp"

Scalar calcPoly6Kernel(const Vec3& r, const Scalar h)
{
    constexpr Scalar pi    = 3.14159265358979323;
    constexpr Scalar coeff = 315.0 / (64.0 * pi);

    const Scalar h_squared = h * h;
    const Scalar r_squared = r.squaredNorm();

    if (r_squared > h_squared) {
        return 0.0;
    }

    const Scalar h_4th_power = h_squared * h_squared;
    const Scalar h_9th_power = h_4th_power * h_4th_power * h;
    const Scalar diff        = h_squared - r_squared;
    const Scalar diff_cubed  = diff * diff * diff;

    return (coeff / h_9th_power) * diff_cubed;
}

Vec3 calcGradPoly6Kernel(const Vec3& r, const Scalar h)
{
    constexpr Scalar pi    = 3.14159265358979323;
    constexpr Scalar coeff = 945.0 / (32.0 * pi);

    const Scalar h_squared = h * h;
    const Scalar r_squared = r.squaredNorm();

    if (r_squared > h_squared) {
        return Vec3::Zero();
    }

    const Scalar h_4th_power = h_squared * h_squared;
    const Scalar h_9th_power = h_4th_power * h_4th_power * h;

    const Scalar diff         = h_squared - r_squared;
    const Scalar diff_squared = diff * diff;

    return -r * (coeff / h_9th_power) * diff_squared;
}

Scalar calcSpikyKernel(const Vec3& r, const Scalar h)
{
    constexpr Scalar pi    = 3.14159265358979323;
    constexpr Scalar coeff = 15.0 / pi;

    const Scalar r_norm = r.norm();

    if (r_norm > h) {
        return 0.0;
    }

    const Scalar h_cubed     = h * h * h;
    const Scalar h_6th_power = h_cubed * h_cubed;
    const Scalar diff        = h - r_norm;
    const Scalar diff_cubed  = diff * diff * diff;

    return (coeff / h_6th_power) * diff_cubed;
}

Vec3 calcGradSpikyKernel(const Vec3& r, const Scalar h)
{
    constexpr Scalar pi    = 3.14159265358979323;
    constexpr Scalar coeff = 45.0 / pi;

    const Scalar r_norm = r.norm();

    if (r_norm > h) {
        return Vec3::Zero();
    }

    const Scalar h_cubed      = h * h * h;
    const Scalar h_6th_power  = h_cubed * h_cubed;
    const Scalar diff         = h - r_norm;
    const Scalar diff_squared = diff * diff;

    return -r * (coeff / (h_6th_power * std::max(r_norm, 1e-24))) * diff_squared;
}
