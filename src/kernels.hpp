#ifndef PBF_KERNELS_HPP
#define PBF_KERNELS_HPP

#include "types.hpp"

Scalar calcPoly6Kernel(const Vec3& r, const Scalar h);
Vec3   calcGradPoly6Kernel(const Vec3& r, const Scalar h);

Scalar calcSpikyKernel(const Vec3& r, const Scalar h);
Vec3   calcGradSpikyKernel(const Vec3& r, const Scalar h);

#endif
