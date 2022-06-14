#include "kernels.hpp"
#include "particle.hpp"
#include "particles-alembic-manager.hpp"
#include "types.hpp"
#include <Eigen/Core>
#include <cstdint>
#include <iostream>
#include <map>
#include <parallel-util.hpp>
#include <timer.hpp>
#include <vector>

constexpr auto calcKernel     = calcPoly6Kernel;
constexpr auto calcGradKernel = calcGradSpikyKernel;

constexpr int n_x = 100;
constexpr int n_y = 100;
constexpr int n_z = 100;
const Vec3    grid_center(0.0, 1.0, 0.0);
const Vec3    grid_numbers(n_x, n_y, n_z);

std::tuple<int, int, int> calcCellIndex(const Scalar radius, const MatX& positions, const int target_particle_index)
{
    const Vec3 pos            = positions.col(target_particle_index);
    const Vec3 grid_coord_pos = (pos - grid_center) * (1.0 / radius) + 0.5 * grid_numbers;

    const int i_x = std::floor(grid_coord_pos[0]);
    const int i_y = std::floor(grid_coord_pos[1]);
    const int i_z = std::floor(grid_coord_pos[2]);

    assert(i_x >= 0);
    assert(i_y >= 0);
    assert(i_z >= 0);
    assert(i_x < n_x);
    assert(i_y < n_y);
    assert(i_z < n_z);

    return std::tuple<int, int, int>{i_x, i_y, i_z};
}

int convertIndex(const std::tuple<int, int, int>& index)
{
    return std::get<0>(index) + n_x * std::get<1>(index) + n_x * n_y * std::get<2>(index);
}

std::unordered_map<int, std::vector<int>> constructGridCells(const Scalar radius, const MatX& positions)
{
    std::unordered_map<int, std::vector<int>> grid_cells;

    for (int i = 0; i < positions.cols(); ++i)
    {
        const int cell_index = convertIndex(calcCellIndex(radius, positions, i));

        if (grid_cells.find(cell_index) == grid_cells.end())
        {
            grid_cells[cell_index] = std::vector<int>{i};
        }
        else
        {
            grid_cells[cell_index].push_back(i);
        }
    }

    return grid_cells;
}

std::vector<std::vector<int>> findNeighborParticlesNaive(const Scalar radius, const std::vector<Particle>& particles)
{
    const int    num_particles  = particles.size();
    const Scalar radius_squared = radius * radius;

    std::vector<std::vector<int>> neighbors_list(num_particles);

    for (int i = 0; i < num_particles; ++i)
    {
        for (int j = 0; j < num_particles; ++j)
        {
            const Scalar squared_dist = (particles[i].p - particles[j].p).squaredNorm();

            if (squared_dist < radius_squared)
            {
                neighbors_list[i].push_back(j);
            }
        }
    }

    return neighbors_list;
}

std::vector<std::vector<int>> findNeighborParticles(const Scalar radius, const std::vector<Particle>& particles)
{
    const int    num_particles  = particles.size();
    const Scalar radius_squared = radius * radius;

    MatX positions(3, num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        positions.col(i) = particles[i].p;
    }

    auto grid_cells = constructGridCells(radius, positions);

    std::vector<std::vector<int>> neighbors_list(num_particles);

    for (int i = 0; i < num_particles; ++i)
    {
        // Determine the cell that the target particle belongs to
        const auto target_cell_index = calcCellIndex(radius, positions, i);

        // Visit the 26 neighbor cells and the cell itself (27 in total)
        for (int x : {-1, 0, 1})
        {
            for (int y : {-1, 0, 1})
            {
                for (int z : {-1, 0, 1})
                {
                    const int i_x = std::get<0>(target_cell_index) + x;
                    const int i_y = std::get<1>(target_cell_index) + y;
                    const int i_z = std::get<2>(target_cell_index) + z;

                    const int cell_index = convertIndex({i_x, i_y, i_z});

                    const auto& list = grid_cells[cell_index];

#if 1
                    neighbors_list[i].insert(neighbors_list[i].end(), list.begin(), list.end());
#else
                    for (const int& index : list)
                    {
                        const Scalar squared_dist = (particles[i].p - particles[index].p).squaredNorm();

                        if (squared_dist < 1.2 * radius_squared)
                        {
                            neighbors_list[i].push_back(index);
                        }
                    }
#endif
                }
            }
        }
    }

    return neighbors_list;
}

Scalar calcDensity(const int                            target_index,
                   const std::vector<Particle>&         particles,
                   const std::vector<std::vector<int>>& neighbor_list,
                   const Scalar                         radius)
{
    const auto& p_target = particles[target_index];

    Scalar density = 0.0;
    for (int neighbor_index : neighbor_list[target_index])
    {
        const auto& p = particles[neighbor_index];

        density += p.m * calcKernel(p_target.p - p.p, radius);
    }

    return density;
}

Scalar calcConstraint(const int                            target_index,
                      const std::vector<Particle>&         particles,
                      const std::vector<std::vector<int>>& neighbor_list,
                      const Scalar                         rest_density,
                      const Scalar                         radius)
{
    const Scalar density = calcDensity(target_index, particles, neighbor_list, radius);

    return (density / rest_density) - 1.0;
}

Vec3 calcGradConstraint(const int                            target_index,
                        const int                            var_index,
                        const std::vector<Particle>&         particles,
                        const std::vector<std::vector<int>>& neighbor_list,
                        const Scalar                         rest_density,
                        const Scalar                         radius)
{
    const auto& p_target = particles[target_index];

    if (target_index == var_index)
    {
        Vec3 sum = Vec3::Zero();
        for (int neighbor_index : neighbor_list[target_index])
        {
            const auto& p = particles[neighbor_index];

            sum += p.m * calcGradKernel(p_target.p - p.p, radius);
        }

        return sum / rest_density;
    }
    else
    {
        const auto& p = particles[var_index];

        return -p.m * calcGradKernel(p_target.p - p.p, radius) / rest_density;
    }
}

void printAverageNumNeighbors(const std::vector<std::vector<int>>& neighbors_list)
{
    VecX nums(neighbors_list.size());
    for (int i = 0; i < neighbors_list.size(); ++i)
    {
        nums[i] = neighbors_list[i].size();
    }
    std::cout << "Average(#neighbors): " << nums.mean() << std::endl;
}

void printAverageDensity(const std::vector<Particle>&         particles,
                         const std::vector<std::vector<int>>& neighbors_list,
                         const Scalar                         radius)
{
    VecX buffer(particles.size());
    for (int i = 0; i < particles.size(); ++i)
    {
        buffer[i] = calcDensity(i, particles, neighbors_list, radius);
    }
    std::cout << "Average(density): " << buffer.mean() << std::endl;
}

void step(const Scalar dt, std::vector<Particle>& particles)
{
    constexpr int    num_iters    = 2;
    constexpr Scalar radius       = 0.10;
    constexpr Scalar rest_density = 1000.0;
    constexpr Scalar epsilon_cfm  = 1e+05;
    constexpr Scalar damping      = 0.999;

    const int num_particles = particles.size();

    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].v = particles[i].v + dt * Vec3(0.0, -9.8, 0.0);
        particles[i].p = particles[i].x + dt * particles[i].v;
    }

    const auto neighbors_list = findNeighborParticles(radius, particles);

    printAverageNumNeighbors(neighbors_list);
    printAverageDensity(particles, neighbors_list, radius);

    for (int k = 0; k < num_iters; ++k)
    {
        // Calculate lambda
        VecX lambda(num_particles);

        const auto calc_lambda = [&](const int i)
        {
            const auto&  p         = particles[i];
            const Scalar numerator = calcConstraint(i, particles, neighbors_list, rest_density, radius);

            Scalar denominator = 0.0;
            for (int neighbor_index : neighbors_list[i])
            {
                const Vec3 grad =
                    calcGradConstraint(i, neighbor_index, particles, neighbors_list, rest_density, radius);

                // Note: In Eq.12, the inverse mass is dropped for simplicity
                denominator += (1.0 / particles[neighbor_index].m) * grad.squaredNorm();
            }
            // Note: Add an epsilon value for relaxation (see Eq.11)
            // TODO: Confirm this equation
            denominator += epsilon_cfm;

            lambda[i] = -numerator / denominator;
        };
        parallelutil::parallel_for(num_particles, calc_lambda);

        // Calculate delta p (note: Jacobi style)
        MatX delta_p(3, num_particles);

        const auto calc_delta_p = [&](const int i)
        {
            const auto& p             = particles[i];
            const int   num_neighbors = neighbors_list[i].size();

            // Calculate the artificial tensile pressure correction constant
            constexpr Scalar corr_n = 4.0;
            constexpr Scalar corr_h = 0.30;
            const Scalar     corr_k = p.m * 1e-04;
            const Scalar     corr_w = calcKernel(corr_h * radius * Vec3::UnitX(), radius);

            // Calculate the sum of pressure effect (Eq.12)
            MatX buffer(3, num_neighbors);
            for (int j = 0; j < num_neighbors; ++j)
            {
                const int neighbor_index = neighbors_list[i][j];

                // Calculate the artificial tensile pressure correction
                const Scalar kernel_val = calcKernel(p.p - particles[neighbor_index].p, radius);
                const Scalar ratio      = kernel_val / corr_w;
                const Scalar corr_coeff = -corr_k * std::pow(ratio, corr_n);

                const Scalar coeff = particles[neighbor_index].m * (lambda[i] + lambda[neighbor_index] + corr_coeff);

                buffer.col(j) = coeff * calcGradKernel(p.p - particles[neighbor_index].p, radius);
            }
            const Vec3 sum = buffer.rowwise().sum();

            // Calculate delta p of this particle
            delta_p.col(i) = (1.0 / p.m) * (1.0 / rest_density) * sum;
        };
        parallelutil::parallel_for(num_particles, calc_delta_p);

        // Apply delta p (note: Jacobi style)
        const auto apply_delta_p = [&](const int i) { particles[i].p += delta_p.col(i); };
        parallelutil::parallel_for(num_particles, apply_delta_p);

        for (int i = 0; i < num_particles; ++i)
        {
            auto& p = particles[i];

            // Detect and handle environmental collisions (in a very naive way)
            p.p = p.p.cwiseMax(Vec3(-1.0, 0.0, -1.0));
            p.p = p.p.cwiseMin(Vec3(+1.0, 8.0, +1.0));
        }
    }

    // Update positions and velocities
    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].v = damping * (particles[i].p - particles[i].x) / dt;
        particles[i].x = particles[i].p;
    }
}

void test_kernel()
{
    constexpr auto   calcFunc     = calcPoly6Kernel;
    constexpr auto   calcGradFunc = calcGradPoly6Kernel;
    constexpr Scalar epsilon      = 1e-08;
    constexpr Scalar radius       = 0.2;

    const Vec3 x = radius * Vec3::Random();

    const Scalar val = calcFunc(x, radius);

    Vec3 numerical_grad;
    for (int i = 0; i < 3; ++i)
    {
        Vec3 eps = Vec3::Zero();
        eps[i]   = epsilon;

        numerical_grad[i] = (calcFunc(x + eps, radius) - val) / epsilon;
    }

    const Vec3 grad = calcGradFunc(x, radius);

    assert((grad - numerical_grad).norm() < 1e-04);
}

int main()
{
    std::vector<Particle> particles;

    constexpr Scalar dt         = 1.0 / 60.0;
    constexpr int    num_frames = 240;

    constexpr int x_size        = 60;
    constexpr int y_size        = 60;
    constexpr int z_size        = 30;
    constexpr int num_particles = x_size * y_size * z_size;

    std::cout << "#particles: " << num_particles << std::endl;

    // Generate and initialize particles
    particles.resize(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].i = i;
        particles[i].m = 3000.0 / static_cast<Scalar>(num_particles);
        particles[i].x = Vec3(0.5, 1.5, 0.5).cwiseProduct(Vec3::Random()) + Vec3(-0.5, 2.5, 0.0);
        particles[i].v = Vec3::Zero();
    }

    // Instantiate an alembic manager and submit the initial status
    ParticlesAlembicManager alembic_manager("./test.abc", dt, "Fluid", &particles);
    alembic_manager.submitCurrentStatus();

    // Simulate particles
    constexpr int    num_substeps = 5;
    constexpr Scalar sub_dt       = dt / static_cast<Scalar>(num_substeps);
    for (int t = 0; t < num_frames; ++t)
    {
        // Instantiate timer
        const auto timer = timer::Timer("Frame #" + std::to_string(t));

        // Step the simulation time
        for (int k = 0; k < num_substeps; ++k)
        {
            step(sub_dt, particles);

            // Note: A dirty solution for managing a bad initial state
            if (t < 5)
            {
                for (auto& p : particles)
                {
                    p.v = Vec3::Zero();
                }
            }
        }

        // Write the current status
        alembic_manager.submitCurrentStatus();
    }

    return 0;
}
