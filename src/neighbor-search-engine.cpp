#include "neighbor-search-engine.hpp"

UniformGridNeighborSearchEngine::UniformGridNeighborSearchEngine(const Scalar                 radius,
                                                                 const std::vector<Particle>& particles)
    : NeighborSearchEngineBase(radius, particles)
{
}

void UniformGridNeighborSearchEngine::searchNeighbors()
{
    const int    num_particles  = m_particles.size();
    const Scalar radius_squared = m_radius * m_radius;

    MatX positions(3, num_particles);
    for (int i = 0; i < num_particles; ++i) {
        positions.col(i) = m_particles[i].p;
    }

    auto grid_cells = constructGridCells(m_radius, positions);

    m_neighbors_list = std::vector<std::vector<int>>(num_particles);

    for (int i = 0; i < num_particles; ++i) {
        // Determine the cell that the target particle belongs to
        const auto target_cell_index = calcGridIndex(m_radius, positions, i);

        // Visit the 26 neighbor cells and the cell itself (27 in total)
        for (int x : {-1, 0, 1}) {
            for (int y : {-1, 0, 1}) {
                for (int z : {-1, 0, 1}) {
                    const int i_x = std::get<0>(target_cell_index) + x;
                    const int i_y = std::get<1>(target_cell_index) + y;
                    const int i_z = std::get<2>(target_cell_index) + z;

                    const int cell_index = convertGridIndexToArrayIndex({i_x, i_y, i_z});

                    const auto& list = grid_cells[cell_index];

                    // Register particles as neighbors if they are within the range
                    for (const int index : list) {
                        const Scalar squared_dist = (m_particles[i].p - m_particles[index].p).squaredNorm();

                        if (squared_dist < radius_squared) {
                            m_neighbors_list[i].push_back(index);
                        }
                    }
                }
            }
        }
    }
}

UniformGridNeighborSearchEngine::GridIndex
UniformGridNeighborSearchEngine::calcGridIndex(const Scalar radius,
                                               const MatX&  positions,
                                               const int    target_particle_index)
{
    const auto pos            = positions.col(target_particle_index);
    const Vec3 grid_coord_pos = (pos - k_grid_center) * (1.0 / radius) + 0.5 * Vec3(k_n_x, k_n_y, k_n_z);

    const int i_x = std::floor(grid_coord_pos[0]);
    const int i_y = std::floor(grid_coord_pos[1]);
    const int i_z = std::floor(grid_coord_pos[2]);

    assert(i_x >= 0);
    assert(i_y >= 0);
    assert(i_z >= 0);
    assert(i_x < n_x);
    assert(i_y < n_y);
    assert(i_z < n_z);

    return GridIndex{i_x, i_y, i_z};
}

int UniformGridNeighborSearchEngine::convertGridIndexToArrayIndex(const GridIndex& index)
{
    return std::get<0>(index) + k_n_x * std::get<1>(index) + k_n_x * k_n_y * std::get<2>(index);
}

std::unordered_map<int, std::vector<int>> UniformGridNeighborSearchEngine::constructGridCells(const Scalar radius,
                                                                                              const MatX&  positions)
{
    std::unordered_map<int, std::vector<int>> grid_cells;

    for (int i = 0; i < positions.cols(); ++i) {
        const int array_index = convertGridIndexToArrayIndex(calcGridIndex(radius, positions, i));

        if (grid_cells.find(array_index) == grid_cells.end()) {
            grid_cells[array_index] = std::vector<int>{i};
        } else {
            grid_cells[array_index].push_back(i);
        }
    }

    return grid_cells;
}
