#ifndef PBF_NEIGHBOR_SEARCH_ENGINE_HPP
#define PBF_NEIGHBOR_SEARCH_ENGINE_HPP

#include "particle.hpp"
#include "types.hpp"
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

/// \brief The base class for neighborhood search engines
class NeighborSearchEngineBase {
public:
    NeighborSearchEngineBase(const Scalar radius, const std::vector<Particle>& particles)
        : m_radius(radius), m_particles(particles)
    {
    }

    virtual void searchNeighbors() = 0;

    const std::vector<int>& retrieveNeighbors(const int particle_index) const
    {
        return m_neighbors_list[particle_index];
    };

    int getNumParticles() const { return m_neighbors_list.size(); }

protected:
    std::vector<std::vector<int>> m_neighbors_list;

    const Scalar                 m_radius;
    const std::vector<Particle>& m_particles;
};

/// \brief Neighborhood search engine based on a simple uniform grid construction
class UniformGridNeighborSearchEngine : public NeighborSearchEngineBase {
public:
    UniformGridNeighborSearchEngine(const Scalar radius, const std::vector<Particle>& particles);

    void searchNeighbors() override;

private:
    using GridIndex = std::tuple<int, int, int>;

    GridIndex calcGridIndex(const Scalar radius, const MatX& positions, const int target_particle_index);

    int convertGridIndexToArrayIndex(const GridIndex& index);

    std::unordered_map<int, std::vector<int>> constructGridCells(const Scalar radius, const MatX& positions);

    const int  k_n_x         = 128; // Sufficiently large value
    const int  k_n_y         = 128; // Sufficiently large value
    const int  k_n_z         = 128; // Sufficiently large value
    const Vec3 k_grid_center = Vec3(0.0, 1.0, 0.0);
};

#endif
