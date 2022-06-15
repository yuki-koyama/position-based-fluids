#ifndef PBF_PARTICLES_ALEMBIC_MANAGER_HPP
#define PBF_PARTICLES_ALEMBIC_MANAGER_HPP

#include "alembic-manager-base.hpp"
#include "particle.hpp"

class ParticlesAlembicManager : public AlembicManagerBase
{
public:
    ParticlesAlembicManager(const std::string&           output_file_path,
                            const double                 delta_time,
                            const std::string&           object_name,
                            const std::vector<Particle>* particles_ptr)
        : AlembicManagerBase(output_file_path, delta_time, object_name), m_particles_ptr(particles_ptr)
    {
        m_position_buffer.resize(3, getNumVerts());
    }

private:
    const float* getAlembicPositionArray() override
    {
        for (std::size_t index = 0; index < getNumVerts(); ++index)
        {
            m_position_buffer.col(index) = (*m_particles_ptr)[index].x.cast<float>();
        }
        return m_position_buffer.data();
    }

    std::int32_t getNumVerts() const override { return m_particles_ptr->size(); }

    void submitCurrentStatusFirstTime() override
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const std::size_t               num_counts = getNumVerts();
        const std::vector<std::int32_t> counts(num_counts, 1);

        std::vector<uint64_t> index_buffer(num_counts);
        for (std::size_t elem_index = 0; elem_index < num_counts; ++elem_index)
        {
            index_buffer[elem_index] = elem_index;
        }

        const float* position_array = getAlembicPositionArray();

        const V3fArraySample        position_array_sample(reinterpret_cast<const V3f*>(position_array), getNumVerts());
        const UInt64ArraySample     index_array_sample(index_buffer.data(), getNumVerts());
        const OPointsSchema::Sample sample(position_array_sample, index_array_sample);

        m_points.getSchema().set(sample);
    }

    const std::vector<Particle>* m_particles_ptr;

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> m_position_buffer;
};

#endif
