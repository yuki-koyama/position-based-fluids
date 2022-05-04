#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcGeom/All.h>
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

class AlembicManagerBase
{
public:
    AlembicManagerBase(const std::string& output_file_path, const double delta_time, const std::string& object_name)
        : m_archive(Alembic::AbcCoreOgawa::WriteArchive(), output_file_path.c_str())
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const TimeSampling  time_sampling(delta_time, 0);
        const std::uint32_t time_sampling_index = m_archive.addTimeSampling(time_sampling);

        m_points = OPoints(OObject(m_archive, kTop), object_name.c_str());
        m_points.getSchema().setTimeSampling(time_sampling_index);
    }

    void submitCurrentStatus()
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        if (m_is_first)
        {
            submitCurrentStatusFirstTime();

            m_is_first = false;

            return;
        }

        const OPointsSchema::Sample sample(V3fArraySample((const V3f*) getAlembicPositionArray(), getNumVerts()));

        m_points.getSchema().set(sample);
    }

protected:
    virtual const float* getAlembicPositionArray()      = 0;
    virtual std::int32_t getNumVerts() const            = 0;
    virtual void         submitCurrentStatusFirstTime() = 0;

    bool m_is_first = true;

    Alembic::Abc::OArchive    m_archive;
    Alembic::AbcGeom::OPoints m_points;
};

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

        const OPointsSchema::Sample sample(V3fArraySample((const V3f*) getAlembicPositionArray(), getNumVerts()),
                                           UInt64ArraySample(index_buffer.data(), getNumVerts()));

        m_points.getSchema().set(sample);
    }

    const std::vector<Particle>* m_particles_ptr;

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> m_position_buffer;
};

int main()
{
    std::vector<Particle> particles;

    constexpr Scalar dt = 1.0 / 60.0;

    constexpr int x_size = 6;
    constexpr int y_size = 6;
    constexpr int z_size = 6;

    // Generate and initialize particles
    constexpr int num_particles = x_size * y_size * z_size;
    particles.resize(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        particles[i].i = i;
        particles[i].m = 1.0 / static_cast<Scalar>(num_particles);
        particles[i].x = Vec3::Random();
        particles[i].v = Vec3::Zero();
    }

    ParticlesAlembicManager alembic_manager("./test.abc", dt, "fluid", &particles);
    alembic_manager.submitCurrentStatus();

    // Print particle status
    for (auto& particle : particles)
    {
        std::cout << particle.x.transpose() << std::endl;
    }

    return 0;
}
