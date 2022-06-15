#ifndef PBF_ALEMBIC_MANAGER_BASE_HPP
#define PBF_ALEMBIC_MANAGER_BASE_HPP

#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcGeom/All.h>
#include <cstdint>
#include <string>

class AlembicManagerBase {
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

        if (m_is_first) {
            submitCurrentStatusFirstTime();

            m_is_first = false;

            return;
        }

        const float* position_array = getAlembicPositionArray();

        const V3fArraySample        position_array_sample(reinterpret_cast<const V3f*>(position_array), getNumVerts());
        const OPointsSchema::Sample sample(position_array_sample);

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

#endif
