#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN
class AreaLight : public Emitter{
public:
	AreaLight(const PropertyList &props){
		m_radiance = props.getColor("radiance", Color3f(10.0f));
	}

	Color3f sample(Sampler *sampler, EmitterSampleRecord &record) const{
		if(record.sampleMesh == nullptr)
			throw NoriException("There is no mesh attached !");
		record.sampleMesh->sampleMesh(sampler, record.samplePoint, record.sampleNormal);
		// std::cout << record.samplePoint.toString() << std::endl;
		auto out_dir = record.surfacePoint - record.samplePoint;
		record.wo = out_dir.normalized();
		record.distance2 = out_dir.squaredNorm();
		record.cosTheta = record.wo.dot(record.sampleNormal);
		record.pdf = pdf(record);
		
		if(record.pdf == 0.0f)
			return 0.0f;
		return eval(record) / record.pdf;
	}

	Color3f eval(const EmitterSampleRecord &record) const{
		if(record.cosTheta <= 0.0f)
			return Color3f(0.0f);
		return m_radiance;
	}

	float pdf(const EmitterSampleRecord &record) const{
		if(record.sampleMesh == nullptr)
			throw NoriException("There is no mesh attached !");	
		float sample_pdf = record.sampleMesh->getPdf() * record.distance2 / std::abs(record.cosTheta);
		if(isnan(sample_pdf) || std::abs(sample_pdf) == INFINITY)	
			return  0.0f;
		return sample_pdf;
	}

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "AreaEmitter[\n"
            "  radiance = %s\n"
            "]", m_radiance.toString());
    }
private:
	Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END