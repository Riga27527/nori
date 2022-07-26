#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator{
public:
	WhittedIntegrator(const PropertyList &props) {
	}
	
	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const{
		Intersection its;
    	if(!scene->rayIntersect(ray, its))
           return Color3f(0.0f);

        if(!its.mesh->getBSDF()->isDiffuse()){
            if(sampler->next1D() >= 0.95f)
                return Color3f(0.0f);
            BSDFQueryRecord bRec(its.toLocal(-ray.d));
            its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            const Ray3f ray_out = Ray3f(its.p, its.toWorld(bRec.wo));
            return Li(scene, sampler, ray_out) / 0.95f;
        }
        
        Color3f Le(0.0f);
        auto normal = its.shFrame.n;
        if(its.mesh->isEmitter()){
            EmitterSampleRecord rec(its.mesh, -ray.d, normal);
            Le += its.mesh->getEmitter()->eval(rec);
        }

        float weights = 1.0f / scene->getLights().size();
        auto Lr = Color3f(0.0f);
        for(auto light : scene->getLights()){
            // sample the area lights
            EmitterSampleRecord rec(light, its.p);
            Color3f throughput = light->getEmitter()->sample(sampler, rec);
            // must be attention of the shadowRay range !!!!
            Ray3f shadowRay = Ray3f(rec.surfacePoint, -rec.wo, Epsilon, std::sqrt(rec.distance2) - Epsilon);
            float cosTheta = normal.dot(-rec.wo);

            if(!scene->rayIntersect(shadowRay)){
                // eval brdf
                BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(-rec.wo), ESolidAngle);
                auto brdf = its.mesh->getBSDF()->eval(bRec);
                Lr += brdf * throughput * weights * std::max(cosTheta, 0.0f); 
                // Lr += Color3f(1.0f);
            }
        }

        return Le + Lr;
	}

	std::string toString() const {
        return "WhittedIntegrator[]";
    }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END