#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class EmitterSamplingIntegrator : public Integrator{
public:
	EmitterSamplingIntegrator(const PropertyList &props) {
	}
	
	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const{
        Intersection its;
        Color3f Lo(0.0f);
        Color3f throughput(1.0f);
        Ray3f ray_out(ray);
        float is_delta = 1.0f;      // handle double-counting

        for(unsigned i=0; ;++i){
            if(!scene->rayIntersect(ray_out, its))
               break;

            // russian roulette            
            if(i >= 3){
                float rr = std::min(throughput.maxCoeff(), 0.99f);
                if(rr < sampler->next1D())
                    break;
                throughput /= rr;
            }

            // emitter
            if(its.mesh->isEmitter()){
                auto normal = its.shFrame.n; 
                EmitterSampleRecord recE(its.mesh, -ray_out.d, normal);
                Lo += throughput * its.mesh->getEmitter()->eval(recE) * is_delta;
            }

            if(its.mesh->getBSDF()->isDiffuse()){
                // direct
                for(auto light : scene->getLights()){
                    EmitterSampleRecord rec(light, its.p);
                    Color3f L = light->getEmitter()->sample(sampler, rec);
                    Ray3f shadowRay = Ray3f(rec.surfacePoint, -rec.wo, Epsilon, std::sqrt(rec.distance2) - Epsilon);
                    
                    BSDFQueryRecord bRec(its.toLocal(-ray_out.d), its.toLocal(-rec.wo), ESolidAngle);
                    Color3f brdf = its.mesh->getBSDF()->eval(bRec);

                    if(!scene->rayIntersect(shadowRay)){    
                        Lo += L * brdf * throughput * Frame::cosTheta(bRec.wo);              
                    }   
                }
                is_delta = 0.0f;                
            }else{
                is_delta = 1.0f;
            }

            // indirect
            BSDFQueryRecord inRec(its.toLocal(-ray_out.d));
            throughput *= its.mesh->getBSDF()->sample(inRec, sampler->next2D());
            ray_out = Ray3f(its.p, its.toWorld(inRec.wo));                 
        }
        return Lo;
	}

	std::string toString() const {
        return "EmitterSamplingIntegrator[]";
    }
};

NORI_REGISTER_CLASS(EmitterSamplingIntegrator, "path_ems");
NORI_NAMESPACE_END