#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class MISSamplingIntegrator : public Integrator{
public:
	MISSamplingIntegrator(const PropertyList &props) {
	}
	
	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const{
        Intersection its;
        Color3f Lo(0.0f);
        Color3f throughput(1.0f);
        Ray3f ray_out(ray);
        float light_pdf = 1.0f, brdf_pdf = 1.0f;

        if(!scene->rayIntersect(ray_out, its))
            return Color3f(0.0f);

        if(its.mesh->isEmitter()){
            auto normal = its.shFrame.n; 
            EmitterSampleRecord recE(its.mesh, -ray_out.d, normal);
            Lo += its.mesh->getEmitter()->eval(recE);
        }

        for(unsigned i=0; ;++i){
            // russian roulette            
            if(i >= 3){
                float rr = std::min(throughput.maxCoeff(), 0.99f);
                if(rr <= sampler->next1D())
                    break;
                throughput /= rr;
            }

            // direct
            for(auto light : scene->getLights()){
                EmitterSampleRecord rec(light, its.p);
                Color3f L = light->getEmitter()->sample(sampler, rec);
                Ray3f shadowRay = Ray3f(rec.surfacePoint, -rec.wo, Epsilon, std::sqrt(rec.distance2) - Epsilon);
                
                BSDFQueryRecord bRec(its.toLocal(-ray_out.d), its.toLocal(-rec.wo), ESolidAngle);
                Color3f brdf = its.mesh->getBSDF()->eval(bRec);

                if(!scene->rayIntersect(shadowRay)){
                    light_pdf = light->getEmitter()->pdf(rec);
                    brdf_pdf = its.mesh->getBSDF()->pdf(bRec);
                    Color3f L_direct = L * brdf * throughput * std::abs(Frame::cosTheta(bRec.wo));   
                    if(light_pdf!=0.0f && brdf_pdf!=0.0f){
                        float mis_light = light_pdf / (light_pdf + brdf_pdf);                        
                        L_direct *= mis_light;   
                    }
                    Lo += L_direct;        
                }   
            }

            // indirect
            BSDFQueryRecord inRec(its.toLocal(-ray_out.d));
            throughput *= its.mesh->getBSDF()->sample(inRec, sampler->next2D());
            ray_out = Ray3f(its.p, its.toWorld(inRec.wo));
            if(scene->rayIntersect(ray_out, its)){
                if(its.mesh->isEmitter()){
                    brdf_pdf = its.mesh->getBSDF()->pdf(inRec);
                    EmitterSampleRecord recL(its.mesh, ray_out.o, its.p, its.shFrame.n);
                    light_pdf = its.mesh->getEmitter()->pdf(recL);
                    Color3f L = its.mesh->getEmitter()->eval(recL);
                    Color3f L_indirect = L * throughput;
                    if(light_pdf!=0.0f && brdf_pdf!=0.0f){
                        float mis_brdf = brdf_pdf / (brdf_pdf + light_pdf);    
                        L_indirect *= mis_brdf;   
                    }
                    Lo += L_indirect;                   
                }
            }else
                break;
        }
        return Lo;
	}

	std::string toString() const {
        return "MISSamplingIntegrator[]";
    }
};

NORI_REGISTER_CLASS(MISSamplingIntegrator, "path_mis");
NORI_NAMESPACE_END