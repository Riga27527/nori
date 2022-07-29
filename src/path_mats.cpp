#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class MatSamplingIntegrator : public Integrator{
public:
	MatSamplingIntegrator(const PropertyList &props) {
	}
	
	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const{
        Intersection its;
        Color3f Lo(0.0f);
        Color3f throughput(1.0f);
        Ray3f ray_out(ray);

        for(unsigned i=0; ;++i){
            if(!scene->rayIntersect(ray_out, its))
               break;

            // russian roulette            
            if(i > 3){
                float rr = std::min(throughput.maxCoeff(), 0.99f);
                if(rr < sampler->next1D())
                    break;
                throughput /= rr;
            }

            // emitter
            if(its.mesh->isEmitter()){
                auto normal = its.shFrame.n; 
                EmitterSampleRecord rec(its.mesh, -ray_out.d, normal);
                Lo += throughput * its.mesh->getEmitter()->eval(rec);
            }

            // BSDF
            BSDFQueryRecord bRec(its.toLocal(-ray_out.d));
            throughput *= its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            ray_out = Ray3f(its.p, its.toWorld(bRec.wo));            
        }
        return Lo;
	}

	std::string toString() const {
        return "MatSamplingIntegrator[]";
    }
};

NORI_REGISTER_CLASS(MatSamplingIntegrator, "path_mats");
NORI_NAMESPACE_END