#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AOIntegrator : public Integrator{
public:
	AOIntegrator(const PropertyList &props) {
	}
	
	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const{
		Intersection its;
    	if(!scene->rayIntersect(ray, its))
           return Color3f(0.0f);
        
        Point3f shadePos = its.p;
        Normal3f shadeNormal = its.shFrame.n.normalized();
        
        auto localDir = Warp::squareToCosineHemisphere(sampler->next2D());
        float pdf = Warp::squareToCosineHemispherePdf(localDir);
        auto worldDir = its.shFrame.toWorld(localDir.normalized()).normalized();
        
        Ray3f rayDir = Ray3f(shadePos, worldDir);
        if(scene->rayIntersect(rayDir))
        	return Color3f(0.0f);

        float cosTheta = shadeNormal.dot(worldDir);
		return Color3f(cosTheta * INV_PI / pdf);		
	}

	std::string toString() const {
        return "AOIntegrator[]";
    }
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END