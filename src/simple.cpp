#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator{
public:
	SimpleIntegrator(const PropertyList &props) 
	: light_position(props.getPoint("position", Point3f(0, 20.0f, 0))), 
	  light_energy(props.getColor("energy", Color3f(1.0e4, 1.0e4, 1.0e4))){
	}
	
	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const{
		Intersection its;
    	if(!scene->rayIntersect(ray, its))
           return Color3f(0.0f);
        
        Point3f shadePos = its.p;
        auto dir = light_position - shadePos;
        Ray3f rayDir = Ray3f(shadePos, dir);
        if(scene->rayIntersect(rayDir))
        	return Color3f(0.0f);
       	return Color3f(1.0f);
/*        Normal3f shadeNormal = its.shFrame.n.normalized();
        float distance = dir.norm();
        float cosTheta = shadeNormal.dot(dir.normalized());
        
        
        Color3f nom = light_energy * INV_PI * INV_PI * std::max(0.0f, cosTheta);
        float denom = 4.0f * distance * distance;

		return nom / denom;
*/		
	}

	std::string toString() const {
        return "SimpleIntegrator[]";
    }
private:
	Point3f light_position;
	Color3f light_energy;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END