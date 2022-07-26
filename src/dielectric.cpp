/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        // throw NoriException("Unimplemented!");
        float cosTheta = Frame::cosTheta(bRec.wi);        
        float fresnel_term = fresnel(cosTheta, m_extIOR, m_intIOR);
        Normal3f localNormal = Normal3f(0.0f, 0.0f, 1.0f);
        bRec.measure = EDiscrete;
        bRec.eta = m_extIOR / m_intIOR;

        if(fresnel_term > sample.x()){
            bRec.wo = reflect(bRec.wi, localNormal);      
        }else{
            if(cosTheta < 0.0f){
                bRec.eta = 1.0f / bRec.eta;
                localNormal = -localNormal;
            }          
            bRec.wo = refract(-bRec.wi, localNormal, bRec.eta);
        }
        return Color3f(1.0f);

        // float cosTheta = Frame::cosTheta(bRec.wi);
        // float extIOR = m_extIOR, intIOR = m_intIOR;
        // Normal3f localNormal = Normal3f(0.0f, 0.0f, 1.0f);
        // if(cosTheta < 0.0f){
        //     std::swap(extIOR, intIOR);
        //     localNormal = -localNormal;
        // }
        // bRec.measure = EDiscrete;
        // bRec.eta = extIOR;
        // float refract_index = m_extIOR / m_intIOR;

        // float sinTheta = Frame::sinTheta(bRec.wi);
        // bool cannot_refract = refract_index * sinTheta > 1.0f;
        // float fresnel_term = fresnel(cosTheta, extIOR, intIOR);
        
        // if(cannot_refract || fresnel_term > sample.x()){
        //     bRec.wo = reflect(bRec.wi, localNormal);
        // }else{
        //     bRec.wo = refract(bRec.wi, localNormal, refract_index);
        // }

        return Color3f(1.0f);
    }

    Vector3f reflect(const Vector3f &wi, const Normal3f &normal) const{
        return (2.0f * wi.dot(normal) * normal - wi).normalized();
    }

    Vector3f refract(const Vector3f &wi_in, const Normal3f &normal, float refract_index) const{
        float cosTheta = -wi_in.dot(normal);
        float discriminant = 1.0f - refract_index * refract_index * (1.0f - cosTheta * cosTheta);
        if(discriminant > 0.0f){
            return refract_index * (wi_in + cosTheta * normal) - normal * std::sqrt(discriminant);
        }else
            return Vector3f(0.0f);
        // float cosTheta = wi_in.dot(normal);
        // Vector3f wo_perp = refract_index * (-wi_in + cosTheta * normal);
        // Vector3f wo_parallel = -std::sqrt(std::abs(1.0f - wo_perp.squaredNorm())) * normal;
        // return (wo_perp + wo_parallel).normalized();
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
