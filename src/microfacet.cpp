/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        m_metalic = propList.getFloat("metalic", 0.1f);
        
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material */
        std::string albedo_tex_filename = propList.getString("albedo_tex", "none");
        if(albedo_tex_filename == "none"){
            Color3f albedo = propList.getColor("albedo", Color3f(0.5f));
            m_kd =  albedo * (1.0f - m_metalic);
            m_ks = lerp(m_metalic, Color3f(0.04f), albedo);
        }
        else{
            m_albedo_tex = Texture(albedo_tex_filename);
            m_hasTexture = true;
        }

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        
    }

    float GeometryTerm(const Vector3f wi, const Vector3f wo, const Vector3f wh, float alpha) const{
        auto G1 = [=](const Vector3f wv, const Vector3f wh){
            float c = wv.dot(wh) / Frame::cosTheta(wv);
            if(c <= 0.0f)
                return 0.0f;

            float b = 1.0f / (alpha * Frame::tanTheta(wv)), b2 = b * b;
            if(b >= 1.6f)
                return 1.0f;
            float numerator = 3.535f * b + 2.181f * b2;
            float denominator = 1.0f + 2.276f * b + 2.577f * b2;
            return numerator / denominator;
        };
        return G1(wi, wh) * G1(wo, wh);
    }

    float Beckmann_NDF(const Vector3f wh, float alpha) const{
        float cosTheta = Frame::cosTheta(wh), tanTheta = Frame::tanTheta(wh);
        float cosTheta3 = cosTheta * cosTheta * cosTheta, tanTheta2 = tanTheta * tanTheta, alpha2 = alpha * alpha;
        float azimuthal = INV_PI;
        float longitude = std::exp(-tanTheta2 / alpha2) / (alpha2 * cosTheta3);
        return azimuthal * longitude;
    }

    // cosTheta is view of normal
    Color3f fresnelSchlick(Color3f F0, float cosTheta) const{
        return F0 + (1.0f - F0) * std::pow(1.0f - cosTheta, 5.0f);
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
    	// throw NoriException("MicrofacetBRDF::eval(): not implemented!");
        if(Frame::cosTheta(bRec.wi) <= 0.0f || Frame::cosTheta(bRec.wo) <= 0.0f )
            return Color3f(0.0f);

        Color3f kd = m_kd, ks = m_ks;
        if(m_hasTexture && bRec.uv_valid){
            Color3f albedo = m_albedo_tex.sample2D(bRec.uv);
            kd = albedo * (1.0f - m_metalic);
            ks = lerp(m_metalic, Color3f(0.04f), albedo);
        }

        Color3f diff = kd * INV_PI;

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        // float VoH = bRec.wi.dot(wh);
        auto D = Beckmann_NDF(wh, m_alpha), G = GeometryTerm(bRec.wi, bRec.wo, wh, m_alpha);
        auto F = fresnelSchlick(ks, Frame::cosTheta(bRec.wi));
        float cosTheta_i = Frame::cosTheta(bRec.wi), cosTheta_o = Frame::cosTheta(bRec.wo);
        Color3f spec = (D * F * G) / (4.f * cosTheta_i * cosTheta_o);

        return diff + spec;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
    	// throw NoriException("MicrofacetBRDF::pdf(): not implemented!");
        if(Frame::cosTheta(bRec.wi) <= 0.0f || Frame::cosTheta(bRec.wo) <= 0.0f)
            return 0.0f;

        Color3f kd = m_kd, ks = m_ks;
        if(m_hasTexture && bRec.uv_valid){
            Color3f albedo = m_albedo_tex.sample2D(bRec.uv);
            kd = albedo * (1.0f - m_metalic);
            ks = lerp(m_metalic, Color3f(0.04f), albedo);
        }
        float ks_max = ks.maxCoeff();

        auto wh = (bRec.wi + bRec.wo).normalized();
        float Jacobian = 0.25 / wh.dot(bRec.wo);
        float diff_pdf = (1.0f - ks_max) * Warp::squareToCosineHemispherePdf(bRec.wo);
        float spec_pdf = ks_max * Warp::squareToBeckmannPdf(wh, m_alpha) * Jacobian;
        return diff_pdf + spec_pdf;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
    	// throw NoriException("MicrofacetBRDF::sample(): not implemented!");
        if(Frame::cosTheta(bRec.wi) <= 0.0f)
            return Color3f(0.0f);

        Color3f kd = m_kd, ks = m_ks;
        if(m_hasTexture && bRec.uv_valid){
            Color3f albedo = m_albedo_tex.sample2D(bRec.uv);
            kd = albedo * (1.0f - m_metalic);
            ks = lerp(m_metalic, Color3f(0.04f), albedo);
        }
        float ks_max = ks.maxCoeff();

        if(_sample.x() > ks_max){
            // diffuse
            Point2f sample = Point2f((_sample.x() - ks_max) / (1.0f - ks_max), _sample.y());
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }else{
            // specular
            Point2f sample = Point2f((_sample.x() / ks_max), _sample.y());
            auto wh = Warp::squareToBeckmann(sample, m_alpha);
            bRec.wo = (2.f * wh.dot(bRec.wi) * wh - bRec.wi).normalized();
        }
        if(Frame::cosTheta(bRec.wo) <= 0.0f)
            return Color3f(0.0f);
        bRec.eta = 1.0f;
        bRec.measure = ESolidAngle;
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }


    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  metalic = %f,\n"
            // "  kd = %s,\n"
            // "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_metalic
            // m_kd.toString(),
            // m_ks
        );
    }

private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_metalic;
    Color3f m_kd;
    Color3f m_ks;
    Texture m_albedo_tex;
    bool m_hasTexture = false;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
