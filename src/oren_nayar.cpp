/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Diffuse / OrenNayar BRDF model
 */
class OrenNayar : public BSDF {
public:
    OrenNayar(const PropertyList &propList) {
        m_type = BSDFType::BSDF_DIFFUSE;
        std::string tex_name = propList.getString("albedo_tex", "none");
        if(tex_name == "none")
            m_albedo = propList.getColor("albedo", Color3f(0.5f));
        else{
            m_hasTexture = true;
            m_texture = Texture(tex_name);
        }
        float sigma = degToRad(propList.getFloat("sigma", 60.f));
        float sigma2 = sigma * sigma;
        A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
        B = 0.45f * sigma2 / (sigma2 + 0.09);
    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        float sinThetaI = Frame::sinTheta(bRec.wi);
        float sinThetaO = Frame::sinTheta(bRec.wo);
        float maxCos = 0;
        if (sinThetaI > 1e-4 && sinThetaO > 1e-4) {
            float sinPhiI = Frame::sinPhi(bRec.wi), cosPhiI = Frame::cosPhi(bRec.wi);
            float sinPhiO = Frame::sinPhi(bRec.wo), cosPhiO = Frame::cosPhi(bRec.wo);
            float dCos = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
            maxCos = std::max(0.f, dCos);
        }

        // Compute sine and tangent terms of Oren-Nayar model
        float sinAlpha, tanBeta;
        float absCosThetaI = std::abs(Frame::cosTheta(bRec.wi)), absCosThetaO = std::abs(Frame::cosTheta(bRec.wi));
        if (absCosThetaO > absCosThetaI) {
            sinAlpha = sinThetaI;
            tanBeta = sinThetaO / absCosThetaO;
        } else {
            sinAlpha = sinThetaO;
            tanBeta = sinThetaI / absCosThetaI;
        }
        float oren = (A + B * maxCos * sinAlpha * tanBeta);

        /* The BRDF is simply the albedo / pi */
        if(m_hasTexture && bRec.uv_valid)
            return m_texture.sample2D(bRec.uv) * INV_PI * oren;
        else
            return m_albedo * INV_PI * oren;
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;


        /* Importance sampling density wrt. solid angles:
           cos(theta) / pi.

           Note that the directions in 'bRec' are in local coordinates,
           so Frame::cosTheta() actually just returns the 'z' component.
        */
        return INV_PI * Frame::cosTheta(bRec.wo);
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */
        bRec.wo = Warp::squareToCosineHemisphere(sample);

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        return eval(bRec) / pdf(bRec) * Frame::cosTheta(bRec.wo);
    }

    bool isDiffuse() const {
        return true;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "OrenNayar[\n"
            "  albedo = %s\n"
            "]", m_albedo.toString());
    }

    EClassType getClassType() const { return EBSDF; }
private:
    Color3f m_albedo;
    float A, B;
    Texture m_texture;
    bool    m_hasTexture = false;
};

NORI_REGISTER_CLASS(OrenNayar, "OrenNayar");
NORI_NAMESPACE_END
