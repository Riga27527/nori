/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    auto tent = [](float Xi){
        Xi = 2.0f * Xi;
        if(Xi < 1)
            return std::sqrt(Xi) - 1.0f;
        else
            return 1.0f - std::sqrt(2.0f-Xi);
    };
    return Point2f(tent(sample.x()), tent(sample.y()));
    // throw NoriException("Warp::squareToTent() is not yet implemented!");
}

float Warp::squareToTentPdf(const Point2f &p) {
    auto pdf = [](float pi){
        return pi >= -1.0f && pi <= 1.0f ? 1.0f - std::abs(pi) : 0.0f;    
    };
    return pdf(p.x()) * pdf(p.y());
    // throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = std::sqrt(sample.x()), theta = 2.0f * M_PI * sample.y();
    return Point2f(r * std::cos(theta), r * std::sin(theta));
    // throw NoriException("Warp::squareToUniformDisk() is not yet implemented!");
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return (p.x() * p.x() + p.y() * p.y()) <= 1.0f ? INV_PI : 0.0f;
    // throw NoriException("Warp::squareToUniformDiskPdf() is not yet implemented!");
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float theta = std::acos(1.0f - 2.0f * sample.x()), phi = 2.0f * M_PI * sample.y();
    return Vector3f(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
    // throw NoriException("Warp::squareToUniformSphere() is not yet implemented!");
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return INV_FOURPI;
    // throw NoriException("Warp::squareToUniformSpherePdf() is not yet implemented!");
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float theta = std::acos(1.0f - sample.x()), phi = 2.0f * M_PI * sample.y();
    return Vector3f(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)); 
    // throw NoriException("Warp::squareToUniformHemisphere() is not yet implemented!");
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return v.z() < 0.0f ? 0.0f : INV_TWOPI;
    // throw NoriException("Warp::squareToUniformHemispherePdf() is not yet implemented!");
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    float theta = std::acos(1.0f - 2.0f * sample.x()) / 2.0f, phi = 2.0f * M_PI * sample.y();
    return Vector3f(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)); 
    // throw NoriException("Warp::squareToCosineHemisphere() is not yet implemented!");
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return v.z() < 0.0f ? 0.0f : v.z() * INV_PI;
    // throw NoriException("Warp::squareToCosineHemispherePdf() is not yet implemented!");
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float phi = 2.0f * M_PI * sample.x();
    float alpha2 = alpha * alpha;
    float tmp = std::sqrt(-alpha2 * std::log(1.0f - sample.y()));
    float theta = std::atan(tmp);
    return Vector3f(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)); 
    // throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if(m.z() <= 0)
        return 0.0f;
    float alpha2 = alpha * alpha;
    float cosTheta = m.z(), cosTheta2 = m.z() * m.z(), sinTheta2 = m.x() * m.x() + m.y() * m.y();
    float tanTheta2 = sinTheta2 / cosTheta2;
    float nom =  std::exp(-tanTheta2 / alpha2);
    float denom = alpha2 * cosTheta2 * cosTheta;
    return INV_PI * (nom / denom);
    // throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
}

NORI_NAMESPACE_END
