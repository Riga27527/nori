/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass of all emitters
 */

struct EmitterSampleRecord
{
    const Mesh* sampleMesh;

    Point3f samplePoint;

    Point3f surfacePoint;

    float pdf;

    Vector3f wo;

    Normal3f sampleNormal;

    float distance2;

    float cosTheta;

    EmitterSampleRecord(const Mesh *mesh) : sampleMesh(mesh){ }
    
    EmitterSampleRecord(const Mesh *mesh, const Point3f &surf_point) : sampleMesh(mesh), surfacePoint(surf_point){ }

    EmitterSampleRecord(const Mesh *mesh, const Vector3f &samp_wo, const Normal3f &samp_normal) 
        : sampleMesh(mesh), wo(samp_wo), sampleNormal(samp_normal), cosTheta(samp_wo.dot(samp_normal)){ }
};

class Emitter : public NoriObject {
public:

    virtual Color3f sample(Sampler *sampler, EmitterSampleRecord &record) const = 0;

    virtual Color3f eval(const EmitterSampleRecord &record)  const = 0;

    virtual float pdf(const EmitterSampleRecord &record) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }
};

NORI_NAMESPACE_END
