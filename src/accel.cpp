/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/
 
#include <nori/accel.h>
#include <Eigen/Geometry>
#include <nori/timer.h>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

OctNode* Accel::buildOctree(const BoundingBox3f& box, const std::vector<uint32_t>& tri_ind, unsigned depth){
    if(tri_ind.size() <= 10 || depth >= octree_maxDepth){
        nodes_num += 1;
        leaf_num += 1;
        tris_per_leaf += tri_ind.size();
        return new OctNode(box, tri_ind);
    }

    // divide bounding box
    std::vector<BoundingBox3f> bbox_list(8);
    // auto center = box.getCenter();
    // for(size_t i=0; i<8; ++i){
    //     auto corner = box.getCorner(i);
    //     decltype(corner) st, ed;
    //     for(size_t j=0; j<3; ++j){
    //         st[j] = std::min(center[j], corner[j]); 
    //         ed[j] = std::max(center[j], corner[j]);    
    //     }
    //     bbox_list[i] = BoundingBox3f(st, ed);
    // }
    Point3f extents = box.getExtents() / 2.0f;
    float dx = extents[0], dy = extents[1], dz = extents[2];
    int n = 0;
    for(float x : {0.f, dx})
        for(float y : {0.f, dy})
            for(float z : {0.f, dz}){
                Point3f ori = box.min + Point3f(x, y, z);
                bbox_list[n++] = BoundingBox3f(ori, ori + extents);
            }

    // judge overlap
    std::vector<uint32_t> tri_list[8];
    for(auto ind : tri_ind){
        for(size_t i=0; i<8; ++i){
            if(bbox_list[i].overlaps(m_mesh->getBoundingBox(ind)))
                tri_list[i].push_back(ind);
        }
    }

    nodes_num += 1;
    OctNode *node = new OctNode(box);
    for(size_t i=0; i<8; ++i)
        if(tri_list[i].size() > 0)
            node->children.push_back(buildOctree(bbox_list[i], tri_list[i], depth + 1));
    
    return node;
}

void Accel::build() {
    /* Nothing to do here for now */
    uint32_t len_tri = m_mesh->getTriangleCount();

    std::vector<uint32_t> vec_tri(len_tri);
    for(uint32_t i=0; i<len_tri; ++i)
        vec_tri[i] = i;
    
    Timer timer_build;
    m_root = buildOctree(m_bbox, vec_tri, 0);
    std::cout << "Octree done. (took " << timer_build.elapsedString() << ")" << std::endl;
    tris_per_leaf /= leaf_num;
    printOctreeInfo();
}

bool Accel::traverseOctree(Ray3f &ray, OctNode* node, Intersection &its, uint32_t &f, bool shadowRay) const{
    bool hit = false;
    if(node->children.size()==0){
        for(auto idx : node->triangle_List){
            float u, v, t;
            if (m_mesh->rayIntersect(idx, ray, u, v, t) && t < ray.maxt) {
                if (shadowRay)
                    return true;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                f = idx;
                hit = true;
            }            
        }
        return hit;
    }

    for(OctNode* child : node->children){
        if(child->oct_bbox.rayIntersect(ray) && traverseOctree(ray, child, its, f, shadowRay)){
            hit = true;
            if(shadowRay)
                return true;
        }
    }
    return hit;
    
    // std::vector<std::pair<OctNode*, float> > dis2ray;
    // for(OctNode* child : node->children){
    //     float near, far;
    //     if(child->oct_bbox.rayIntersect(ray, near, far));
    //         dis2ray.push_back(std::make_pair(child, near));
    // }

    // std::sort(dis2ray.begin(), dis2ray.end(), [](const auto &n1, const auto &n2){
    //     return n1.second < n2.second;
    // });

    // for(auto child : dis2ray){
    //     if(traverseOctree(ray, child.first, its, f, shadowRay))
    //         return true;
    // }
    // return false;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)
    
    if(m_root->oct_bbox.rayIntersect(ray))
        foundIntersection = traverseOctree(ray, m_root, its, f, shadowRay);

    if(shadowRay)
        return foundIntersection;
    /* Brute force search through all triangles */
/*    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        float u, v, t;
        if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
             An intersection was found! Can terminate
               immediately if this is a shadow ray query 
            if (shadowRay)
                return true;
            ray.maxt = its.t = t;
            its.uv = Point2f(u, v);
            its.mesh = m_mesh;
            f = idx;
            foundIntersection = true;
        }
    }*/


    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

