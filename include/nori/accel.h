/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/mesh.h>
#include <vector>

NORI_NAMESPACE_BEGIN

struct OctNode{
    BoundingBox3f oct_bbox;

    std::vector<uint32_t> triangle_List;

    std::vector<OctNode*> children;

    OctNode() = default;

    OctNode(BoundingBox3f bbox) : oct_bbox(bbox){}

    OctNode(BoundingBox3f bbox, std::vector<uint32_t> tri_List) : oct_bbox(bbox), triangle_List(tri_List){}
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }
 
    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

    OctNode* buildOctree(const BoundingBox3f& box, const std::vector<uint32_t>& tri_ind, unsigned depth);

    bool traverseOctree(Ray3f &ray, OctNode* node, Intersection &its, uint32_t &f, bool shadowRay) const;

    void printOctreeInfo() const{
        std::cout << "Nodes_num : " << nodes_num << "\nLeaf_num : " << leaf_num << \
        "\ntris_num_per_leaf : " << tris_per_leaf << std::endl;
    }

private:
    OctNode* buildOctree(const BoundingBox3f& box, const std::vector<uint32_t>& tri_ind);

    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
    OctNode      *m_root = nullptr;
    unsigned nodes_num = 0;
    unsigned leaf_num = 0;
    unsigned tris_per_leaf = 0;
    unsigned octree_maxDepth = 8;
};

NORI_NAMESPACE_END
