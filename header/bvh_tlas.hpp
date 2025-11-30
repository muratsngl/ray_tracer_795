#pragma once
#include "types.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

inline TLAS tlas;
inline std::vector<BVH> bvh_list; 
inline std::vector<BLAS> blas_list;
// NEW: Global list to store and partition BLAS indices, mirroring bvh.prim_indices
inline std::vector<int> tlas_blas_indices; 



inline fl calculate_surface_area(const AABB& bbox) {
    fl dx = bbox.bbox_xmax - bbox.bbox_xmin;
    fl dy = bbox.bbox_ymax - bbox.bbox_ymin;
    fl dz = bbox.bbox_zmax - bbox.bbox_zmin;
    return 2.0f * (dx * dy + dy * dz + dz * dx);
}

// Transform an AABB with a given transformation matrix
// This function transforms all 8 corners of the AABB and computes a new AABB
inline AABB transform_aabb(const AABB& bbox, const Mat4f& transform) {
    // Get the 8 corners of the original AABB
    fl corners_x[8] = {
        bbox.bbox_xmin, bbox.bbox_xmax, bbox.bbox_xmin, bbox.bbox_xmax,
        bbox.bbox_xmin, bbox.bbox_xmax, bbox.bbox_xmin, bbox.bbox_xmax
    };
    fl corners_y[8] = {
        bbox.bbox_ymin, bbox.bbox_ymin, bbox.bbox_ymax, bbox.bbox_ymax,
        bbox.bbox_ymin, bbox.bbox_ymin, bbox.bbox_ymax, bbox.bbox_ymax
    };
    fl corners_z[8] = {
        bbox.bbox_zmin, bbox.bbox_zmin, bbox.bbox_zmin, bbox.bbox_zmin,
        bbox.bbox_zmax, bbox.bbox_zmax, bbox.bbox_zmax, bbox.bbox_zmax
    };
    
    // Initialize the transformed AABB with extreme values
    AABB transformed_bbox;
    transformed_bbox.bbox_xmin = INFINITY;
    transformed_bbox.bbox_ymin = INFINITY;
    transformed_bbox.bbox_zmin = INFINITY;
    transformed_bbox.bbox_xmax = -INFINITY;
    transformed_bbox.bbox_ymax = -INFINITY;
    transformed_bbox.bbox_zmax = -INFINITY;
    
    // Transform each corner and expand the new AABB
    for (int i = 0; i < 8; i++) {
        fl x = corners_x[i];
        fl y = corners_y[i];
        fl z = corners_z[i];
        
 // Apply transformation: P' = M * P
fl transformed_x = transform.m11 * x + transform.m12 * y + transform.m13 * z + transform.m14;
fl transformed_y = transform.m21 * x + transform.m22 * y + transform.m23 * z + transform.m24;
fl transformed_z = transform.m31 * x + transform.m32 * y + transform.m33 * z + transform.m34;
fl transformed_w = transform.m41 * x + transform.m42 * y + transform.m43 * z + transform.m44;
        
        // Homogeneous divide (if w != 1)
        if (transformed_w != 0.0f && transformed_w != 1.0f) {
            transformed_x /= transformed_w;
            transformed_y /= transformed_w;
            transformed_z /= transformed_w;
        }
        
        // Expand the transformed AABB to include this transformed corner
        transformed_bbox.bbox_xmin = std::min(transformed_bbox.bbox_xmin, transformed_x);
        transformed_bbox.bbox_ymin = std::min(transformed_bbox.bbox_ymin, transformed_y);
        transformed_bbox.bbox_zmin = std::min(transformed_bbox.bbox_zmin, transformed_z);
        transformed_bbox.bbox_xmax = std::max(transformed_bbox.bbox_xmax, transformed_x);
        transformed_bbox.bbox_ymax = std::max(transformed_bbox.bbox_ymax, transformed_y);
        transformed_bbox.bbox_zmax = std::max(transformed_bbox.bbox_zmax, transformed_z);
    }
    
    return transformed_bbox;
}

// --- BVH CONSTRUCTION (Unchanged) ---

// Helper function to get primitive centroid
inline void get_primitive_centroid(const BVH& bvh,const Scene& scene, int prim_idx, fl& out_x, fl& out_y, fl& out_z) {
    int prim_ind = bvh.prim_indices[prim_idx];
    
    fl v0_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v0_ind[prim_ind]];
    fl v0_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v0_ind[prim_ind]];
    fl v0_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v0_ind[prim_ind]];
    
    fl v1_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v1_ind[prim_ind]];
    fl v1_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v1_ind[prim_ind]];
    fl v1_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v1_ind[prim_ind]];
    
    fl v2_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v2_ind[prim_ind]];
    fl v2_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v2_ind[prim_ind]];
    fl v2_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v2_ind[prim_ind]];
    
    out_x = (v0_x + v1_x + v2_x) / 3.0f;
    out_y = (v0_y + v1_y + v2_y) / 3.0f;
    out_z = (v0_z + v1_z + v2_z) / 3.0f;

}

inline void update_bbox(BVH&bvh,const Scene& scene, int node_index){
    auto& node = bvh.bvh_node_list[node_index];
    node.bbox.bbox_xmax = -INFINITY;
    node.bbox.bbox_ymax = -INFINITY;
    node.bbox.bbox_zmax = -INFINITY;

    node.bbox.bbox_xmin = INFINITY;
    node.bbox.bbox_ymin = INFINITY;
    node.bbox.bbox_zmin = INFINITY;

    for(int i = node.first_prim;i<node.prim_count+node.first_prim;i++){
        int prim_ind = bvh.prim_indices[i];
        const auto& v0_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v0_ind[prim_ind]];
        const auto& v0_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v0_ind[prim_ind]];
        const auto& v0_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v0_ind[prim_ind]];
        
        const auto& v1_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v1_ind[prim_ind]];
        const auto& v1_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v1_ind[prim_ind]];
        const auto& v1_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v1_ind[prim_ind]];
        
        const auto& v2_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v2_ind[prim_ind]];
        const auto& v2_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v2_ind[prim_ind]];
        const auto& v2_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v2_ind[prim_ind]];
        
        node.bbox.bbox_xmin = std::min({node.bbox.bbox_xmin, v0_x, v1_x, v2_x});
        node.bbox.bbox_ymin = std::min({node.bbox.bbox_ymin, v0_y, v1_y, v2_y});
        node.bbox.bbox_zmin = std::min({node.bbox.bbox_zmin, v0_z, v1_z, v2_z});
        
        node.bbox.bbox_xmax = std::max({node.bbox.bbox_xmax, v0_x, v1_x, v2_x});
        node.bbox.bbox_ymax = std::max({node.bbox.bbox_ymax, v0_y, v1_y, v2_y});
        node.bbox.bbox_zmax = std::max({node.bbox.bbox_zmax, v0_z, v1_z, v2_z});
    }
}

inline void subdivide_bvh(BVH&bvh,const Scene& scene, int node_index){
    auto& node = bvh.bvh_node_list[node_index];
    if(node.prim_count<=2) return;

    const int NUM_SAMPLES = 5;
    fl best_cost= INFINITY;
    int best_axis = -1;
    fl best_split_pos = 0.0f;
    fl parent_surface_area = calculate_surface_area(node.bbox);


    for (int axis = 0; axis < 3; axis++){
        fl min_coord = (axis == 0) ? node.bbox.bbox_xmin : (axis == 1) ? node.bbox.bbox_ymin : node.bbox.bbox_zmin;
        fl max_coord = (axis == 0) ? node.bbox.bbox_xmax : (axis == 1) ? node.bbox.bbox_ymax : node.bbox.bbox_zmax;
        
         fl extent = max_coord - min_coord;
        if (extent < 1e-6f) continue; // Skip degenerate axis
        for (int sample = 1; sample < NUM_SAMPLES; sample++){
            fl t = (fl)sample / (fl)NUM_SAMPLES;
            fl split_pos = min_coord + t * extent;
            
            // Count primitives and compute bounding boxes for left and right partitions
            int left_count = 0;
            int right_count = 0;
            AABB left_bbox, right_bbox;
            left_bbox.bbox_xmax = left_bbox.bbox_ymax = left_bbox.bbox_zmax = -INFINITY;
            left_bbox.bbox_xmin = left_bbox.bbox_ymin = left_bbox.bbox_zmin = INFINITY;
            right_bbox.bbox_xmax = right_bbox.bbox_ymax = right_bbox.bbox_zmax = -INFINITY;
            right_bbox.bbox_xmin = right_bbox.bbox_ymin = right_bbox.bbox_zmin = INFINITY;
            
            for (int i = node.first_prim; i < node.first_prim + node.prim_count; i++){
                fl cx,cy,cz;
                get_primitive_centroid(bvh,scene,i,cx,cy,cz);
                
                fl centroid_coord = (axis == 0) ? cx : (axis == 1) ? cy : cz;
                
                int prim_ind = bvh.prim_indices[i];
                
                AABB prim_bbox;
                
                fl v0_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v0_ind[prim_ind]];
                fl v0_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v0_ind[prim_ind]];
                fl v0_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v0_ind[prim_ind]];
                fl v1_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v1_ind[prim_ind]];
                fl v1_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v1_ind[prim_ind]];
                fl v1_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v1_ind[prim_ind]];
                fl v2_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v2_ind[prim_ind]];
                fl v2_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v2_ind[prim_ind]];
                fl v2_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v2_ind[prim_ind]];
                
                prim_bbox.bbox_xmin = std::min({v0_x, v1_x, v2_x});
                prim_bbox.bbox_ymin = std::min({v0_y, v1_y, v2_y});
                prim_bbox.bbox_zmin = std::min({v0_z, v1_z, v2_z});
                prim_bbox.bbox_xmax = std::max({v0_x, v1_x, v2_x});
                prim_bbox.bbox_ymax = std::max({v0_y, v1_y, v2_y});
                prim_bbox.bbox_zmax = std::max({v0_z, v1_z, v2_z});

                if (centroid_coord < split_pos) {
                    left_count++;
                    left_bbox.bbox_xmin = std::min(left_bbox.bbox_xmin, prim_bbox.bbox_xmin);
                    left_bbox.bbox_ymin = std::min(left_bbox.bbox_ymin, prim_bbox.bbox_ymin);
                    left_bbox.bbox_zmin = std::min(left_bbox.bbox_zmin, prim_bbox.bbox_zmin);
                    left_bbox.bbox_xmax = std::max(left_bbox.bbox_xmax, prim_bbox.bbox_xmax);
                    left_bbox.bbox_ymax = std::max(left_bbox.bbox_ymax, prim_bbox.bbox_ymax);
                    left_bbox.bbox_zmax = std::max(left_bbox.bbox_zmax, prim_bbox.bbox_zmax);
                } else {
                    right_count++;
                    right_bbox.bbox_xmin = std::min(right_bbox.bbox_xmin, prim_bbox.bbox_xmin);
                    right_bbox.bbox_ymin = std::min(right_bbox.bbox_ymin, prim_bbox.bbox_ymin);
                    right_bbox.bbox_zmin = std::min(right_bbox.bbox_zmin, prim_bbox.bbox_zmin);
                    right_bbox.bbox_xmax = std::max(right_bbox.bbox_xmax, prim_bbox.bbox_xmax);
                    right_bbox.bbox_ymax = std::max(right_bbox.bbox_ymax, prim_bbox.bbox_ymax);
                    right_bbox.bbox_zmax = std::max(right_bbox.bbox_zmax, prim_bbox.bbox_zmax);
                }
            }

            if (left_count == 0 || right_count == 0) continue;
            
            fl left_surface_area = calculate_surface_area(left_bbox);
            fl right_surface_area = calculate_surface_area(right_bbox);
            fl cost = left_surface_area * left_count + right_surface_area * right_count;
            
            if (cost < best_cost) {
                best_cost = cost;
                best_axis = axis;
                best_split_pos = split_pos;
            }
        }
    }

    if (best_axis == -1) return;
    
    // Partition primitives using the best split found
    int left = node.first_prim;
    int right = node.first_prim + node.prim_count - 1;
    
    while (left <= right) {
        fl cx, cy, cz;
        get_primitive_centroid(bvh,scene, left, cx, cy, cz);
        fl centroid_coord = (best_axis == 0) ? cx : (best_axis == 1) ? cy : cz;
        
        if (centroid_coord < best_split_pos) {
            left++;
        } else {
            std::swap(bvh.prim_indices[left], bvh.prim_indices[right]);
            right--;
        }
    }
    
    int left_count = left - node.first_prim;
    int right_count = node.prim_count - left_count;
    
    // Avoid empty partitions (shouldn't happen but safety check)
    if (left_count == 0 || right_count == 0) return;
    
    // Create child nodes
    int left_index = bvh.next_node++;
    int right_index = bvh.next_node++;
    
    node.left_child_index = left_index;
    
    // Left child
    bvh.bvh_node_list[left_index].first_prim = node.first_prim;
    bvh.bvh_node_list[left_index].prim_count = left_count;
    bvh.bvh_node_list[left_index].left_child_index = 0;
    
    // Right child
    bvh.bvh_node_list[right_index].first_prim = node.first_prim + left_count;
    bvh.bvh_node_list[right_index].prim_count = right_count;
    bvh.bvh_node_list[right_index].left_child_index = 0;
    
    // Update bounding boxes
    update_bbox(bvh,scene, left_index);
    update_bbox(bvh,scene, right_index);
    
    // Recursively subdivide
    subdivide_bvh(bvh, scene, left_index);
    subdivide_bvh(bvh, scene, right_index);
}

//constructed per mesh, at the parser mesh section. 
//NO SEPERATE CONSTRUCTION FOR SPHERES
inline void construct_BVH(const MeshInfo& mesh,BVH& bvh,const Scene& scene ){
    //populate the prim indices while traversing the mesh::DANGER
    bvh.bvh_node_list = (BVHNode*)malloc((bvh.prim_indices.size()*2-1)*sizeof(BVHNode));

    BVHNode&root = bvh.bvh_node_list[0];
    root.left_child_index = 0;
    root.first_prim=0;root.prim_count = mesh.triangle_count;
    bvh.next_node = 1; // Explicitly set next_node
    update_bbox(bvh,scene,0);
    subdivide_bvh(bvh,scene,0);
}

// --- BLAS CONSTRUCTION (Unchanged) ---

inline void construct_BLAS(const MeshInfo& mesh,BLAS& blas,const Scene& scene){
    const auto& bvh = bvh_list[mesh.bvh_index];
    const auto& bvh_bbox = bvh.bvh_node_list[0].bbox;
    blas.inv_transform = mesh.inverse_transformation_matrix;
    blas.transformation_matrix = mesh.transformation_matrix;
    blas.bbox = transform_aabb(bvh_bbox,mesh.transformation_matrix);
    blas.bvh_index = mesh.bvh_index;
    blas.material_id = mesh.material_id;
    }

// --- NEW TLAS CONSTRUCTION ---

/**
 * @brief Gets the centroid of a BLAS's AABB.
 */
inline void get_blas_centroid(int blas_index, fl& out_x, fl& out_y, fl& out_z) {
    const AABB& bbox = blas_list[blas_index].bbox;
    out_x = (bbox.bbox_xmin + bbox.bbox_xmax) * 0.5f;
    out_y = (bbox.bbox_ymin + bbox.bbox_ymax) * 0.5f;
    out_z = (bbox.bbox_zmin + bbox.bbox_zmax) * 0.5f;
}

/**
 * @brief Updates a TLASNode's AABB to encompass all BLASes in its range.
 * @param node_index The index of the TLAS node to update.
 * @param start The starting index in the global `tlas_blas_indices` vector.
 * @param count The number of BLASes in this node's range.
 */
inline void update_tlas_bbox(int node_index, int start, int count) {
    auto& node = tlas.tlas_root[node_index];
    node.bbox.bbox_xmax = -INFINITY;
    node.bbox.bbox_ymax = -INFINITY;
    node.bbox.bbox_zmax = -INFINITY;
    node.bbox.bbox_xmin = INFINITY;
    node.bbox.bbox_ymin = INFINITY;
    node.bbox.bbox_zmin = INFINITY;

    for (int i = start; i < start + count; i++) {
        int blas_index = tlas_blas_indices[i];
        const AABB& blas_bbox = blas_list[blas_index].bbox;

        node.bbox.bbox_xmin = std::min(node.bbox.bbox_xmin, blas_bbox.bbox_xmin);
        node.bbox.bbox_ymin = std::min(node.bbox.bbox_ymin, blas_bbox.bbox_ymin);
        node.bbox.bbox_zmin = std::min(node.bbox.bbox_zmin, blas_bbox.bbox_zmin);
        node.bbox.bbox_xmax = std::max(node.bbox.bbox_xmax, blas_bbox.bbox_xmax);
        node.bbox.bbox_ymax = std::max(node.bbox.bbox_ymax, blas_bbox.bbox_ymax);
        node.bbox.bbox_zmax = std::max(node.bbox.bbox_zmax, blas_bbox.bbox_zmax);
    }
}

/**
 * @brief Recursively subdivides a TLAS node using SAH, following the BVH procedure.
 * @param node_index The index of the TLAS node to subdivide.
 * @param start The starting index in the global `tlas_blas_indices` vector.
 * @param count The number of BLASes in this node's range.
 */
inline void subdivide_TLAS(int node_index, int start, int count) {
    auto& node = tlas.tlas_root[node_index];

    // Base case: 1 or 2 BLASes. Make a leaf node.
    // This is the primary difference from subdivide_bvh.
    if (count <= 2) {
        node.left_child_index = 0; // Mark as leaf
        node.blas_id_0 = tlas_blas_indices[start];
        node.blas_id_1 = (count == 2) ? tlas_blas_indices[start + 1] : -1;
        return;
    }

    // Recursive case: Find the best split using SAH (same as subdivide_bvh)
    const int NUM_SAMPLES = 5;
    fl best_cost = INFINITY;
    int best_axis = -1;
    fl best_split_pos = 0.0f;
    fl parent_surface_area = calculate_surface_area(node.bbox);

    for (int axis = 0; axis < 3; axis++) {
        fl min_coord = (axis == 0) ? node.bbox.bbox_xmin : (axis == 1) ? node.bbox.bbox_ymin : node.bbox.bbox_zmin;
        fl max_coord = (axis == 0) ? node.bbox.bbox_xmax : (axis == 1) ? node.bbox.bbox_ymax : node.bbox.bbox_zmax;

        fl extent = max_coord - min_coord;
        if (extent < 1e-6f) continue;

        for (int sample = 1; sample < NUM_SAMPLES; sample++) {
            fl t = (fl)sample / (fl)NUM_SAMPLES;
            fl split_pos = min_coord + t * extent;

            int left_count = 0;
            int right_count = 0;
            AABB left_bbox, right_bbox;
            left_bbox.bbox_xmax = left_bbox.bbox_ymax = left_bbox.bbox_zmax = -INFINITY;
            left_bbox.bbox_xmin = left_bbox.bbox_ymin = left_bbox.bbox_zmin = INFINITY;
            right_bbox.bbox_xmax = right_bbox.bbox_ymax = right_bbox.bbox_zmax = -INFINITY;
            right_bbox.bbox_xmin = right_bbox.bbox_ymin = right_bbox.bbox_zmin = INFINITY;

            // Partition test (same logic as BVH, but on BLASes)
            for (int i = start; i < start + count; i++) {
                int blas_index = tlas_blas_indices[i];
                fl cx, cy, cz;
                get_blas_centroid(blas_index, cx, cy, cz);
                fl centroid_coord = (axis == 0) ? cx : (axis == 1) ? cy : cz;
                
                const AABB& blas_bbox = blas_list[blas_index].bbox;

                if (centroid_coord < split_pos) {
                    left_count++;
                    left_bbox.bbox_xmin = std::min(left_bbox.bbox_xmin, blas_bbox.bbox_xmin);
                    left_bbox.bbox_ymin = std::min(left_bbox.bbox_ymin, blas_bbox.bbox_ymin);
                    left_bbox.bbox_zmin = std::min(left_bbox.bbox_zmin, blas_bbox.bbox_zmin);
                    left_bbox.bbox_xmax = std::max(left_bbox.bbox_xmax, blas_bbox.bbox_xmax);
                    left_bbox.bbox_ymax = std::max(left_bbox.bbox_ymax, blas_bbox.bbox_ymax);
                    left_bbox.bbox_zmax = std::max(left_bbox.bbox_zmax, blas_bbox.bbox_zmax);
                } else {
                    right_count++;
                    right_bbox.bbox_xmin = std::min(right_bbox.bbox_xmin, blas_bbox.bbox_xmin);
                    right_bbox.bbox_ymin = std::min(right_bbox.bbox_ymin, blas_bbox.bbox_ymin);
                    right_bbox.bbox_zmin = std::min(right_bbox.bbox_zmin, blas_bbox.bbox_zmin);
                    right_bbox.bbox_xmax = std::max(right_bbox.bbox_xmax, blas_bbox.bbox_xmax);
                    right_bbox.bbox_ymax = std::max(right_bbox.bbox_ymax, blas_bbox.bbox_ymax);
                    right_bbox.bbox_zmax = std::max(right_bbox.bbox_zmax, blas_bbox.bbox_zmax);
                }
            }

            if (left_count == 0 || right_count == 0) continue;
            fl left_surface_area = calculate_surface_area(left_bbox);
            fl right_surface_area = calculate_surface_area(right_bbox);
            fl cost = left_surface_area * left_count + right_surface_area * right_count;

            if (cost < best_cost) {
                best_cost = cost;
                best_axis = axis;
                best_split_pos = split_pos;
            }
        }
    }

    // If SAH fails to find a split, we MUST still subdivide
    // (unlike BVH) because our leaves can't hold > 2 items.
    // We'll fall back to a simple median split.
    if (best_axis == -1) {
        int left = start;
        int right = start + count - 1;
        int mid_point = start + (count / 2);

        // Simple median partition
        while (left < mid_point) {
            left++;
        }
        while (right >= mid_point) {
            std::swap(tlas_blas_indices[left], tlas_blas_indices[right]);
            left++;
            right--;
        }
    } else {
        // Partition BLAS indices using the best split found (same as BVH)
        int left = start;
        int right = start + count - 1;
        while (left <= right) {
            fl cx, cy, cz;
            get_blas_centroid(tlas_blas_indices[left], cx, cy, cz);
            fl centroid_coord = (best_axis == 0) ? cx : (best_axis == 1) ? cy : cz;

            if (centroid_coord < best_split_pos) {
                left++;
            } else {
                std::swap(tlas_blas_indices[left], tlas_blas_indices[right]);
                right--;
            }
        }
        
        int left_count = left - start;
        int right_count = count - left_count;
        
        // If partition failed, fall back to median
        if (left_count == 0 || right_count == 0) {
            // This re-does the partition from the `best_axis == -1` block
            left = start;
            right = start + count - 1;
            int mid_point = start + (count / 2);
            while (left < mid_point) left++;
            while (right >= mid_point) {
                std::swap(tlas_blas_indices[left], tlas_blas_indices[right]);
                left++;
                right--;
            }
        }
    }
    
    // This logic is slightly different from the BVH one to find the partition
    // We just need the final split point from the `left` index.
    int left = start;
    if (best_axis != -1) {
        int right = start + count - 1;
        while (left <= right) {
             fl cx, cy, cz;
            get_blas_centroid(tlas_blas_indices[left], cx, cy, cz);
            fl centroid_coord = (best_axis == 0) ? cx : (best_axis == 1) ? cy : cz;
            if (centroid_coord < best_split_pos) left++;
            else { std::swap(tlas_blas_indices[left], tlas_blas_indices[right]); right--; }
        }
    }
    int left_count = left - start;
    int right_count = count - left_count;
    
    // If SAH failed or partition failed, force a median split
    if (best_axis == -1 || left_count == 0 || right_count == 0) {
        left = start + (count / 2);
        left_count = left - start;
        right_count = count - left_count;
        
        // Need to re-partition the array
        std::nth_element(&tlas_blas_indices[start], 
                         &tlas_blas_indices[left], 
                         &tlas_blas_indices[start + count],
                         [](int a, int b) {
                            fl ax, ay, az, bx, by, bz;
                            get_blas_centroid(a, ax, ay, az);
                            get_blas_centroid(b, bx, by, bz);
                            // Just pick an axis, e.g., x
                            return ax < bx;
                         });
    }

    // Create child nodes
    int left_index = tlas.next_node++;
    int right_index = tlas.next_node++;
    
    node.left_child_index = left_index;

    // Left child
    update_tlas_bbox(left_index, start, left_count);
    
    // Right child
    update_tlas_bbox(right_index, left, right_count);
    
    // Recursively subdivide
    subdivide_TLAS(left_index, start, left_count);
    subdivide_TLAS(right_index, left, right_count);
}

/**
 * @brief Constructs the Top-Level Acceleration Structure (TLAS) from the
 * global `blas_list`.
 */
inline void construct_TLAS(){
    tlas.next_node = 0; // Start with root node 0
    int num_blas = blas_list.size();
    if (num_blas == 0) return;

    // Allocate memory for all potential nodes
    tlas.tlas_root = (TLASNode*)malloc((num_blas * 2 - 1) * sizeof(TLASNode));

    // Initialize the global BLAS indices vector
    tlas_blas_indices.resize(num_blas);
    for (int i = 0; i < num_blas; i++) {
        tlas_blas_indices[i] = i;
    }

    // Initialize the root node (node 0)
    int root_index = tlas.next_node++;
    update_tlas_bbox(root_index, 0, num_blas);

    // Start subdivision
    subdivide_TLAS(root_index, 0, num_blas);
}