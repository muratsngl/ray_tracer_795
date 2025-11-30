#pragma once

#include<types.hpp>
#include<utils.hpp>
#include <algorithm>
#include <cmath>
#include <vector>

inline BVHNode* bvhNode;

inline int next_node_index = 1;

// Helper function to calculate surface area of an AABB
inline fl calculate_surface_area_1(const AABB& bbox) {
    fl dx = bbox.bbox_xmax - bbox.bbox_xmin;
    fl dy = bbox.bbox_ymax - bbox.bbox_ymin;
    fl dz = bbox.bbox_zmax - bbox.bbox_zmin;
    return 2.0f * (dx * dy + dy * dz + dz * dx);
}

// Helper function to get primitive centroid
inline void get_primitive_centroid(const Scene& scene, int prim_idx, fl& out_x, fl& out_y, fl& out_z) {
    int prim_ind = scene.prim_data__.prim_index[prim_idx];
    PrimType pt = scene.prim_data__.primitive_type[prim_idx];
    
    if (pt == TRIANGLE) {
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
    } else { // SPHERE
        out_x = scene.vertex_data__.v_pos_x[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
        out_y = scene.vertex_data__.v_pos_y[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
        out_z = scene.vertex_data__.v_pos_z[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
    }
}

// Helper function to compute AABB for a range of primitives
inline AABB compute_bbox_for_range(const Scene& scene, int start, int count) {
    AABB bbox;
    bbox.bbox_xmax = bbox.bbox_ymax = bbox.bbox_zmax = -INFINITY;
    bbox.bbox_xmin = bbox.bbox_ymin = bbox.bbox_zmin = INFINITY;
    
    for(int i = start; i < start + count; i++) {
        int prim_ind = scene.prim_data__.prim_index[i];
        PrimType prim_type = scene.prim_data__.primitive_type[i];
        
        if(prim_type == TRIANGLE) {
            fl v0_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v0_ind[prim_ind]];
            fl v0_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v0_ind[prim_ind]];
            fl v0_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v0_ind[prim_ind]];
            
            fl v1_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v1_ind[prim_ind]];
            fl v1_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v1_ind[prim_ind]];
            fl v1_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v1_ind[prim_ind]];
            
            fl v2_x = scene.vertex_data__.v_pos_x[scene.triangle_data__.v2_ind[prim_ind]];
            fl v2_y = scene.vertex_data__.v_pos_y[scene.triangle_data__.v2_ind[prim_ind]];
            fl v2_z = scene.vertex_data__.v_pos_z[scene.triangle_data__.v2_ind[prim_ind]];
            
            bbox.bbox_xmin = std::min({bbox.bbox_xmin, v0_x, v1_x, v2_x});
            bbox.bbox_ymin = std::min({bbox.bbox_ymin, v0_y, v1_y, v2_y});
            bbox.bbox_zmin = std::min({bbox.bbox_zmin, v0_z, v1_z, v2_z});
            
            bbox.bbox_xmax = std::max({bbox.bbox_xmax, v0_x, v1_x, v2_x});
            bbox.bbox_ymax = std::max({bbox.bbox_ymax, v0_y, v1_y, v2_y});
            bbox.bbox_zmax = std::max({bbox.bbox_zmax, v0_z, v1_z, v2_z});
        } else { // SPHERE
            fl center_x = scene.vertex_data__.v_pos_x[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
            fl center_y = scene.vertex_data__.v_pos_y[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
            fl center_z = scene.vertex_data__.v_pos_z[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
            
            fl radius = std::sqrt(scene.sphere_data__.sphere_radius_sq[prim_ind]);
            
            bbox.bbox_xmin = std::min(bbox.bbox_xmin, center_x - radius);
            bbox.bbox_ymin = std::min(bbox.bbox_ymin, center_y - radius);
            bbox.bbox_zmin = std::min(bbox.bbox_zmin, center_z - radius);
            
            bbox.bbox_xmax = std::max(bbox.bbox_xmax, center_x + radius);
            bbox.bbox_ymax = std::max(bbox.bbox_ymax, center_y + radius);
            bbox.bbox_zmax = std::max(bbox.bbox_zmax, center_z + radius);
        }
    }
    
    return bbox;
}

//IN ORDER TO UNDERSTAND THE LOGIC WITHOUT GETTING LOST IN THE SIMD REALM I AM WRITING THE CODE IN SINGLE OPERATION FASHION
void update_bbox(const Scene& scene,int node_index){
    auto& node = bvhNode[node_index];
    node.bbox.bbox_xmax = -INFINITY;
    node.bbox.bbox_ymax = -INFINITY;
    node.bbox.bbox_zmax = -INFINITY;

    node.bbox.bbox_xmin = INFINITY;
    node.bbox.bbox_ymin = INFINITY;
    node.bbox.bbox_zmin = INFINITY;

    
    for(int i = node.first_prim;i<node.prim_count+node.first_prim;i++){
        const auto& prim_ind = scene.prim_data__.prim_index[i];
        const auto& prim_type = scene.prim_data__.primitive_type[i]; 

        switch(prim_type)
        {
            case TRIANGLE:
            {
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
                break;
            }
            case SPHERE:
            {
                const auto& center_x = scene.vertex_data__.v_pos_x[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
                const auto& center_y = scene.vertex_data__.v_pos_y[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
                const auto& center_z = scene.vertex_data__.v_pos_z[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
                
                const auto radius = std::sqrt(scene.sphere_data__.sphere_radius_sq[prim_ind]);
                
                node.bbox.bbox_xmin = std::min(node.bbox.bbox_xmin, center_x - radius);
                node.bbox.bbox_ymin = std::min(node.bbox.bbox_ymin, center_y - radius);
                node.bbox.bbox_zmin = std::min(node.bbox.bbox_zmin, center_z - radius);
                
                node.bbox.bbox_xmax = std::max(node.bbox.bbox_xmax, center_x + radius);
                node.bbox.bbox_ymax = std::max(node.bbox.bbox_ymax, center_y + radius);
                node.bbox.bbox_zmax = std::max(node.bbox.bbox_zmax, center_z + radius);
                
                break;
            }
        }
    } 

};


void subdivide(Scene& scene, int node_index){
    auto& node = bvhNode[node_index];
    if (node.prim_count <= 2) return; // leaf node with <= 2 primitives

    // Surface Area Heuristic (SAH) based split selection
    const int NUM_SAMPLES = 5;
    fl best_cost = INFINITY;
    int best_axis = -1;
    fl best_split_pos = 0.0f;
    
    fl parent_surface_area = calculate_surface_area_1(node.bbox);
    
    // Try each axis
   
    for (int axis = 0; axis < 3; axis++) {
        fl min_coord = (axis == 0) ? node.bbox.bbox_xmin : (axis == 1) ? node.bbox.bbox_ymin : node.bbox.bbox_zmin;
        fl max_coord = (axis == 0) ? node.bbox.bbox_xmax : (axis == 1) ? node.bbox.bbox_ymax : node.bbox.bbox_zmax;
        
        fl extent = max_coord - min_coord;
        if (extent < 1e-6f) continue; // Skip degenerate axis
        
        // Test NUM_SAMPLES split positions along this axis
        for (int sample = 1; sample < NUM_SAMPLES; sample++) {
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
            
            for (int i = node.first_prim; i < node.first_prim + node.prim_count; i++) {
                fl cx, cy, cz;
                get_primitive_centroid(scene, i, cx, cy, cz);
                fl centroid_coord = (axis == 0) ? cx : (axis == 1) ? cy : cz;
                
                int prim_ind = scene.prim_data__.prim_index[i];
                PrimType prim_type = scene.prim_data__.primitive_type[i];
                
                // Compute primitive bounds
                AABB prim_bbox;
                if (prim_type == TRIANGLE) {
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
                } else { // SPHERE
                    fl center_x = scene.vertex_data__.v_pos_x[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
                    fl center_y = scene.vertex_data__.v_pos_y[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
                    fl center_z = scene.vertex_data__.v_pos_z[scene.sphere_data__.sphere_center_vertex_id[prim_ind]];
                    fl radius = std::sqrt(scene.sphere_data__.sphere_radius_sq[prim_ind]);
                    
                    prim_bbox.bbox_xmin = center_x - radius;
                    prim_bbox.bbox_ymin = center_y - radius;
                    prim_bbox.bbox_zmin = center_z - radius;
                    prim_bbox.bbox_xmax = center_x + radius;
                    prim_bbox.bbox_ymax = center_y + radius;
                    prim_bbox.bbox_zmax = center_z + radius;
                }
                
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
            
            // Skip if partition is empty
            if (left_count == 0 || right_count == 0) continue;
            
            // Calculate SAH cost: C = C_trav + (SA_left / SA_parent) * N_left * C_isect + (SA_right / SA_parent) * N_right * C_isect
            // Simplified: cost = left_area * left_count + right_area * right_count
            fl left_surface_area = calculate_surface_area_1(left_bbox);
            fl right_surface_area = calculate_surface_area_1(right_bbox);
            fl cost = left_surface_area * left_count + right_surface_area * right_count;
            
            if (cost < best_cost) {
                best_cost = cost;
                best_axis = axis;
                best_split_pos = split_pos;
            }
        }
    }
    
    // If no good split found, make this a leaf
    if (best_axis == -1) return;
    
    // Partition primitives using the best split found
    int left = node.first_prim;
    int right = node.first_prim + node.prim_count - 1;
    
    while (left <= right) {
        fl cx, cy, cz;
        get_primitive_centroid(scene, left, cx, cy, cz);
        fl centroid_coord = (best_axis == 0) ? cx : (best_axis == 1) ? cy : cz;
        
        if (centroid_coord < best_split_pos) {
            left++;
        } else {
            std::swap(scene.prim_data__.prim_index[left], scene.prim_data__.prim_index[right]);
            std::swap(scene.prim_data__.primitive_type[left], scene.prim_data__.primitive_type[right]);
            right--;
        }
    }
    
    int left_count = left - node.first_prim;
    int right_count = node.prim_count - left_count;
    
    // Avoid empty partitions (shouldn't happen but safety check)
    if (left_count == 0 || right_count == 0) return;
    
    // Create child nodes
    int left_index = next_node_index++;
    int right_index = next_node_index++;
    
    node.left_child_index = left_index;
    
    // Left child
    bvhNode[left_index].first_prim = node.first_prim;
    bvhNode[left_index].prim_count = left_count;
    bvhNode[left_index].left_child_index = 0;
    
    // Right child
    bvhNode[right_index].first_prim = node.first_prim + left_count;
    bvhNode[right_index].prim_count = right_count;
    bvhNode[right_index].left_child_index = 0;
    
    // Update bounding boxes
    update_bbox(scene, left_index);
    update_bbox(scene, right_index);
    
    // Recursively subdivide
    subdivide(scene, left_index);
    subdivide(scene, right_index);
}

void construct_bvh_node(Scene& scene)
{
    bvhNode = (BVHNode*)malloc(scene.primitive_count*sizeof(BVHNode)*2-1);
    next_node_index = 1;

    BVHNode& root = bvhNode[0];
    root.left_child_index = 0;
    root.first_prim=0,root.prim_count = scene.primitive_count;
    update_bbox(scene,0);
    subdivide(scene,0);
}



