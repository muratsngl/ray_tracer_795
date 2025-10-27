#pragma once
#ifndef RTFUNC
#define RTFUNC


#include <xsimd/xsimd.hpp>
#include "utils.hpp" // Assuming this contains operator overloads for Vec3f
#include "types.hpp"

namespace xs = xsimd;

// Acronym type aliases for xsimd batches
using f_batch = xs::batch<fl>;
using i_batch = xs::batch<int>;
using b_batch = xs::batch_bool<fl>;

void inline intersect_planes(const f_batch& dir_x,
    const f_batch& dir_y,
    const f_batch& dir_z,
    const f_batch& orig_x,
    const f_batch& orig_y,
    const f_batch& orig_z,
    f_batch& t_min,
    i_batch& mat_id,
    // Add normal output batches
    f_batch& norm_x_out,
    f_batch& norm_y_out,
    f_batch& norm_z_out,
    const PlaneData& planes, const VertexData& vertices, const Scene& scene, int plane_ID) {

    // 1. Get plane data (normal N, point P0) and broadcast to 8-wide batches
    auto norm_x = xs::broadcast(planes.plane_norm_x[plane_ID]);
    auto norm_y = xs::broadcast(planes.plane_norm_y[plane_ID]);
    auto norm_z = xs::broadcast(planes.plane_norm_z[plane_ID]);

    int p_v_id = planes.plane_point_vertex_id[plane_ID];
    auto p0_x = xs::broadcast(vertices.v_pos_x[p_v_id]);
    auto p0_y = xs::broadcast(vertices.v_pos_y[p_v_id]);
    auto p0_z = xs::broadcast(vertices.v_pos_z[p_v_id]);

    // 2. Calculate denominator: denom = dot(RayDirection, PlaneNormal)
    auto denom = (dir_x * norm_x) + (dir_y * norm_y) + (dir_z * norm_z);

    // 3. Epsilon check for parallel rays
    const fl epsilon = scene.intersection_test_epsilon;
    b_batch active_mask = (xs::abs(denom) > epsilon);

    if (xs::none(active_mask)) {
        return; // All 8 rays are parallel to the plane
    }

    // 4. Calculate numerator: numer = dot(PlanePoint - RayOrigin, PlaneNormal)
    auto p0_minus_o_x = p0_x - orig_x;
    auto p0_minus_o_y = p0_y - orig_y;
    auto p0_minus_o_z = p0_z - orig_z;

    auto numer = (p0_minus_o_x * norm_x) + (p0_minus_o_y * norm_y) + (p0_minus_o_z * norm_z);

    // 5. Calculate t = numer / denom
    auto t_new = numer / denom;

    // 6. Create final hit mask
    b_batch final_hit_mask = active_mask & (t_new > epsilon) & (t_new < t_min);

    if (xs::none(final_hit_mask)) {
        return; // No new valid hits
    }

    // --- Update Hit Information ---

    // 7. Update normals for the rays that hit
    // The normal is constant for the entire plane
    norm_x_out = xs::select(final_hit_mask, norm_x, norm_x_out);
    norm_y_out = xs::select(final_hit_mask, norm_y, norm_y_out);
    norm_z_out = xs::select(final_hit_mask, norm_z, norm_z_out);

    // 8. Update t_min for the rays that hit
    t_min = xs::select(final_hit_mask, t_new, t_min);

    // 9. Update material ID for the rays that hit
    auto int_mask = xs::batch_bool_cast<int32_t>(final_hit_mask);
    i_batch plane_mat_id_batch = xs::broadcast(planes.plane_material_id[plane_ID]);
    mat_id = xs::select(int_mask, plane_mat_id_batch, mat_id);
}

/**
 * @brief Intersects a packet of 8 rays with a single sphere (SoA).
 */
void inline intersect_spheres(
    const f_batch& dir_x,
    const f_batch& dir_y,
    const f_batch& dir_z,
    const f_batch& orig_x,
    const f_batch& orig_y,
    const f_batch& orig_z,
    f_batch& t_min,
    i_batch& mat_id,
    // Add normal output batches
    f_batch& norm_x_out,
    f_batch& norm_y_out,
    f_batch& norm_z_out,
    const SphereData& spheres, const VertexData& vertices, const Scene& scene, const int Sphere_ID) {

    int center_vertex_id = spheres.sphere_center_vertex_id[Sphere_ID];
    auto sphere_center_x = xs::broadcast(vertices.v_pos_x[center_vertex_id]);
    auto sphere_center_y = xs::broadcast(vertices.v_pos_y[center_vertex_id]);
    auto sphere_center_z = xs::broadcast(vertices.v_pos_z[center_vertex_id]);
    auto sphere_radius_sq = xs::broadcast(spheres.sphere_radius_sq[Sphere_ID]);

    auto oc_x = orig_x - sphere_center_x;
    auto oc_y = orig_y - sphere_center_y;
    auto oc_z = orig_z - sphere_center_z;

    auto b_half = (dir_x * oc_x) + (dir_y * oc_y) + (dir_z * oc_z);

    auto c = (oc_x * oc_x) + (oc_y * oc_y) + (oc_z * oc_z) - sphere_radius_sq;

    auto discriminant = (b_half * b_half) - c;

    b_batch hit_mask = (discriminant >= 0.f);

    if (xs::none(hit_mask)) {
        return;
    }

    auto sqrt_discriminant = xs::sqrt(discriminant);

    auto t0 = -b_half - sqrt_discriminant;
    auto t1 = -b_half + sqrt_discriminant;

    auto t_smaller = xs::min(t0, t1);
    auto t_larger = xs::max(t0, t1);

    const fl epsilon = scene.intersection_test_epsilon;

    auto t_new = xs::select(t_smaller > epsilon, t_smaller, t_larger);

    b_batch final_hit_mask = hit_mask & (t_new > epsilon) & (t_new < t_min);

    if (xs::none(final_hit_mask)) {
        return;
    }

    // --- Update Hit Information ---
    
    // 1. Calculate hit point: P = O + t*D
    auto hit_x = orig_x + t_new * dir_x;
    auto hit_y = orig_y + t_new * dir_y;
    auto hit_z = orig_z + t_new * dir_z;

    // 2. Calculate normal: N = (P - C) / R
    auto n_x = hit_x - sphere_center_x;
    auto n_y = hit_y - sphere_center_y;
    auto n_z = hit_z - sphere_center_z;

    // 3. Normalize
    auto radius = xs::sqrt(sphere_radius_sq);
    f_batch epsilon_batch(1e-8f);
    b_batch radius_valid = (radius > epsilon_batch);
    auto inv_radius = xs::select(radius_valid, 1.0f / radius, f_batch(0.0f));
    
    auto new_norm_x = n_x * inv_radius;
    auto new_norm_y = n_y * inv_radius;
    auto new_norm_z = n_z * inv_radius;

    // 4. Conditionally store the new normal
    norm_x_out = xs::select(final_hit_mask, new_norm_x, norm_x_out);
    norm_y_out = xs::select(final_hit_mask, new_norm_y, norm_y_out);
    norm_z_out = xs::select(final_hit_mask, new_norm_z, norm_z_out);

    // 5. Update t_min
    t_min = xs::select(final_hit_mask, t_new, t_min);

    // 6. Update material ID
    auto int_mask = xs::batch_bool_cast<int32_t>(final_hit_mask);
    i_batch sphere_id_batch = xs::broadcast(spheres.sphere_mat_id[Sphere_ID]);
    mat_id = xs::select(int_mask, sphere_id_batch, mat_id);
}

/**
 * @brief Intersects a packet of 8 rays with a single triangle (SoA).
 * Uses Möller–Trumbore intersection algorithm.
 */
// ... (in rt_functions.hpp)

/**
 * @brief Intersects a packet of 8 rays with a single triangle (SoA).
 * Uses Möller–Trumbore intersection algorithm.
 */
void inline intersect_triangles(
    const f_batch& dir_x,
    const f_batch& dir_y,
    const f_batch& dir_z,
    const f_batch& orig_x,
    const f_batch& orig_y,
    const f_batch& orig_z,
    f_batch& t_min,
    i_batch& mat_id,
    // Add normal output batches
    f_batch& norm_x_out,
    f_batch& norm_y_out,
    f_batch& norm_z_out,
    const TriangleData& triangles, const VertexData& vertices, const Scene& scene, const int Triangle_ID) {


    // Get vertex indices
    int i0 = triangles.v0_ind[Triangle_ID];
    int i1 = triangles.v1_ind[Triangle_ID];
    int i2 = triangles.v2_ind[Triangle_ID];

    // Get triangle normal for backface culling
    auto tri_norm_x = xs::broadcast(triangles.tri_norm_x[Triangle_ID]);
    auto tri_norm_y = xs::broadcast(triangles.tri_norm_y[Triangle_ID]);
    auto tri_norm_z = xs::broadcast(triangles.tri_norm_z[Triangle_ID]);
    
    // Backface culling: dot(normal, ray_direction) should be negative (facing camera)
    auto dot_normal_ray = (tri_norm_x * dir_x) + (tri_norm_y * dir_y) + (tri_norm_z * dir_z);
    b_batch backface_mask = (dot_normal_ray < 0.0f);
    
    if (xs::none(backface_mask)) {
        return; // All rays are hitting backfaces, skip this triangle
    }

    // Get vertex positions
    Vec3f a = { vertices.v_pos_x[i0], vertices.v_pos_y[i0], vertices.v_pos_z[i0] };
    Vec3f b = { vertices.v_pos_x[i1], vertices.v_pos_y[i1], vertices.v_pos_z[i1] };
    Vec3f c = { vertices.v_pos_x[i2], vertices.v_pos_y[i2], vertices.v_pos_z[i2] };

    // Calculate edges
    Vec3f e1_scalar = b - a;
    Vec3f e2_scalar = c - a;

    // Broadcast edges to 8-wide batches
    auto e1_x = xs::broadcast(e1_scalar.x);
    auto e1_y = xs::broadcast(e1_scalar.y);
    auto e1_z = xs::broadcast(e1_scalar.z);

    auto e2_x = xs::broadcast(e2_scalar.x);
    auto e2_y = xs::broadcast(e2_scalar.y);
    auto e2_z = xs::broadcast(e2_scalar.z);

    // pvec = cross(dir, e2)
    auto pvec_x = dir_y * e2_z - dir_z * e2_y;
    auto pvec_y = dir_z * e2_x - dir_x * e2_z;
    auto pvec_z = dir_x * e2_y - dir_y * e2_x;

    // det = dot(e1, pvec)
    auto det = e1_x * pvec_x + e1_y * pvec_y + e1_z * pvec_z;

    // Epsilon check - use absolute value to allow both front and back facing triangles
    const fl epsilon = 1e-12f;
    
    // --- FIX 1: Define a slightly larger epsilon for barycentric coordinates ---
    // ################### THIS IS THE CHANGE ###################
    const fl bary_epsilon = 1e-5f; // MUCH more forgiving (was 1e-5f)
    // ########################################################
    f_batch zero_batch(0.0f);
    f_batch one_batch(1.0f);
    // --- END FIX 1 ---

    // Check determinant using abs() to disable backface culling
    b_batch active_mask = (det > epsilon);

    if (xs::none(active_mask)) {
        return;
    }

    // Safe division - only divide where det is valid
    auto invDet = xs::select(active_mask, 1.0f / det, f_batch(0.0f));

    // tvec = orig - a
    auto a_x = xs::broadcast(a.x);
    auto a_y = xs::broadcast(a.y);
    auto a_z = xs::broadcast(a.z);
    auto tvec_x = orig_x - a_x;
    auto tvec_y = orig_y - a_y;
    auto tvec_z = orig_z - a_z;

    // u = dot(tvec, pvec) * invDet
    auto u = (tvec_x * pvec_x + tvec_y * pvec_y + tvec_z * pvec_z) * invDet;

    // Mask out rays with invalid barycentric coord u
    // --- FIX 2: Use bary_epsilon for the 'u' check ---
    active_mask = active_mask & (u >= zero_batch - bary_epsilon) & (u <= one_batch + bary_epsilon);
    // --- END FIX 2 ---
    if (xs::none(active_mask)) return;

    // qvec = cross(tvec, e1)
    auto qvec_x = tvec_y * e1_z - tvec_z * e1_y;
    auto qvec_y = tvec_z * e1_x - tvec_x * e1_z;
    auto qvec_z = tvec_x * e1_y - tvec_y * e1_x;

    // v = dot(dir, qvec) * invDet
    auto v = (dir_x * qvec_x + dir_y * qvec_y + dir_z * qvec_z) * invDet;

    // Mask out rays with invalid barycentric coord v
    // --- FIX 3: Use bary_epsilon for the 'v' and 'u+v' check ---
    active_mask = active_mask & (v >= zero_batch - bary_epsilon) & (u + v <= one_batch + bary_epsilon);
    // --- END FIX 3 ---
    if (xs::none(active_mask)) return;

    // t = dot(e2, qvec) * invDet
    auto t_new = (e2_x * qvec_x + e2_y * qvec_y + e2_z * qvec_z) * invDet;

    // Final mask: active, t > scene epsilon, and t < t_min
    // --- FIX 4: Check 't_new > epsilon', not 'fabs(t_new) > 1e-7' ---
    b_batch final_hit_mask = active_mask & (t_new > epsilon) & (t_new < t_min);
    // --- END FIX 4 ---

    if (xs::none(final_hit_mask)) {
        return;
    }

    // --- Update Hit Information ---
    // Use flat triangle normal (no interpolation)
    
    // 1. Get the precomputed triangle normal
    auto new_norm_x = tri_norm_x;
    auto new_norm_y = tri_norm_y;
    auto new_norm_z = tri_norm_z;
    
    // 2. Conditionally store the triangle normal
    norm_x_out = xs::select(final_hit_mask, new_norm_x, norm_x_out);
    norm_y_out = xs::select(final_hit_mask, new_norm_y, norm_y_out);
    norm_z_out = xs::select(final_hit_mask, new_norm_z, norm_z_out);

    // 3. Update t_min
    t_min = xs::select(final_hit_mask, t_new, t_min);

    // 4. Update material ID
    auto int_mask = xs::batch_bool_cast<int32_t>(final_hit_mask);
    auto mat_id_batch = xs::broadcast(triangles.triangle_material_id[Triangle_ID]);
    mat_id = xs::select(int_mask, mat_id_batch, mat_id);
}

// ... (rest of rt_functions.hpp)
/**
 * @brief Finds the closest intersection for a packet of 8 rays against all objects in the scene.
 * @param ray_pack (Input/Output) The ray packet. t_min, mat_id, and hit_norm_x/y/z will be updated.
 * @param scene The scene data.
 * @param vertices The vertex data.
 */
void return_closest_hit(RP8& ray_pack, const Scene& scene, const VertexData& vertices) {
    // 1. Load the UNNORMALIZED data
    auto dir_x = xs::load(&ray_pack.d_x[0]);
    auto dir_y = xs::load(&ray_pack.d_y[0]);
    auto dir_z = xs::load(&ray_pack.d_z[0]);

    // 2. Calculate length-squared for all 8 vectors
    auto len_sq = (dir_x * dir_x) + (dir_y * dir_y) + (dir_z * dir_z);

    // 3. Calculate length      
    auto len = xs::sqrt(len_sq);

    // 4. Create a mask to prevent division by zero and normalize
    b_batch valid_mask = (len > 1e-8f);
    f_batch zero_batch(0.0f);
    dir_x = xs::select(valid_mask, dir_x / len, zero_batch);
    dir_y = xs::select(valid_mask, dir_y / len, zero_batch);
    dir_z = xs::select(valid_mask, dir_z / len, zero_batch);

    // 5. Load ray origins
    f_batch orig_x = xs::load(&ray_pack.o_x[0]);
    f_batch orig_y = xs::load(&ray_pack.o_y[0]);
    f_batch orig_z = xs::load(&ray_pack.o_z[0]);

    // 6. Load current intersection state (t_min, material id, and normals)
    f_batch t_min = xs::load(&ray_pack.t_min[0]);
    i_batch hit_id = xs::load(&ray_pack.mat_id[0]);
    f_batch norm_x = xs::load(&ray_pack.hit_norm_x[0]);
    f_batch norm_y = xs::load(&ray_pack.hit_norm_y[0]);
    f_batch norm_z = xs::load(&ray_pack.hit_norm_z[0]);


    // --- Iterate over all scene objects (SoA) ---

    // 7. Intersect Spheres
    int num_spheres = scene.sphere_data__.sphere_id.size();
    for (int i = 0; i < num_spheres; ++i) {
        intersect_spheres(dir_x, dir_y, dir_z, orig_x, orig_y, orig_z,
            t_min, hit_id, norm_x, norm_y, norm_z, // Pass normals
            scene.sphere_data__, vertices, scene, i);
    }

    // 8. Intersect Triangles
    int num_triangles = scene.triangle_data__.triangle_id.size();
    for (int i = 0; i < num_triangles; ++i) {
        intersect_triangles(dir_x, dir_y, dir_z, orig_x, orig_y, orig_z,
            t_min, hit_id, norm_x, norm_y, norm_z, // Pass normals
            scene.triangle_data__, vertices, scene, i);
    }

    // 9. Intersect Planes (Now using the SoA data)
    int num_planes = scene.plane_data__.plane_id.size();
    for (int i = 0; i < num_planes; ++i) {
        intersect_planes(dir_x, dir_y, dir_z, orig_x, orig_y, orig_z,
            t_min, hit_id, norm_x, norm_y, norm_z, // Pass normals
            scene.plane_data__, vertices, scene, i);
    }
    

    // 10. Store the final results back into the ray packet
    xs::store(&ray_pack.t_min[0], t_min);
    xs::store(&ray_pack.mat_id[0], hit_id);
    xs::store(&ray_pack.hit_norm_x[0], norm_x);
    xs::store(&ray_pack.hit_norm_y[0], norm_y);
    xs::store(&ray_pack.hit_norm_z[0], norm_z);

    // 11. Calculate and store hit positions: P = O + t * D
    auto hit_pos_x = orig_x + t_min * dir_x;
    auto hit_pos_y = orig_y + t_min * dir_y;
    auto hit_pos_z = orig_z + t_min * dir_z;
    
    xs::store(&ray_pack.hit_pos_x[0], hit_pos_x);
    xs::store(&ray_pack.hit_pos_y[0], hit_pos_y);
    xs::store(&ray_pack.hit_pos_z[0], hit_pos_z);

    
}

void inline return_any_hit_shadow(RP8& ray_pack, const Scene& scene, xsimd::batch_bool<fl>& in_light, int light_index, const f_batch& light_distance){
    
    // --- Step 1: Load and Normalize Ray Packet Data ---
    
    // Load UNNORMALIZED ray directions (from surface to light)
    auto dir_x = xs::load(&ray_pack.d_x[0]);
    auto dir_y = xs::load(&ray_pack.d_y[0]);
    auto dir_z = xs::load(&ray_pack.d_z[0]);

    // Normalize the directions
    auto len_sq = (dir_x * dir_x) + (dir_y * dir_y) + (dir_z * dir_z);
    auto len = xs::sqrt(len_sq);
    b_batch valid_mask = (len > 1e-8f);
    f_batch zero_batch(0.0f);
    dir_x = xs::select(valid_mask, dir_x / len, zero_batch);
    dir_y = xs::select(valid_mask, dir_y / len, zero_batch);
    dir_z = xs::select(valid_mask, dir_z / len, zero_batch);

    // Load ray origins
    f_batch orig_x = xs::load(&ray_pack.o_x[0]);
    f_batch orig_y = xs::load(&ray_pack.o_y[0]);
    f_batch orig_z = xs::load(&ray_pack.o_z[0]);

    // Use the light distance as the maximum t-value.
    // We are only looking for hits *between* 0 and light_distance.
    const f_batch& t_max = light_distance;
    
    // Get scene epsilon for intersection tests
    const fl epsilon = scene.intersection_test_epsilon;

    // This mask tracks which rays are *occluded* (in shadow).
    // We start by assuming all rays are in light (not occluded).
    b_batch is_occluded = 1.0f<0.f;

    const VertexData& vertices = scene.vertex_data__;

    // --- Step 2: Intersect Spheres (Any-Hit) ---
    int num_spheres = scene.sphere_data__.sphere_id.size();
    for (int i = 0; i < num_spheres; ++i) {
        // Find which rays are still active (not yet occluded)
        b_batch active_rays = !is_occluded;
        if (xs::none(active_rays)) break; // All rays are in shadow, we can stop.

        int center_vertex_id = scene.sphere_data__.sphere_center_vertex_id[i];
        auto sphere_center_x = xs::broadcast(vertices.v_pos_x[center_vertex_id]);
        auto sphere_center_y = xs::broadcast(vertices.v_pos_y[center_vertex_id]);
        auto sphere_center_z = xs::broadcast(vertices.v_pos_z[center_vertex_id]);
        auto sphere_radius_sq = xs::broadcast(scene.sphere_data__.sphere_radius_sq[i]);

        auto oc_x = orig_x - sphere_center_x;
        auto oc_y = orig_y - sphere_center_y;
        auto oc_z = orig_z - sphere_center_z;

        auto b_half = (dir_x * oc_x) + (dir_y * oc_y) + (dir_z * oc_z);
        auto c = (oc_x * oc_x) + (oc_y * oc_y) + (oc_z * oc_z) - sphere_radius_sq;
        auto discriminant = (b_half * b_half) - c;

        b_batch hit_mask = (discriminant >= 0.f);
        if (xs::none(hit_mask)) continue;

        auto sqrt_discriminant = xs::sqrt(discriminant);
        auto t0 = -b_half - sqrt_discriminant;
        auto t1 = -b_half + sqrt_discriminant;

        auto t_smaller = xs::min(t0, t1);
        auto t_larger = xs::max(t0, t1);

        auto t_new = xs::select(t_smaller > epsilon, t_smaller, t_larger);
        
        // A ray is newly occluded if it's active AND it hit
        // AND the hit is in front of the ray (t > epsilon)
        // AND the hit is *before* the light (t < t_max).
        b_batch new_occlusion = active_rays & hit_mask & (t_new > epsilon) & (t_new < t_max);

        is_occluded = is_occluded | new_occlusion;
    }

    // --- Step 3: Intersect Triangles (Any-Hit, No Backface Culling) ---
    int num_triangles = scene.triangle_data__.triangle_id.size();
    for (int i = 0; i < num_triangles; ++i) {
        b_batch active_rays = !is_occluded;
        if (xs::none(active_rays)) break;

        // Get vertex indices
        int i0 = scene.triangle_data__.v0_ind[i];
        int i1 = scene.triangle_data__.v1_ind[i];
        int i2 = scene.triangle_data__.v2_ind[i];

        // Get vertex positions
        Vec3f a = { vertices.v_pos_x[i0], vertices.v_pos_y[i0], vertices.v_pos_z[i0] };
        Vec3f b = { vertices.v_pos_x[i1], vertices.v_pos_y[i1], vertices.v_pos_z[i1] };
        Vec3f c = { vertices.v_pos_x[i2], vertices.v_pos_y[i2], vertices.v_pos_z[i2] };

        // Calculate edges
        auto e1_x = xs::broadcast(b.x - a.x);
        auto e1_y = xs::broadcast(b.y - a.y);
        auto e1_z = xs::broadcast(b.z - a.z);
        auto e2_x = xs::broadcast(c.x - a.x);
        auto e2_y = xs::broadcast(c.y - a.y);
        auto e2_z = xs::broadcast(c.z - a.z);

        // pvec = cross(dir, e2)
        auto pvec_x = dir_y * e2_z - dir_z * e2_y;
        auto pvec_y = dir_z * e2_x - dir_x * e2_z;
        auto pvec_z = dir_x * e2_y - dir_y * e2_x;

        // det = dot(e1, pvec)
        auto det = e1_x * pvec_x + e1_y * pvec_y + e1_z * pvec_z;

        // ##################################################
        // ## NO BACKFACE CULLING ##
        // Check determinant using abs() to allow all hits
        b_batch det_mask = (xs::abs(det) > epsilon);
        // ##################################################

        if (xs::none(det_mask)) continue;

        auto invDet = xs::select(det_mask, 1.0f / det, zero_batch);

        // tvec = orig - a
        auto a_x = xs::broadcast(a.x);
        auto a_y = xs::broadcast(a.y);
        auto a_z = xs::broadcast(a.z);
        auto tvec_x = orig_x - a_x;
        auto tvec_y = orig_y - a_y;
        auto tvec_z = orig_z - a_z;

        // u = dot(tvec, pvec) * invDet
        auto u = (tvec_x * pvec_x + tvec_y * pvec_y + tvec_z * pvec_z) * invDet;

        const fl bary_epsilon = 1e-5f;
        f_batch one_batch(1.0f);
        b_batch bary_mask = (u >= zero_batch - bary_epsilon) & (u <= one_batch + bary_epsilon);
        if (xs::none(bary_mask)) continue;

        // qvec = cross(tvec, e1)
        auto qvec_x = tvec_y * e1_z - tvec_z * e1_y;
        auto qvec_y = tvec_z * e1_x - tvec_x * e1_z;
        auto qvec_z = tvec_x * e1_y - tvec_y * e1_x;

        // v = dot(dir, qvec) * invDet
        auto v = (dir_x * qvec_x + dir_y * qvec_y + dir_z * qvec_z) * invDet;

        bary_mask = bary_mask & (v >= zero_batch - bary_epsilon) & (u + v <= one_batch + bary_epsilon);
        if (xs::none(bary_mask)) continue;

        // t = dot(e2, qvec) * invDet
        auto t_new = (e2_x * qvec_x + e2_y * qvec_y + e2_z * qvec_z) * invDet;

        // A ray is newly occluded if it's active AND passes all tests
        // AND the hit is in front of the ray (t > epsilon)
        // AND the hit is *before* the light (t < t_max).
        b_batch new_occlusion = active_rays & det_mask & bary_mask & (t_new > epsilon) & (t_new < t_max);

        is_occluded = is_occluded | new_occlusion;
    }

    // --- Step 4: Intersect Planes (Any-Hit) ---
    int num_planes = scene.plane_data__.plane_id.size();
    for (int i = 0; i < num_planes; ++i) {
        b_batch active_rays = !is_occluded;
        if (xs::none(active_rays)) break;

        auto norm_x = xs::broadcast(scene.plane_data__.plane_norm_x[i]);
        auto norm_y = xs::broadcast(scene.plane_data__.plane_norm_y[i]);
        auto norm_z = xs::broadcast(scene.plane_data__.plane_norm_z[i]);

        int p_v_id = scene.plane_data__.plane_point_vertex_id[i];
        auto p0_x = xs::broadcast(vertices.v_pos_x[p_v_id]);
        auto p0_y = xs::broadcast(vertices.v_pos_y[p_v_id]);
        auto p0_z = xs::broadcast(vertices.v_pos_z[p_v_id]);

        auto denom = (dir_x * norm_x) + (dir_y * norm_y) + (dir_z * norm_z);
        b_batch parallel_mask = (xs::abs(denom) > epsilon);
        if (xs::none(parallel_mask)) continue;

        auto p0_minus_o_x = p0_x - orig_x;
        auto p0_minus_o_y = p0_y - orig_y;
        auto p0_minus_o_z = p0_z - orig_z;
        auto numer = (p0_minus_o_x * norm_x) + (p0_minus_o_y * norm_y) + (p0_minus_o_z * norm_z);
        auto t_new = numer / denom;
        
        // A ray is newly occluded if it's active AND not parallel
        // AND the hit is in front of the ray (t > epsilon)
        // AND the hit is *before* the light (t < t_max).
        b_batch new_occlusion = active_rays & parallel_mask & (t_new > epsilon) & (t_new < t_max);

        is_occluded = is_occluded | new_occlusion;
    }

    // --- Step 5: Set Final Output ---
    // 'in_light' is the inverse of 'is_occluded'
    in_light = !is_occluded;
}


//MASKED CODE 
// --- MASKED SIMD FUNCTIONS ---
// These functions operate only on the lanes specified by the active_mask.

void inline intersect_planes_masked(const f_batch& dir_x,
    const f_batch& dir_y,
    const f_batch& dir_z,
    const f_batch& orig_x,
    const f_batch& orig_y,
    const f_batch& orig_z,
    f_batch& t_min,
    i_batch& mat_id,
    f_batch& norm_x_out,
    f_batch& norm_y_out,
    f_batch& norm_z_out,
    const PlaneData& planes, const VertexData& vertices, const Scene& scene, int plane_ID,
    const b_batch& active_mask) { // New parameter

    // 1. Early exit if no rays in this packet are active
    if (xs::none(active_mask)) {
        return;
    }

    // 2. Get plane data
    auto norm_x = xs::broadcast(planes.plane_norm_x[plane_ID]);
    auto norm_y = xs::broadcast(planes.plane_norm_y[plane_ID]);
    auto norm_z = xs::broadcast(planes.plane_norm_z[plane_ID]);

    int p_v_id = planes.plane_point_vertex_id[plane_ID];
    auto p0_x = xs::broadcast(vertices.v_pos_x[p_v_id]);
    auto p0_y = xs::broadcast(vertices.v_pos_y[p_v_id]);
    auto p0_z = xs::broadcast(vertices.v_pos_z[p_v_id]);

    // 3. Calculate denominator
    auto denom = (dir_x * norm_x) + (dir_y * norm_y) + (dir_z * norm_z);

    // 4. Epsilon check: combine incoming mask with parallel check
    const fl epsilon = scene.intersection_test_epsilon;
    b_batch parallel_mask = (xs::abs(denom) > epsilon);
    b_batch plane_active_mask = active_mask & parallel_mask; // Only test active, non-parallel rays

    if (xs::none(plane_active_mask)) {
        return; // All active rays are parallel to the plane
    }

    // 5. Calculate numerator
    auto p0_minus_o_x = p0_x - orig_x;
    auto p0_minus_o_y = p0_y - orig_y;
    auto p0_minus_o_z = p0_z - orig_z;

    auto numer = (p0_minus_o_x * norm_x) + (p0_minus_o_y * norm_y) + (p0_minus_o_z * norm_z);

    // 6. Calculate t
    auto t_new = numer / denom;

    // 7. Create final hit mask (must be a subset of plane_active_mask)
    b_batch final_hit_mask = plane_active_mask & (t_new > epsilon) & (t_new < t_min);

    if (xs::none(final_hit_mask)) {
        return; // No new valid hits
    }

    // --- Update Hit Information ---

    // 8. Update normals
    norm_x_out = xs::select(final_hit_mask, norm_x, norm_x_out);
    norm_y_out = xs::select(final_hit_mask, norm_y, norm_y_out);
    norm_z_out = xs::select(final_hit_mask, norm_z, norm_z_out);

    // 9. Update t_min
    t_min = xs::select(final_hit_mask, t_new, t_min);

    // 10. Update material ID
    auto int_mask = xs::batch_bool_cast<int32_t>(final_hit_mask);
    i_batch plane_mat_id_batch = xs::broadcast(planes.plane_material_id[plane_ID]);
    mat_id = xs::select(int_mask, plane_mat_id_batch, mat_id);
}

void inline intersect_spheres_masked(
    const f_batch& dir_x,
    const f_batch& dir_y,
    const f_batch& dir_z,
    const f_batch& orig_x,
    const f_batch& orig_y,
    const f_batch& orig_z,
    f_batch& t_min,
    i_batch& mat_id,
    f_batch& norm_x_out,
    f_batch& norm_y_out,
    f_batch& norm_z_out,
    const SphereData& spheres, const VertexData& vertices, const Scene& scene, const int Sphere_ID,
    const b_batch& active_mask) { // New parameter

    // 1. Early exit if no rays in this packet are active
    if (xs::none(active_mask)) {
        return;
    }

    int center_vertex_id = spheres.sphere_center_vertex_id[Sphere_ID];
    auto sphere_center_x = xs::broadcast(vertices.v_pos_x[center_vertex_id]);
    auto sphere_center_y = xs::broadcast(vertices.v_pos_y[center_vertex_id]);
    auto sphere_center_z = xs::broadcast(vertices.v_pos_z[center_vertex_id]);
    auto sphere_radius_sq = xs::broadcast(spheres.sphere_radius_sq[Sphere_ID]);

    auto oc_x = orig_x - sphere_center_x;
    auto oc_y = orig_y - sphere_center_y;
    auto oc_z = orig_z - sphere_center_z;

    auto b_half = (dir_x * oc_x) + (dir_y * oc_y) + (dir_z * oc_z);

    auto c = (oc_x * oc_x) + (oc_y * oc_y) + (oc_z * oc_z) - sphere_radius_sq;

    auto discriminant = (b_half * b_half) - c;

    // 2. Combine incoming mask with discriminant check
    b_batch discriminant_mask = (discriminant >= 0.f);
    b_batch sphere_active_mask = active_mask & discriminant_mask;

    if (xs::none(sphere_active_mask)) {
        return;
    }

    auto sqrt_discriminant = xs::sqrt(discriminant);

    auto t0 = -b_half - sqrt_discriminant;
    auto t1 = -b_half + sqrt_discriminant;

    auto t_smaller = xs::min(t0, t1);
    auto t_larger = xs::max(t0, t1);

    const fl epsilon = scene.intersection_test_epsilon;

    auto t_new = xs::select(t_smaller > epsilon, t_smaller, t_larger);

    // 3. Create final hit mask (must be a subset of sphere_active_mask)
    b_batch final_hit_mask = sphere_active_mask & (t_new > epsilon) & (t_new < t_min);

    if (xs::none(final_hit_mask)) {
        return;
    }

    // --- Update Hit Information ---
    
    // 1. Calculate hit point: P = O + t*D
    auto hit_x = orig_x + t_new * dir_x;
    auto hit_y = orig_y + t_new * dir_y;
    auto hit_z = orig_z + t_new * dir_z;

    // 2. Calculate normal: N = (P - C) / R
    auto n_x = hit_x - sphere_center_x;
    auto n_y = hit_y - sphere_center_y;
    auto n_z = hit_z - sphere_center_z;

    // 3. Normalize
    auto radius = xs::sqrt(sphere_radius_sq);
    f_batch epsilon_batch(1e-8f);
    b_batch radius_valid = (radius > epsilon_batch);
    auto inv_radius = xs::select(radius_valid, 1.0f / radius, f_batch(0.0f));
    
    auto new_norm_x = n_x * inv_radius;
    auto new_norm_y = n_y * inv_radius;
    auto new_norm_z = n_z * inv_radius;

    // 4. Conditionally store the new normal
    norm_x_out = xs::select(final_hit_mask, new_norm_x, norm_x_out);
    norm_y_out = xs::select(final_hit_mask, new_norm_y, norm_y_out);
    norm_z_out = xs::select(final_hit_mask, new_norm_z, norm_z_out);

    // 5. Update t_min
    t_min = xs::select(final_hit_mask, t_new, t_min);

    // 6. Update material ID
    auto int_mask = xs::batch_bool_cast<int32_t>(final_hit_mask);
    i_batch sphere_id_batch = xs::broadcast(spheres.sphere_mat_id[Sphere_ID]);
    mat_id = xs::select(int_mask, sphere_id_batch, mat_id);
}

void inline intersect_triangles_masked(
    const f_batch& dir_x,
    const f_batch& dir_y,
    const f_batch& dir_z,
    const f_batch& orig_x,
    const f_batch& orig_y,
    const f_batch& orig_z,
    f_batch& t_min,
    i_batch& mat_id,
    f_batch& norm_x_out,
    f_batch& norm_y_out,
    f_batch& norm_z_out,
    const TriangleData& triangles, const VertexData& vertices, const Scene& scene, const int Triangle_ID,
    const b_batch& active_mask) { // New parameter

    // 1. Early exit if no rays in this packet are active
    if (xs::none(active_mask)) {
        return;
    }

    // Get vertex indices
    int i0 = triangles.v0_ind[Triangle_ID];
    int i1 = triangles.v1_ind[Triangle_ID];
    int i2 = triangles.v2_ind[Triangle_ID];

    // Get triangle normal for backface culling
    auto tri_norm_x = xs::broadcast(triangles.tri_norm_x[Triangle_ID]);
    auto tri_norm_y = xs::broadcast(triangles.tri_norm_y[Triangle_ID]);
    auto tri_norm_z = xs::broadcast(triangles.tri_norm_z[Triangle_ID]);
    
    // 2. Backface culling: combine incoming mask with backface check
    auto dot_normal_ray = (tri_norm_x * dir_x) + (tri_norm_y * dir_y) + (tri_norm_z * dir_z);
    b_batch backface_mask = (dot_normal_ray < 0.0f);
    b_batch tri_active_mask = active_mask & backface_mask;
    
    if (xs::none(tri_active_mask)) {
        return; // All active rays are hitting backfaces
    }

    // Get vertex positions
    Vec3f a = { vertices.v_pos_x[i0], vertices.v_pos_y[i0], vertices.v_pos_z[i0] };
    Vec3f b = { vertices.v_pos_x[i1], vertices.v_pos_y[i1], vertices.v_pos_z[i1] };
    Vec3f c = { vertices.v_pos_x[i2], vertices.v_pos_y[i2], vertices.v_pos_z[i2] };

    // Calculate edges
    Vec3f e1_scalar = b - a;
    Vec3f e2_scalar = c - a;

    // Broadcast edges
    auto e1_x = xs::broadcast(e1_scalar.x);
    auto e1_y = xs::broadcast(e1_scalar.y);
    auto e1_z = xs::broadcast(e1_scalar.z);
    auto e2_x = xs::broadcast(e2_scalar.x);
    auto e2_y = xs::broadcast(e2_scalar.y);
    auto e2_z = xs::broadcast(e2_scalar.z);

    // pvec = cross(dir, e2)
    auto pvec_x = dir_y * e2_z - dir_z * e2_y;
    auto pvec_y = dir_z * e2_x - dir_x * e2_z;
    auto pvec_z = dir_x * e2_y - dir_y * e2_x;

    // det = dot(e1, pvec)
    auto det = e1_x * pvec_x + e1_y * pvec_y + e1_z * pvec_z;

    // Epsilon checks
    const fl epsilon = 1e-12f;
    const fl bary_epsilon = 1e-5f;
    f_batch zero_batch(0.0f);
    f_batch one_batch(1.0f);

    // 3. Check determinant (and combine with active mask)
    b_batch det_mask = (det > epsilon);
    tri_active_mask = tri_active_mask & det_mask;

    if (xs::none(tri_active_mask)) {
        return;
    }

    // Safe division
    auto invDet = xs::select(tri_active_mask, 1.0f / det, zero_batch); // Correct use of zero_batch

    // tvec = orig - a
    auto a_x = xs::broadcast(a.x);
    auto a_y = xs::broadcast(a.y);
    auto a_z = xs::broadcast(a.z);
    auto tvec_x = orig_x - a_x;
    auto tvec_y = orig_y - a_y;
    auto tvec_z = orig_z - a_z;

    // u = dot(tvec, pvec) * invDet
    auto u = (tvec_x * pvec_x + tvec_y * pvec_y + tvec_z * pvec_z) * invDet;

    // 4. Mask out rays with invalid barycentric coord u
    tri_active_mask = tri_active_mask & (u >= zero_batch - bary_epsilon) & (u <= one_batch + bary_epsilon);
    if (xs::none(tri_active_mask)) return;

    // qvec = cross(tvec, e1)
    auto qvec_x = tvec_y * e1_z - tvec_z * e1_y;
    auto qvec_y = tvec_z * e1_x - tvec_x * e1_z;
    auto qvec_z = tvec_x * e1_y - tvec_y * e1_x;

    // v = dot(dir, qvec) * invDet
    auto v = (dir_x * qvec_x + dir_y * qvec_y + dir_z * qvec_z) * invDet;

    // 5. Mask out rays with invalid barycentric coord v
    tri_active_mask = tri_active_mask & (v >= zero_batch - bary_epsilon) & (u + v <= one_batch + bary_epsilon);
    if (xs::none(tri_active_mask)) return;

    // t = dot(e2, qvec) * invDet
    auto t_new = (e2_x * qvec_x + e2_y * qvec_y + e2_z * qvec_z) * invDet;

    // 6. Final mask: must be a subset of tri_active_mask
    b_batch final_hit_mask = tri_active_mask & (t_new > epsilon) & (t_new < t_min);

    if (xs::none(final_hit_mask)) {
        return;
    }

    // --- Update Hit Information ---
    
    // 1. Get the precomputed triangle normal
    auto new_norm_x = tri_norm_x;
    auto new_norm_y = tri_norm_y;
    auto new_norm_z = tri_norm_z;
    
    // 2. Conditionally store the triangle normal
    norm_x_out = xs::select(final_hit_mask, new_norm_x, norm_x_out);
    norm_y_out = xs::select(final_hit_mask, new_norm_y, norm_y_out);
    norm_z_out = xs::select(final_hit_mask, new_norm_z, norm_z_out);

    // 3. Update t_min
    t_min = xs::select(final_hit_mask, t_new, t_min);

    // 4. Update material ID
    auto int_mask = xs::batch_bool_cast<int32_t>(final_hit_mask);
    auto mat_id_batch = xs::broadcast(triangles.triangle_material_id[Triangle_ID]);
    mat_id = xs::select(int_mask, mat_id_batch, mat_id);
}

void return_closest_hit_masked(RP8& ray_pack, const Scene& scene, const VertexData& vertices,
                               const b_batch& active_mask) { // New parameter
    
    // 1. Early exit if no rays are active
    if (xs::none(active_mask)) {
        return;
    }

    // 2. Load the UNNORMALIZED data
    auto dir_x = xs::load(&ray_pack.d_x[0]);
    auto dir_y = xs::load(&ray_pack.d_y[0]);
    auto dir_z = xs::load(&ray_pack.d_z[0]);

    // 3. Calculate length-squared
    auto len_sq = (dir_x * dir_x) + (dir_y * dir_y) + (dir_z * dir_z);

    // 4. Calculate length      
    auto len = xs::sqrt(len_sq);

    // 5. Create a mask to prevent division by zero and normalize (only for active rays)
    b_batch valid_len_mask = (len > 1e-8f);
    b_batch norm_mask = active_mask & valid_len_mask; // Only normalize active rays with valid length

    f_batch zero_batch(0.0f);
    dir_x = xs::select(norm_mask, dir_x / len, zero_batch); // Correct use of zero_batch
    dir_y = xs::select(norm_mask, dir_y / len, zero_batch); // Correct use of zero_batch
    dir_z = xs::select(norm_mask, dir_z / len, zero_batch); // Correct use of zero_batch

    // 6. Load ray origins
    f_batch orig_x = xs::load(&ray_pack.o_x[0]);
    f_batch orig_y = xs::load(&ray_pack.o_y[0]);
    f_batch orig_z = xs::load(&ray_pack.o_z[0]);

    // 7. Load current intersection state
    f_batch t_min = xs::load(&ray_pack.t_min[0]);
    i_batch hit_id = xs::load(&ray_pack.mat_id[0]);
    f_batch norm_x = xs::load(&ray_pack.hit_norm_x[0]);
    f_batch norm_y = xs::load(&ray_pack.hit_norm_y[0]);
    f_batch norm_z = xs::load(&ray_pack.hit_norm_z[0]);


    // --- Iterate over all scene objects (SoA) ---

    // 8. Intersect Spheres
    int num_spheres = scene.sphere_data__.sphere_id.size();
    for (int i = 0; i < num_spheres; ++i) {
        intersect_spheres_masked(dir_x, dir_y, dir_z, orig_x, orig_y, orig_z,
            t_min, hit_id, norm_x, norm_y, norm_z,
            scene.sphere_data__, vertices, scene, i,
            active_mask); // Pass the mask
    }

    // 9. Intersect Triangles
    int num_triangles = scene.triangle_data__.triangle_id.size();
    for (int i = 0; i < num_triangles; ++i) {
        intersect_triangles_masked(dir_x, dir_y, dir_z, orig_x, orig_y, orig_z,
            t_min, hit_id, norm_x, norm_y, norm_z,
            scene.triangle_data__, vertices, scene, i,
            active_mask); // Pass the mask
    }

    // 10. Intersect Planes
    int num_planes = scene.plane_data__.plane_id.size();
    for (int i = 0; i < num_planes; ++i) {
        intersect_planes_masked(dir_x, dir_y, dir_z, orig_x, orig_y, orig_z,
            t_min, hit_id, norm_x, norm_y, norm_z,
            scene.plane_data__, vertices, scene, i,
            active_mask); // Pass the mask
    }
    

    // 11. Store the final results back into the ray packet
    xs::store(&ray_pack.t_min[0], t_min);
    xs::store(&ray_pack.mat_id[0], hit_id);
    xs::store(&ray_pack.hit_norm_x[0], norm_x);
    xs::store(&ray_pack.hit_norm_y[0], norm_y);
    xs::store(&ray_pack.hit_norm_z[0], norm_z);

    // 12. Calculate and store hit positions: P = O + t * D
    f_batch infinity_batch = xs::broadcast(__builtin_inff());
    // Only calculate for active rays that *actually* hit something
    b_batch final_hit_mask = (t_min < infinity_batch) & active_mask; 

    // Use select to avoid (inf * 0 = NaN) for rays that missed
    auto hit_pos_x = xs::select(final_hit_mask, orig_x + t_min * dir_x, zero_batch);
    auto hit_pos_y = xs::select(final_hit_mask, orig_y + t_min * dir_y, zero_batch);
    auto hit_pos_z = xs::select(final_hit_mask, orig_z + t_min * dir_z, zero_batch);
    
    xs::store(&ray_pack.hit_pos_x[0], hit_pos_x);
    xs::store(&ray_pack.hit_pos_y[0], hit_pos_y);
    xs::store(&ray_pack.hit_pos_z[0], hit_pos_z);
}

void inline return_any_hit_shadow_masked(RP8& ray_pack, const Scene& scene, xsimd::batch_bool<fl>& in_light, int light_index, const f_batch& light_distance,
                                         const b_batch& active_mask) { // New parameter
    
    // 1. Early exit if no rays are active
    if (xs::none(active_mask)) {
        in_light = false; // Set all to "not occluded"
        return;
    }
    
    // --- Step 1: Load and Normalize Ray Packet Data ---
    
    auto dir_x = xs::load(&ray_pack.d_x[0]);
    auto dir_y = xs::load(&ray_pack.d_y[0]);
    auto dir_z = xs::load(&ray_pack.d_z[0]);

    // Normalize the directions (only for active rays)
    auto len_sq = (dir_x * dir_x) + (dir_y * dir_y) + (dir_z * dir_z);
    auto len = xs::sqrt(len_sq);
    b_batch valid_len_mask = (len > 1e-8f);
    b_batch norm_mask = active_mask & valid_len_mask;

    f_batch zero_batch(0.0f);
    dir_x = xs::select(norm_mask, dir_x / len, zero_batch);
    dir_y = xs::select(norm_mask, dir_y / len, zero_batch);
    dir_z = xs::select(norm_mask, dir_z / len, zero_batch);

    f_batch orig_x = xs::load(&ray_pack.o_x[0]);
    f_batch orig_y = xs::load(&ray_pack.o_y[0]);
    f_batch orig_z = xs::load(&ray_pack.o_z[0]);

    const f_batch& t_max = light_distance;
    const fl epsilon = scene.intersection_test_epsilon;

    // This mask tracks which rays are *occluded* (in shadow).
    b_batch is_occluded = false; // All false

    const VertexData& vertices = scene.vertex_data__;

    // --- Step 2: Intersect Spheres (Any-Hit) ---
    int num_spheres = scene.sphere_data__.sphere_id.size();
    for (int i = 0; i < num_spheres; ++i) {
        // Find which rays are still active AND not yet occluded
        b_batch sphere_active_mask = !is_occluded & active_mask;
        if (xs::none(sphere_active_mask)) break;

        int center_vertex_id = scene.sphere_data__.sphere_center_vertex_id[i];
        auto sphere_center_x = xs::broadcast(vertices.v_pos_x[center_vertex_id]);
        auto sphere_center_y = xs::broadcast(vertices.v_pos_y[center_vertex_id]);
        auto sphere_center_z = xs::broadcast(vertices.v_pos_z[center_vertex_id]);
        auto sphere_radius_sq = xs::broadcast(scene.sphere_data__.sphere_radius_sq[i]);

        auto oc_x = orig_x - sphere_center_x;
        auto oc_y = orig_y - sphere_center_y;
        auto oc_z = orig_z - sphere_center_z;

        auto b_half = (dir_x * oc_x) + (dir_y * oc_y) + (dir_z * oc_z);
        auto c = (oc_x * oc_x) + (oc_y * oc_y) + (oc_z * oc_z) - sphere_radius_sq;
        auto discriminant = (b_half * b_half) - c;

        b_batch hit_mask = (discriminant >= 0.f);
        if (xs::none(hit_mask)) continue;

        auto sqrt_discriminant = xs::sqrt(discriminant);
        auto t0 = -b_half - sqrt_discriminant;
        auto t1 = -b_half + sqrt_discriminant;

        auto t_smaller = xs::min(t0, t1);
        auto t_larger = xs::max(t0, t1);

        auto t_new = xs::select(t_smaller > epsilon, t_smaller, t_larger);
        
        b_batch new_occlusion = sphere_active_mask & hit_mask & (t_new > epsilon) & (t_new < t_max);

        is_occluded = is_occluded | new_occlusion;
    }

    // --- Step 3: Intersect Triangles (Any-Hit, No Backface Culling) ---
    int num_triangles = scene.triangle_data__.triangle_id.size();
    for (int i = 0; i < num_triangles; ++i) {
        b_batch tri_active_mask = !is_occluded & active_mask;
        if (xs::none(tri_active_mask)) break;

        // Get vertex indices
        int i0 = scene.triangle_data__.v0_ind[i];
        int i1 = scene.triangle_data__.v1_ind[i];
        int i2 = scene.triangle_data__.v2_ind[i];

        // Get vertex positions
        Vec3f a = { vertices.v_pos_x[i0], vertices.v_pos_y[i0], vertices.v_pos_z[i0] };
        Vec3f b = { vertices.v_pos_x[i1], vertices.v_pos_y[i1], vertices.v_pos_z[i1] };
        Vec3f c = { vertices.v_pos_x[i2], vertices.v_pos_y[i2], vertices.v_pos_z[i2] };

        // ... (Omitting the rest of the MT algorithm for brevity, it's identical to legacy) ...
        auto e1_x = xs::broadcast(b.x - a.x);
        auto e1_y = xs::broadcast(b.y - a.y);
        auto e1_z = xs::broadcast(b.z - a.z);
        auto e2_x = xs::broadcast(c.x - a.x);
        auto e2_y = xs::broadcast(c.y - a.y);
        auto e2_z = xs::broadcast(c.z - a.z);
        auto pvec_x = dir_y * e2_z - dir_z * e2_y;
        auto pvec_y = dir_z * e2_x - dir_x * e2_z;
        auto pvec_z = dir_x * e2_y - dir_y * e2_x;
        auto det = e1_x * pvec_x + e1_y * pvec_y + e1_z * pvec_z;
        b_batch det_mask = (xs::abs(det) > epsilon);
        if (xs::none(det_mask)) continue;
        auto invDet = xs::select(det_mask, 1.0f / det, zero_batch);
        auto a_x = xs::broadcast(a.x);
        auto a_y = xs::broadcast(a.y);
        auto a_z = xs::broadcast(a.z);
        auto tvec_x = orig_x - a_x;
        auto tvec_y = orig_y - a_y;
        auto tvec_z = orig_z - a_z;
        auto u = (tvec_x * pvec_x + tvec_y * pvec_y + tvec_z * pvec_z) * invDet;
        const fl bary_epsilon = 1e-5f;
        f_batch one_batch(1.0f);
        b_batch bary_mask = (u >= zero_batch - bary_epsilon) & (u <= one_batch + bary_epsilon);
        if (xs::none(bary_mask)) continue;
        auto qvec_x = tvec_y * e1_z - tvec_z * e1_y;
        auto qvec_y = tvec_z * e1_x - tvec_x * e1_z;
        auto qvec_z = tvec_x * e1_y - tvec_y * e1_x;
        auto v = (dir_x * qvec_x + dir_y * qvec_y + dir_z * qvec_z) * invDet;
        bary_mask = bary_mask & (v >= zero_batch - bary_epsilon) & (u + v <= one_batch + bary_epsilon);
        if (xs::none(bary_mask)) continue;
        auto t_new = (e2_x * qvec_x + e2_y * qvec_y + e2_z * qvec_z) * invDet;

        b_batch new_occlusion = tri_active_mask & det_mask & bary_mask & (t_new > epsilon) & (t_new < t_max);

        is_occluded = is_occluded | new_occlusion;
    }

    // --- Step 4: Intersect Planes (Any-Hit) ---
    int num_planes = scene.plane_data__.plane_id.size();
    for (int i = 0; i < num_planes; ++i) {
        b_batch plane_active_mask = !is_occluded & active_mask;
        if (xs::none(plane_active_mask)) break;

        auto norm_x = xs::broadcast(scene.plane_data__.plane_norm_x[i]);
        auto norm_y = xs::broadcast(scene.plane_data__.plane_norm_y[i]);
        auto norm_z = xs::broadcast(scene.plane_data__.plane_norm_z[i]);

        int p_v_id = scene.plane_data__.plane_point_vertex_id[i];
        auto p0_x = xs::broadcast(vertices.v_pos_x[p_v_id]);
        auto p0_y = xs::broadcast(vertices.v_pos_y[p_v_id]);
        auto p0_z = xs::broadcast(vertices.v_pos_z[p_v_id]);

        auto denom = (dir_x * norm_x) + (dir_y * norm_y) + (dir_z * norm_z);
        b_batch parallel_mask = (xs::abs(denom) > epsilon);
        if (xs::none(parallel_mask)) continue;

        auto p0_minus_o_x = p0_x - orig_x;
        auto p0_minus_o_y = p0_y - orig_y;
        auto p0_minus_o_z = p0_z - orig_z;
        auto numer = (p0_minus_o_x * norm_x) + (p0_minus_o_y * norm_y) + (p0_minus_o_z * norm_z);
        auto t_new = numer / denom;
        
        b_batch new_occlusion = plane_active_mask & parallel_mask & (t_new > epsilon) & (t_new < t_max);

        is_occluded = is_occluded | new_occlusion;
    }

    // --- Step 5: Set Final Output ---
    // 'in_light' should be true ONLY for rays that were in the active_mask
    // AND were not occluded.
    in_light = !is_occluded & active_mask;
}





#endif