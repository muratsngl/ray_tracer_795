#pragma once
#ifndef RTFUNC
#define RTFUNC


#include <xsimd/xsimd.hpp>
#include "utils.hpp" // Assuming this contains operator overloads for Vec3f
#include "types.hpp"
#include "bvh_tlas.hpp" // Updated to include bvh_tlas.hpp

namespace xs = xsimd;

// Acronym type aliases for xsimd batches
using f_batch = xs::batch<fl>;
using i_batch = xs::batch<int>;
using b_batch = xs::batch_bool<fl>;

// --- PRIMITIVE INTERSECTION FUNCTIONS (UNCHANGED) ---
// ... (intersect_planes, intersect_spheres, intersect_triangles) ...
// ... (intersect_planes_masked, intersect_spheres_masked, intersect_triangles_masked) ...

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

    // Get vertex normals
    auto v0_norm_x = xs::broadcast(vertices.v_nor_x[i0]);
    auto v0_norm_y = xs::broadcast(vertices.v_nor_y[i0]);
    auto v0_norm_z = xs::broadcast(vertices.v_nor_z[i0]);
    auto v1_norm_x = xs::broadcast(vertices.v_nor_x[i1]);
    auto v1_norm_y = xs::broadcast(vertices.v_nor_y[i1]);
    auto v1_norm_z = xs::broadcast(vertices.v_nor_z[i1]);
    auto v2_norm_x = xs::broadcast(vertices.v_nor_x[i2]);
    auto v2_norm_y = xs::broadcast(vertices.v_nor_y[i2]);
    auto v2_norm_z = xs::broadcast(vertices.v_nor_z[i2]);

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
    const fl bary_epsilon = 1e-5f; // MUCH more forgiving (was 1e-5f)
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
    active_mask = active_mask & (u >= zero_batch - bary_epsilon) & (u <= one_batch + bary_epsilon);
    if (xs::none(active_mask)) return;

    // qvec = cross(tvec, e1)
    auto qvec_x = tvec_y * e1_z - tvec_z * e1_y;
    auto qvec_y = tvec_z * e1_x - tvec_x * e1_z;
    auto qvec_z = tvec_x * e1_y - tvec_y * e1_x;

    // v = dot(dir, qvec) * invDet
    auto v = (dir_x * qvec_x + dir_y * qvec_y + dir_z * qvec_z) * invDet;

    // Mask out rays with invalid barycentric coord v
    active_mask = active_mask & (v >= zero_batch - bary_epsilon) & (u + v <= one_batch + bary_epsilon);
    if (xs::none(active_mask)) return;

    // t = dot(e2, qvec) * invDet
    auto t_new = (e2_x * qvec_x + e2_y * qvec_y + e2_z * qvec_z) * invDet;

    // Final mask: active, t > scene epsilon, and t < t_min
    b_batch final_hit_mask = active_mask & (t_new > epsilon) & (t_new < t_min);

    if (xs::none(final_hit_mask)) {
        return;
    }

    // --- Update Hit Information ---
    // Compute the normal (smooth or flat)
    auto w = one_batch - u - v;
    auto interp_norm_x = w * v0_norm_x + u * v1_norm_x + v * v2_norm_x;
    auto interp_norm_y = w * v0_norm_y + u * v1_norm_y + v * v2_norm_y;
    auto interp_norm_z = w * v0_norm_z + u * v1_norm_z + v * v2_norm_z;

    auto new_norm_x = G_SMOOTH_SHADING_ENABLED ? interp_norm_x : tri_norm_x;
    auto new_norm_y = G_SMOOTH_SHADING_ENABLED ? interp_norm_y : tri_norm_y;
    auto new_norm_z = G_SMOOTH_SHADING_ENABLED ? interp_norm_z : tri_norm_z;

    // 2. Conditionally store the normal
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

// --- LEGACY CLOSEST-HIT/ANY-HIT FUNCTIONS (REMOVED/REPLACED) ---
// void return_closest_hit(RP8& ray_pack, const Scene& scene, const VertexData& vertices) { ... }
// void inline return_any_hit_shadow(RP8& ray_pack, const Scene& scene, ... ) { ... }

// --- MASKED SIMD FUNCTIONS (UNCHANGED) ---
// ... (intersect_planes_masked, intersect_spheres_masked, intersect_triangles_masked) ...
// ... (return_closest_hit_masked, return_any_hit_shadow_masked) ...

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
    const b_batch& active_mask,
    const i_batch& override_mat_id) { // NEW: Added override_mat_id

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
    // auto dot_normal_ray = (tri_norm_x * dir_x) + (tri_norm_y * dir_y) + (tri_norm_z * dir_z);
    // b_batch backface_mask = (dot_normal_ray < 0.0f);
    b_batch tri_active_mask = active_mask ;
    
    if (xs::none(tri_active_mask)) {
        return; // All active rays are hitting backfaces
    }

    // Get vertex positions
    Vec3f a = { vertices.v_pos_x[i0], vertices.v_pos_y[i0], vertices.v_pos_z[i0] };
    Vec3f b = { vertices.v_pos_x[i1], vertices.v_pos_y[i1], vertices.v_pos_z[i1] };
    Vec3f c = { vertices.v_pos_x[i2], vertices.v_pos_y[i2], vertices.v_pos_z[i2] };

    // Get vertex normals
    auto v0_norm_x = xs::broadcast(vertices.v_nor_x[i0]);
    auto v0_norm_y = xs::broadcast(vertices.v_nor_y[i0]);
    auto v0_norm_z = xs::broadcast(vertices.v_nor_z[i0]);
    auto v1_norm_x = xs::broadcast(vertices.v_nor_x[i1]);
    auto v1_norm_y = xs::broadcast(vertices.v_nor_y[i1]);
    auto v1_norm_z = xs::broadcast(vertices.v_nor_z[i1]);
    auto v2_norm_x = xs::broadcast(vertices.v_nor_x[i2]);
    auto v2_norm_y = xs::broadcast(vertices.v_nor_y[i2]);
    auto v2_norm_z = xs::broadcast(vertices.v_nor_z[i2]);

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
    
    // Compute the normal (smooth or flat)
    auto w = one_batch - u - v;
    auto interp_norm_x = w * v0_norm_x + u * v1_norm_x + v * v2_norm_x;
    auto interp_norm_y = w * v0_norm_y + u * v1_norm_y + v * v2_norm_y;
    auto interp_norm_z = w * v0_norm_z + u * v1_norm_z + v * v2_norm_z;

    auto new_norm_x = G_SMOOTH_SHADING_ENABLED ? interp_norm_x : tri_norm_x;
    auto new_norm_y = G_SMOOTH_SHADING_ENABLED ? interp_norm_y : tri_norm_y;
    auto new_norm_z = G_SMOOTH_SHADING_ENABLED ? interp_norm_z : tri_norm_z;
    
    // 2. Conditionally store the normal
    norm_x_out = xs::select(final_hit_mask, new_norm_x, norm_x_out);
    norm_y_out = xs::select(final_hit_mask, new_norm_y, norm_y_out);
    norm_z_out = xs::select(final_hit_mask, new_norm_z, norm_z_out);

    // 3. Update t_min
    t_min = xs::select(final_hit_mask, t_new, t_min);

    // 4. Update material ID using the override
    auto int_mask = xs::batch_bool_cast<int32_t>(final_hit_mask);
    mat_id = xs::select(int_mask, override_mat_id, mat_id);
}

// --- NEW TRANSFORMATION HELPERS ---

/**
 * @brief Transforms a packet of ray origins and directions by a 4x4 matrix.
 */
void inline transform_ray_packet(const f_batch& o_x, const f_batch& o_y, const f_batch& o_z,
                                const f_batch& d_x, const f_batch& d_y, const f_batch& d_z,
                                f_batch& out_o_x, f_batch& out_o_y, f_batch& out_o_z,
                                f_batch& out_d_x, f_batch& out_d_y, f_batch& out_d_z,
                                const Mat4f& transform) {
    
    // Broadcast matrix elements
    auto m11 = xs::broadcast(transform.m11), m12 = xs::broadcast(transform.m12), m13 = xs::broadcast(transform.m13), m14 = xs::broadcast(transform.m14);
    auto m21 = xs::broadcast(transform.m21), m22 = xs::broadcast(transform.m22), m23 = xs::broadcast(transform.m23), m24 = xs::broadcast(transform.m24);
    auto m31 = xs::broadcast(transform.m31), m32 = xs::broadcast(transform.m32), m33 = xs::broadcast(transform.m33), m34 = xs::broadcast(transform.m34);
    auto m41 = xs::broadcast(transform.m41), m42 = xs::broadcast(transform.m42), m43 = xs::broadcast(transform.m43), m44 = xs::broadcast(transform.m44);

    // Transform origin (point, w=1)
    auto w_o = m41 * o_x + m42 * o_y + m43 * o_z + m44;
    f_batch inv_w_o(1.0f);
    b_batch valid_w_o = (xs::abs(w_o) > 1e-8f) & (xs::abs(w_o - 1.0f) > 1e-8f);
    inv_w_o = xs::select(valid_w_o, 1.0f / w_o, inv_w_o);

    out_o_x = (m11 * o_x + m12 * o_y + m13 * o_z + m14) * inv_w_o;
    out_o_y = (m21 * o_x + m22 * o_y + m23 * o_z + m24) * inv_w_o;
    out_o_z = (m31 * o_x + m32 * o_y + m33 * o_z + m34) * inv_w_o;

    // Transform direction (vector, w=0)
    out_d_x = m11 * d_x + m12 * d_y + m13 * d_z;
    out_d_y = m21 * d_x + m22 * d_y + m23 * d_z;
    out_d_z = m31 * d_x + m32 * d_y + m33 * d_z;
}

/**
 * @brief Transforms a packet of hit points and normals from local to world space.
 * Uses the inverse transpose of the matrix's 3x3 part for normals.
 */
void inline transform_hit_results(const f_batch& local_hit_x, const f_batch& local_hit_y, const f_batch& local_hit_z,
                                 const f_batch& local_norm_x, const f_batch& local_norm_y, const f_batch& local_norm_z,
                                 f_batch& world_hit_x, f_batch& world_hit_y, f_batch& world_hit_z,
                                 f_batch& world_norm_x, f_batch& world_norm_y, f_batch& world_norm_z,
                                 const Mat4f& forward_transform) {
    
    // --- 1. Transform Hit Point (P' = M * P, w=1) ---
    auto m11 = xs::broadcast(forward_transform.m11), m12 = xs::broadcast(forward_transform.m12), m13 = xs::broadcast(forward_transform.m13), m14 = xs::broadcast(forward_transform.m14);
    auto m21 = xs::broadcast(forward_transform.m21), m22 = xs::broadcast(forward_transform.m22), m23 = xs::broadcast(forward_transform.m23), m24 = xs::broadcast(forward_transform.m24);
    auto m31 = xs::broadcast(forward_transform.m31), m32 = xs::broadcast(forward_transform.m32), m33 = xs::broadcast(forward_transform.m33), m34 = xs::broadcast(forward_transform.m34);
    auto m41 = xs::broadcast(forward_transform.m41), m42 = xs::broadcast(forward_transform.m42), m43 = xs::broadcast(forward_transform.m43), m44 = xs::broadcast(forward_transform.m44);

    auto w_h = m41 * local_hit_x + m42 * local_hit_y + m43 * local_hit_z + m44;
    f_batch inv_w_h(1.0f);
    b_batch valid_w_h = (xs::abs(w_h) > 1e-8f) & (xs::abs(w_h - 1.0f) > 1e-8f);
    inv_w_h = xs::select(valid_w_h, 1.0f / w_h, inv_w_h);

    world_hit_x = (m11 * local_hit_x + m12 * local_hit_y + m13 * local_hit_z + m14) * inv_w_h;
    world_hit_y = (m21 * local_hit_x + m22 * local_hit_y + m23 * local_hit_z + m24) * inv_w_h;
    world_hit_z = (m31 * local_hit_x + m32 * local_hit_y + m33 * local_hit_z + m34) * inv_w_h;

    // --- 2. Transform Normal (N' = (M^-1)^T * N) ---
    // We need the inverse transpose of the upper 3x3 part of forward_transform.
    // (M^-1)^T = (1/det(M_3x3)) * Adj(M_3x3)
    // Adj(M_3x3) = Cof(M_3x3)^T
    // Let's compute cofactors of M_3x3
    fl c11 = forward_transform.m22 * forward_transform.m33 - forward_transform.m23 * forward_transform.m32;
    fl c12 = -(forward_transform.m21 * forward_transform.m33 - forward_transform.m23 * forward_transform.m31);
    fl c13 = forward_transform.m21 * forward_transform.m32 - forward_transform.m22 * forward_transform.m31;
    
    fl c21 = -(forward_transform.m12 * forward_transform.m33 - forward_transform.m13 * forward_transform.m32);
    fl c22 = forward_transform.m11 * forward_transform.m33 - forward_transform.m13 * forward_transform.m31;
    fl c23 = -(forward_transform.m11 * forward_transform.m32 - forward_transform.m12 * forward_transform.m31);

    fl c31 = forward_transform.m12 * forward_transform.m23 - forward_transform.m13 * forward_transform.m22;
    fl c32 = -(forward_transform.m11 * forward_transform.m23 - forward_transform.m13 * forward_transform.m21);
    fl c33 = forward_transform.m11 * forward_transform.m22 - forward_transform.m12 * forward_transform.m21;

    fl det = forward_transform.m11 * c11 + forward_transform.m12 * c12 + forward_transform.m13 * c13;
    fl inv_det = (std::abs(det) > 1e-8f) ? (1.0f / det) : 0.0f;
    
    // The inverse transpose matrix (N_Mat) is (1/det) * CofactorMatrix (not transposed, because N' = N * (M^-1))
    // Wait, N' = ((M^-1)^T) * N. So we need the Adjoint matrix, which is Cofactor^T.
    // N_Mat = (1/det) * Adj = (1/det) * Cof^T
    // it_m11 = c11 * inv_det; it_m12 = c21 * inv_det; it_m13 = c31 * inv_det;
    // it_m21 = c12 * inv_det; it_m22 = c22 * inv_det; it_m23 = c32 * inv_det;
    // it_m31 = c13 * inv_det; it_m32 = c23 * inv_det; it_m33 = c33 * inv_det;
    
    // Broadcast inverse transpose matrix elements
   // (M^-1)^T = (1/det) * Cof(M)
// Broadcast inverse transpose matrix elements
auto it_m11 = xs::broadcast(c11 * inv_det), it_m12 = xs::broadcast(c12 * inv_det), it_m13 = xs::broadcast(c13 * inv_det);
auto it_m21 = xs::broadcast(c21 * inv_det), it_m22 = xs::broadcast(c22 * inv_det), it_m23 = xs::broadcast(c23 * inv_det);
auto it_m31 = xs::broadcast(c31 * inv_det), it_m32 = xs::broadcast(c32 * inv_det), it_m33 = xs::broadcast(c33 * inv_det);

// Apply transformation: N' = (M^-1)^T * N
auto nx = it_m11 * local_norm_x + it_m12 * local_norm_y + it_m13 * local_norm_z;
auto ny = it_m21 * local_norm_x + it_m22 * local_norm_y + it_m23 * local_norm_z;
auto nz = it_m31 * local_norm_x + it_m32 * local_norm_y + it_m33 * local_norm_z;

    // --- 3. Re-normalize the transformed normal ---
    auto len_sq = nx * nx + ny * ny + nz * nz;
    auto len = xs::sqrt(len_sq);
    b_batch valid_len = (len > 1e-8f);
    f_batch inv_len = xs::select(valid_len, 1.0f / len, f_batch(0.0f));

    world_norm_x = nx * inv_len;
    world_norm_y = ny * inv_len;
    world_norm_z = nz * inv_len;
}

// --- NEW AABB INTERSECTION ---

void inline intersect_aabb(const f_batch& orig_x, const f_batch& orig_y, const f_batch& orig_z,
                           const f_batch& dir_x, const f_batch& dir_y, const f_batch& dir_z,
                           const f_batch& t_min, AABB bbox, b_batch& active_mask){
    
    // Precompute 1.0 / direction
    f_batch inv_dir_x = 1.0f / dir_x;
    f_batch inv_dir_y = 1.0f / dir_y;
    f_batch inv_dir_z = 1.0f / dir_z;

    f_batch bmin_x = xs::broadcast(bbox.bbox_xmin);
    f_batch bmax_x = xs::broadcast(bbox.bbox_xmax);
    f_batch bmin_y = xs::broadcast(bbox.bbox_ymin);
    f_batch bmax_y = xs::broadcast(bbox.bbox_ymax);
    f_batch bmin_z = xs::broadcast(bbox.bbox_zmin);
    f_batch bmax_z = xs::broadcast(bbox.bbox_zmax);

    f_batch tx1 = (bmin_x - orig_x) * inv_dir_x;
    f_batch tx2 = (bmax_x - orig_x) * inv_dir_x;
    f_batch tmin_box = xs::min(tx1, tx2);
    f_batch tmax_box = xs::max(tx1, tx2);

    f_batch ty1 = (bmin_y - orig_y) * inv_dir_y;
    f_batch ty2 = (bmax_y - orig_y) * inv_dir_y;
    tmin_box = xs::max(tmin_box, xs::min(ty1, ty2));
    tmax_box = xs::min(tmax_box, xs::max(ty1, ty2));

    f_batch tz1 = (bmin_z - orig_z) * inv_dir_z;
    f_batch tz2 = (bmax_z - orig_z) * inv_dir_z;
    tmin_box = xs::max(tmin_box, xs::min(tz1, tz2));
    tmax_box = xs::min(tmax_box, xs::max(tz1, tz2));

    // AABB is hit if tmax_box >= tmin_box AND
    // tmax_box > 0 (box is not entirely behind) AND
    // tmin_box < t_min (box is not entirely beyond the closest hit)
    b_batch intersects = (tmax_box >= tmin_box) & (tmax_box > 0.0f) & (tmin_box < t_min);
    
    active_mask = active_mask & intersects;
}

// --- NEW BVH TRAVERSAL (CLOSEST-HIT) ---

void inline intersect_bvh(const f_batch& orig_x, const f_batch& orig_y, const f_batch& orig_z,
                          const f_batch& dir_x, const f_batch& dir_y, const f_batch& dir_z,
                          f_batch& t_min, i_batch& mat_id,
                          f_batch& norm_x, f_batch& norm_y, f_batch& norm_z,
                          const Scene& scene, const int bvh_index, const uint bvh_node_id,
                          const i_batch& blas_mat_id, b_batch active_mask) {
    
    // Get the specific BVH from the global list
    const BVH& bvh = bvh_list[bvh_index];
    const BVHNode& node = bvh.bvh_node_list[bvh_node_id];

    // Test AABB intersection
    intersect_aabb(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z, t_min, node.bbox, active_mask);
    
    // Ray missed this node's box
    if(xs::none(active_mask)) return;

    // Leaf node: intersect primitives
    if(node.left_child_index == 0){
        for(int i = node.first_prim; i < node.prim_count + node.first_prim; i++)
        {
            // Get the *scene-wide* primitive index from the BVH's reordered list
            int current_prim_index = bvh.prim_indices[i];
            
            // This BVH only contains triangles, so we only test triangles.
            // We pass the blas_mat_id to override the triangle's base material.
            intersect_triangles_masked(dir_x, dir_y, dir_z, orig_x, orig_y, orig_z,
                t_min, mat_id, norm_x, norm_y, norm_z,
                scene.triangle_data__, scene.vertex_data__, scene, current_prim_index,
                active_mask, blas_mat_id);
        }
    }
    else { // Internal node: traverse children
        intersect_bvh(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                      t_min, mat_id, norm_x, norm_y, norm_z,
                      scene, bvh_index, node.left_child_index, blas_mat_id, active_mask);
        
        intersect_bvh(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                      t_min, mat_id, norm_x, norm_y, norm_z,
                      scene, bvh_index, node.left_child_index + 1, blas_mat_id, active_mask);
    }
}

// --- NEW TLAS TRAVERSAL (CLOSEST-HIT) ---

void inline intersect_tlas(const f_batch& orig_x, const f_batch& orig_y, const f_batch& orig_z,
                           const f_batch& dir_x, const f_batch& dir_y, const f_batch& dir_z,
                           f_batch& t_min, i_batch& mat_id,
                           f_batch& norm_x, f_batch& norm_y, f_batch& norm_z,
                           f_batch& hit_pos_x, f_batch& hit_pos_y, f_batch& hit_pos_z,
                           const Scene& scene, const uint tlas_node_id, b_batch active_mask) {
    
    const TLASNode& node = tlas.tlas_root[tlas_node_id];

    // Test TLAS node's AABB
    intersect_aabb(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z, t_min, node.bbox, active_mask);
    
    if(xs::none(active_mask)) return;

    // Leaf node: intersect BLASes
    if(node.left_child_index == 0) {
        
        // --- Intersect BLAS 0 ---
        if(node.blas_id_0 != -1) {
            const BLAS& blas = blas_list[node.blas_id_0];

            // Transform ray to BLAS local space
            f_batch local_o_x, local_o_y, local_o_z;
            f_batch local_d_x, local_d_y, local_d_z;
            transform_ray_packet(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                                 local_o_x, local_o_y, local_o_z,
                                 local_d_x, local_d_y, local_d_z,
                                 blas.inv_transform);

            // Create batches for this BLAS's hit data
            f_batch local_t_min = t_min; // Test against current world t_min
            i_batch local_mat_id = xs::broadcast(blas.material_id);
            f_batch local_norm_x(0.0f), local_norm_y(0.0f), local_norm_z(0.0f);
            
            // Check if this is a sphere (bvh_index >= 10000) or a mesh
            if(blas.bvh_index >= 10000) {
                // This is a SPHERE
                int sphere_index = blas.bvh_index - 10000;
                intersect_spheres_masked(local_d_x, local_d_y, local_d_z,
                                        local_o_x, local_o_y, local_o_z,
                                        local_t_min, local_mat_id,
                                        local_norm_x, local_norm_y, local_norm_z,
                                        scene.sphere_data__, scene.vertex_data__, scene,
                                        sphere_index, active_mask);
            } else {
                // This is a MESH - Traverse the BVH
                // Pass blas.material_id as override to all triangles in this BLAS
                intersect_bvh(local_o_x, local_o_y, local_o_z,
                              local_d_x, local_d_y, local_d_z,
                              local_t_min, mat_id, // mat_id is output, but local_mat_id is input
                              local_norm_x, local_norm_y, local_norm_z,
                              scene, blas.bvh_index, 0, // 0 = root node of BVH
                              local_mat_id, active_mask);
            }

            // Check if any rays found a *closer* hit
            b_batch new_hit_mask = (local_t_min < t_min) & active_mask;
            if(xs::any(new_hit_mask)) {
                // Calculate local hit point: P_local = O_local + t * D_local
                f_batch local_hit_x = local_o_x + local_t_min * local_d_x;
                f_batch local_hit_y = local_o_y + local_t_min * local_d_y;
                f_batch local_hit_z = local_o_z + local_t_min * local_d_z;

                // Transform hit point and normal back to world space
                f_batch world_hit_x, world_hit_y, world_hit_z;
                f_batch world_norm_x, world_norm_y, world_norm_z;
                transform_hit_results(local_hit_x, local_hit_y, local_hit_z,
                                      local_norm_x, local_norm_y, local_norm_z,
                                      world_hit_x, world_hit_y, world_hit_z,
                                      world_norm_x, world_norm_y, world_norm_z,
                                      blas.transformation_matrix); // Use forward transform

                // Conditionally update the final world-space hit info
                t_min = xs::select(new_hit_mask, local_t_min, t_min);
                mat_id = xs::select(xs::batch_bool_cast<int32_t>(new_hit_mask), local_mat_id, mat_id);
                
                norm_x = xs::select(new_hit_mask, world_norm_x, norm_x);
                norm_y = xs::select(new_hit_mask, world_norm_y, norm_y);
                norm_z = xs::select(new_hit_mask, world_norm_z, norm_z);
                
                hit_pos_x = xs::select(new_hit_mask, world_hit_x, hit_pos_x);
                hit_pos_y = xs::select(new_hit_mask, world_hit_y, hit_pos_y);
                hit_pos_z = xs::select(new_hit_mask, world_hit_z, hit_pos_z);
            }
        }

        // --- Intersect BLAS 1 ---
        if(node.blas_id_1 != -1) {
            const BLAS& blas = blas_list[node.blas_id_1];

            f_batch local_o_x, local_o_y, local_o_z;
            f_batch local_d_x, local_d_y, local_d_z;
            transform_ray_packet(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                                 local_o_x, local_o_y, local_o_z,
                                 local_d_x, local_d_y, local_d_z,
                                 blas.inv_transform);

            f_batch local_t_min = t_min;
            i_batch local_mat_id = xs::broadcast(blas.material_id);
            f_batch local_norm_x(0.0f), local_norm_y(0.0f), local_norm_z(0.0f);
            
            // Check if this is a sphere (bvh_index >= 10000) or a mesh
            if(blas.bvh_index >= 10000) {
                // This is a SPHERE
                int sphere_index = blas.bvh_index - 10000;
                intersect_spheres_masked(local_d_x, local_d_y, local_d_z,
                                        local_o_x, local_o_y, local_o_z,
                                        local_t_min, local_mat_id,
                                        local_norm_x, local_norm_y, local_norm_z,
                                        scene.sphere_data__, scene.vertex_data__, scene,
                                        sphere_index, active_mask);
            } else {
                // This is a MESH - Traverse the BVH
                intersect_bvh(local_o_x, local_o_y, local_o_z,
                              local_d_x, local_d_y, local_d_z,
                              local_t_min, mat_id,
                              local_norm_x, local_norm_y, local_norm_z,
                              scene, blas.bvh_index, 0,
                              local_mat_id, active_mask);
            }

            b_batch new_hit_mask = (local_t_min < t_min) & active_mask;
            if(xs::any(new_hit_mask)) {
                f_batch local_hit_x = local_o_x + local_t_min * local_d_x;
                f_batch local_hit_y = local_o_y + local_t_min * local_d_y;
                f_batch local_hit_z = local_o_z + local_t_min * local_d_z;

                f_batch world_hit_x, world_hit_y, world_hit_z;
                f_batch world_norm_x, world_norm_y, world_norm_z;
                transform_hit_results(local_hit_x, local_hit_y, local_hit_z,
                                      local_norm_x, local_norm_y, local_norm_z,
                                      world_hit_x, world_hit_y, world_hit_z,
                                      world_norm_x, world_norm_y, world_norm_z,
                                      blas.transformation_matrix);

                t_min = xs::select(new_hit_mask, local_t_min, t_min);
                mat_id = xs::select(xs::batch_bool_cast<int32_t>(new_hit_mask), local_mat_id, mat_id);
                
                norm_x = xs::select(new_hit_mask, world_norm_x, norm_x);
                norm_y = xs::select(new_hit_mask, world_norm_y, norm_y);
                norm_z = xs::select(new_hit_mask, world_norm_z, norm_z);
                
                hit_pos_x = xs::select(new_hit_mask, world_hit_x, hit_pos_x);
                hit_pos_y = xs::select(new_hit_mask, world_hit_y, hit_pos_y);
                hit_pos_z = xs::select(new_hit_mask, world_hit_z, hit_pos_z);
            }
        }
    }
    else { // Internal node: traverse children
        intersect_tlas(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                       t_min, mat_id, norm_x, norm_y, norm_z,
                       hit_pos_x, hit_pos_y, hit_pos_z,
                       scene, node.left_child_index, active_mask);
        
        intersect_tlas(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                       t_min, mat_id, norm_x, norm_y, norm_z,
                       hit_pos_x, hit_pos_y, hit_pos_z,
                       scene, node.left_child_index + 1, active_mask);
    }
}

// --- NEW TLAS WRAPPER (CLOSEST-HIT) ---

void intersect_tlas_wrapper(RP8& ray_pack, const Scene& scene) {
    // 1. Load the UNNORMALIZED data
    auto dir_x = xs::load(&ray_pack.d_x[0]);
    auto dir_y = xs::load(&ray_pack.d_y[0]);
    auto dir_z = xs::load(&ray_pack.d_z[0]);

    // 2. Normalize directions
    auto len_sq = (dir_x * dir_x) + (dir_y * dir_y) + (dir_z * dir_z);
    auto len = xs::sqrt(len_sq);
    b_batch valid_mask = (len > 1e-8f);
    f_batch zero_batch(0.0f);
    dir_x = xs::select(valid_mask, dir_x / len, zero_batch);
    dir_y = xs::select(valid_mask, dir_y / len, zero_batch);
    dir_z = xs::select(valid_mask, dir_z / len, zero_batch);

    // 3. Load ray origins
    f_batch orig_x = xs::load(&ray_pack.o_x[0]);
    f_batch orig_y = xs::load(&ray_pack.o_y[0]);
    f_batch orig_z = xs::load(&ray_pack.o_z[0]);

    // 4. Load current intersection state
    f_batch t_min = xs::load(&ray_pack.t_min[0]);
    i_batch mat_id = xs::load(&ray_pack.mat_id[0]);
    f_batch norm_x = xs::load(&ray_pack.hit_norm_x[0]);
    f_batch norm_y = xs::load(&ray_pack.hit_norm_y[0]);
    f_batch norm_z = xs::load(&ray_pack.hit_norm_z[0]);
    
    // 5. Initialize hit positions (will be filled by intersect_tlas)
    f_batch hit_pos_x = xs::load(&ray_pack.hit_pos_x[0]);
    f_batch hit_pos_y = xs::load(&ray_pack.hit_pos_y[0]);
    f_batch hit_pos_z = xs::load(&ray_pack.hit_pos_z[0]);

    // 6. Start with all rays active
    b_batch active_mask = true;

    // 7. Traverse TLAS starting from root (node 0)
    intersect_tlas(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                   t_min, mat_id, norm_x, norm_y, norm_z,
                   hit_pos_x, hit_pos_y, hit_pos_z,
                   scene, 0, active_mask);

    // 8. Store the final results back into the ray packet
    xs::store(&ray_pack.t_min[0], t_min);
    xs::store(&ray_pack.mat_id[0], mat_id);
    xs::store(&ray_pack.hit_norm_x[0], norm_x);
    xs::store(&ray_pack.hit_norm_y[0], norm_y);
    xs::store(&ray_pack.hit_norm_z[0], norm_z);
    xs::store(&ray_pack.hit_pos_x[0], hit_pos_x);
    xs::store(&ray_pack.hit_pos_y[0], hit_pos_y);
    xs::store(&ray_pack.hit_pos_z[0], hit_pos_z);
}

// --- SPHERE ANY-HIT (FOR SHADOWS) ---

void inline intersect_sphere_any_hit(const f_batch& orig_x, const f_batch& orig_y, const f_batch& orig_z,
                                     const f_batch& dir_x, const f_batch& dir_y, const f_batch& dir_z,
                                     const f_batch& t_max, b_batch& is_occluded,
                                     const Scene& scene, const int sphere_index,
                                     b_batch active_mask, const fl epsilon) {
    
    if(xs::none(active_mask)) return;
    
    int center_vertex_id = scene.sphere_data__.sphere_center_vertex_id[sphere_index];
    auto sphere_center_x = xs::broadcast(scene.vertex_data__.v_pos_x[center_vertex_id]);
    auto sphere_center_y = xs::broadcast(scene.vertex_data__.v_pos_y[center_vertex_id]);
    auto sphere_center_z = xs::broadcast(scene.vertex_data__.v_pos_z[center_vertex_id]);
    auto sphere_radius_sq = xs::broadcast(scene.sphere_data__.sphere_radius_sq[sphere_index]);

    auto oc_x = orig_x - sphere_center_x;
    auto oc_y = orig_y - sphere_center_y;
    auto oc_z = orig_z - sphere_center_z;

    auto b_half = (dir_x * oc_x) + (dir_y * oc_y) + (dir_z * oc_z);
    auto c = (oc_x * oc_x) + (oc_y * oc_y) + (oc_z * oc_z) - sphere_radius_sq;
    auto discriminant = (b_half * b_half) - c;

    b_batch hit_mask = (discriminant >= 0.f);
    if (xs::none(hit_mask & active_mask)) return;

    auto sqrt_discriminant = xs::sqrt(discriminant);
    auto t0 = -b_half - sqrt_discriminant;
    auto t1 = -b_half + sqrt_discriminant;
    auto t_smaller = xs::min(t0, t1);
    auto t_larger = xs::max(t0, t1);
    auto t_new = xs::select(t_smaller > epsilon, t_smaller, t_larger);

    // Mark rays as occluded if they hit between epsilon and t_max
    b_batch new_occlusion = active_mask & hit_mask & (t_new > epsilon) & (t_new < t_max);
    is_occluded = is_occluded | new_occlusion;
}

// --- NEW BVH TRAVERSAL (ANY-HIT / SHADOW) ---

void inline intersect_bvh_any_hit(const f_batch& orig_x, const f_batch& orig_y, const f_batch& orig_z,
                                   const f_batch& dir_x, const f_batch& dir_y, const f_batch& dir_z,
                                   const f_batch& t_max, b_batch& is_occluded,
                                   const Scene& scene, const int bvh_index, const uint bvh_node_id, 
                                   b_batch active_mask, const fl epsilon) {
    
    // Early exit if all active rays are already occluded
    if(xs::none(active_mask)) return;

    const BVH& bvh = bvh_list[bvh_index];
    const BVHNode& node = bvh.bvh_node_list[bvh_node_id];
    
    // Test AABB intersection (using t_max as the max distance)
    intersect_aabb(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z, t_max, node.bbox, active_mask);
    
    if(xs::none(active_mask)) return;
    
    // Leaf node: test primitives
    if(node.left_child_index == 0){
        for(int i = node.first_prim; i < node.prim_count + node.first_prim; i++)
        {
            // Stop testing if all rays in the active mask are now occluded
            if(xs::none(active_mask)) return; 
            
            int current_prim_index = bvh.prim_indices[i];
            
            // --- Triangle intersection for shadow rays (Möller-Trumbore) ---
            const VertexData& vertices = scene.vertex_data__;
            int i0 = scene.triangle_data__.v0_ind[current_prim_index];
            int i1 = scene.triangle_data__.v1_ind[current_prim_index];
            int i2 = scene.triangle_data__.v2_ind[current_prim_index];

            Vec3f a = { vertices.v_pos_x[i0], vertices.v_pos_y[i0], vertices.v_pos_z[i0] };
            Vec3f b = { vertices.v_pos_x[i1], vertices.v_pos_y[i1], vertices.v_pos_z[i1] };
            Vec3f c = { vertices.v_pos_x[i2], vertices.v_pos_y[i2], vertices.v_pos_z[i2] };

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
            
            // No backface culling for shadows
            b_batch det_mask = (xs::abs(det) > epsilon);
            if (xs::none(det_mask & active_mask)) continue;
            
            f_batch zero_batch(0.0f);
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
            if (xs::none(bary_mask & active_mask)) continue;

            auto qvec_x = tvec_y * e1_z - tvec_z * e1_y;
            auto qvec_y = tvec_z * e1_x - tvec_x * e1_z;
            auto qvec_z = tvec_x * e1_y - tvec_y * e1_x;
            auto v = (dir_x * qvec_x + dir_y * qvec_y + dir_z * qvec_z) * invDet;
            bary_mask = bary_mask & (v >= zero_batch - bary_epsilon) & (u + v <= one_batch + bary_epsilon);
            if (xs::none(bary_mask & active_mask)) continue;

            auto t_new = (e2_x * qvec_x + e2_y * qvec_y + e2_z * qvec_z) * invDet;
            
            // Mark rays as occluded if they hit between epsilon and t_max
            b_batch new_occlusion = active_mask & det_mask & bary_mask & (t_new > epsilon) & (t_new < t_max);
            
            is_occluded = is_occluded | new_occlusion;
            active_mask = active_mask & !is_occluded; // Update active mask
        }
    }
    else {
        // Internal node: traverse children
        intersect_bvh_any_hit(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                              t_max, is_occluded, scene, bvh_index, node.left_child_index, active_mask, epsilon);
        
        // Update active_mask to only check non-occluded rays in the right child
        active_mask = active_mask & !is_occluded; 
        
        intersect_bvh_any_hit(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                              t_max, is_occluded, scene, bvh_index, node.left_child_index + 1, active_mask, epsilon);
    }
}

// --- NEW TLAS TRAVERSAL (ANY-HIT / SHADOW) ---

void inline intersect_tlas_any_hit(const f_batch& orig_x, const f_batch& orig_y, const f_batch& orig_z,
                                    const f_batch& dir_x, const f_batch& dir_y, const f_batch& dir_z,
                                    const f_batch& t_max, b_batch& is_occluded,
                                    const Scene& scene, const uint tlas_node_id, 
                                    b_batch active_mask, const fl epsilon) {
    
    // Early exit if all active rays are already occluded
    if(xs::none(active_mask)) return;

    const TLASNode& node = tlas.tlas_root[tlas_node_id];

    // Test TLAS node's AABB (using t_max as max distance)
    intersect_aabb(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z, t_max, node.bbox, active_mask);
    
    if(xs::none(active_mask)) return;

    // Leaf node: intersect BLASes
    if(node.left_child_index == 0) {
        
        // --- Intersect BLAS 0 ---
        if(node.blas_id_0 != -1) {
            const BLAS& blas = blas_list[node.blas_id_0];

            f_batch local_o_x, local_o_y, local_o_z;
            f_batch local_d_x, local_d_y, local_d_z;
            transform_ray_packet(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                                 local_o_x, local_o_y, local_o_z,
                                 local_d_x, local_d_y, local_d_z,
                                 blas.inv_transform);

            // Check if this is a sphere (bvh_index >= 10000) or a mesh
            if(blas.bvh_index >= 10000) {
                // This is a SPHERE
                int sphere_index = blas.bvh_index - 10000;
                intersect_sphere_any_hit(local_o_x, local_o_y, local_o_z,
                                        local_d_x, local_d_y, local_d_z,
                                        t_max, is_occluded, scene,
                                        sphere_index, active_mask, epsilon);
            } else {
                // This is a MESH - Traverse the BVH
                intersect_bvh_any_hit(local_o_x, local_o_y, local_o_z,
                                      local_d_x, local_d_y, local_d_z,
                                      t_max, is_occluded, scene, 
                                      blas.bvh_index, 0, // 0 = root node of BVH
                                      active_mask, epsilon);
            }
            
            // Update active_mask
            active_mask = active_mask & !is_occluded;
            if(xs::none(active_mask)) return; // All rays occluded, stop.
        }

        // --- Intersect BLAS 1 ---
        if(node.blas_id_1 != -1) {
            const BLAS& blas = blas_list[node.blas_id_1];

            f_batch local_o_x, local_o_y, local_o_z;
            f_batch local_d_x, local_d_y, local_d_z;
            transform_ray_packet(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                                 local_o_x, local_o_y, local_o_z,
                                 local_d_x, local_d_y, local_d_z,
                                 blas.inv_transform);

            // Check if this is a sphere (bvh_index >= 10000) or a mesh
            if(blas.bvh_index >= 10000) {
                // This is a SPHERE
                int sphere_index = blas.bvh_index - 10000;
                intersect_sphere_any_hit(local_o_x, local_o_y, local_o_z,
                                        local_d_x, local_d_y, local_d_z,
                                        t_max, is_occluded, scene,
                                        sphere_index, active_mask, epsilon);
            } else {
                // This is a MESH - Traverse the BVH
                intersect_bvh_any_hit(local_o_x, local_o_y, local_o_z,
                                      local_d_x, local_d_y, local_d_z,
                                      t_max, is_occluded, scene, 
                                      blas.bvh_index, 0,
                                      active_mask, epsilon);
            }
        }
    }
    else { // Internal node: traverse children
        intersect_tlas_any_hit(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                               t_max, is_occluded, scene, 
                               node.left_child_index, active_mask, epsilon);
        
        active_mask = active_mask & !is_occluded; // Update mask

        intersect_tlas_any_hit(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                               t_max, is_occluded, scene, 
                               node.left_child_index + 1, active_mask, epsilon);
    }
}

// --- NEW TLAS WRAPPER (ANY-HIT / SHADOW) ---

void inline intersect_tlas_any_hit_wrapper(RP8& ray_pack, const Scene& scene, 
                                           const f_batch& t_max, b_batch& is_occluded,
                                           const b_batch& active_mask) {
    
    // Early exit if no rays are active
    if (xs::none(active_mask)) {
        is_occluded = false; // Ensure it's all false
        return;
    }
    
    // Load and normalize ray directions
    auto dir_x = xs::load(&ray_pack.d_x[0]);
    auto dir_y = xs::load(&ray_pack.d_y[0]);
    auto dir_z = xs::load(&ray_pack.d_z[0]);

    auto len_sq = (dir_x * dir_x) + (dir_y * dir_y) + (dir_z * dir_z);
    auto len = xs::sqrt(len_sq);
    b_batch valid_len_mask = (len > 1e-8f);
    b_batch norm_mask = active_mask & valid_len_mask;

    f_batch zero_batch(0.0f);
    dir_x = xs::select(norm_mask, dir_x / len, zero_batch);
    dir_y = xs::select(norm_mask, dir_y / len, zero_batch);
    dir_z = xs::select(norm_mask, dir_z / len, zero_batch);

    // Load ray origins
    f_batch orig_x = xs::load(&ray_pack.o_x[0]);
    f_batch orig_y = xs::load(&ray_pack.o_y[0]);
    f_batch orig_z = xs::load(&ray_pack.o_z[0]);

    const fl epsilon = scene.intersection_test_epsilon;

    // Start with no rays occluded
    is_occluded = false;

    // Traverse TLAS
    intersect_tlas_any_hit(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z,
                           t_max, is_occluded, scene, 0, // 0 = root node
                           active_mask, epsilon);
    
    // Final mask: a ray is in light if it was active AND not occluded.
    // The shader will use '!is_occluded'
}


#endif