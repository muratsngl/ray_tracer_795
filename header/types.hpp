#pragma once

#include <vector>
#include <string>


typedef float fl;

// 3D vector of floats.
struct Vec3f {
    fl x, y, z;
};

// 3D vector of integers.
struct Vec3i {
    int v0, v1, v2;
};

// 2D vector of integers.
struct Vec2i {
    int x, y;
};

// 4D vector of floats.
struct Vec4f {
    fl x, y, z, w;
};

struct Mat3f{
    fl m11=0,m12=0,m13=0,
       m21=0,m22=0,m23=0,
       m31=0,m32=0,m33=0;
};

// Holds all properties related to a single material.
struct Material {
    std::string type;
    Vec3f ambient_reflectance;
    Vec3f diffuse_reflectance;
    Vec3f specular_reflectance;
    Vec3f mirror_reflectance;
    Vec3f absorption_coefficient;
    fl phong_exponent;
    fl refraction_index;
    fl absorption_index;
    int id;
};

struct ColorBlock{
    unsigned char rgb[24]={};
};

struct ColorBlockFl{
    fl rgb[24]={};
};

// Defines the properties of a single camera.
struct Camera {
    std::string image_name;
    Vec4f near_plane; // Left, Right, Bottom, Top
    Vec3f position;
    Vec3f gaze;
    Vec3f up;
    Vec2i image_resolution;
    fl near_distance;
    int id;
};

//Ray Structure
struct Ray{
    Vec3f origin;
    Vec3f direction;
};

//RAY PACKETS FOR PARALLELIZATION USING SIMD
struct alignas(32) RP8{
    fl o_x[8] = {0.f};
    fl o_y[8] = {0.f};
    fl o_z[8] = {0.f};
    fl d_x[8] = {0.f};
    fl d_y[8] = {0.f};
    fl d_z[8] = {0.f};

    fl t_min[8] = {__builtin_inff(), __builtin_inff(), __builtin_inff(), __builtin_inff(), __builtin_inff(), __builtin_inff(), __builtin_inff(), __builtin_inff()};
    fl hit_norm_x[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    fl hit_norm_y[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    fl hit_norm_z[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    fl hit_pos_x[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    fl hit_pos_y[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    fl hit_pos_z[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    int mat_id[8] = {-1, -1, -1, -1, -1, -1, -1, -1};

};

struct HitRecord{
    Vec3f hit_point;
    Vec3f normal;
    fl t;
    int material_id;
};

// Represents an infinite plane object (still using AoS - not many planes in typical scenes).
struct Plane {
    Vec3f normal;
    int id;
    int material_id;
    int point_vertex_id; // ID of a vertex that lies on the plane
};

// SoA (Structure of Arrays) data structures
struct VertexData{
    std::vector<fl> v_pos_x;
    std::vector<fl> v_pos_y;
    std::vector<fl> v_pos_z;

    std::vector<fl> v_nor_x;
    std::vector<fl> v_nor_y;
    std::vector<fl> v_nor_z;
};

struct PointLightData{
    std::vector<fl> pl_pos_x;
    std::vector<fl> pl_pos_y;
    std::vector<fl> pl_pos_z;
    
    std::vector<fl> pl_intensity_r;
    std::vector<fl> pl_intensity_g;
    std::vector<fl> pl_intensity_b;

    std::vector<int> pl_id;
};

struct SphereData{
    std::vector<int> sphere_id;
    std::vector<int> sphere_mat_id;
    std::vector<int> sphere_center_vertex_id;
    std::vector<fl> sphere_radius_sq;
};

struct TriangleData{
    std::vector<int> v0_ind;
    std::vector<int> v1_ind;
    std::vector<int> v2_ind;

    std::vector<fl> tri_norm_x;
    std::vector<fl> tri_norm_y;
    std::vector<fl> tri_norm_z;

    std::vector<int> triangle_id;
    std::vector<int> triangle_material_id;
};

struct PlaneData{
    std::vector<fl> plane_norm_x;
    std::vector<fl> plane_norm_y;
    std::vector<fl> plane_norm_z;

    std::vector<int> plane_point_vertex_id;
    
    std::vector<int> plane_id;

    std::vector<int> plane_material_id;
};

// The main struct that holds all scene data, parsed from the JSON file.
struct Scene {
    // SoA (Structure of Arrays) for data-oriented design
    VertexData vertex_data__;
    PointLightData point_light_data__;
    SphereData sphere_data__;
    TriangleData triangle_data__;
    PlaneData plane_data__;

    
    std::vector<Camera> cameras;
    std::vector<Material> materials;
    
    // General scene settings
    Vec3f background_color;
    Vec3f ambient_light;
    fl shadow_ray_epsilon;
    fl intersection_test_epsilon;
    int max_recursion_depth;
};




