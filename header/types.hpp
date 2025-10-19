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


struct Vertex{
    Vec3f position;
    //will be used for smooth shading
    Vec3f normal;
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

// Represents a point light source in the scene.
struct PointLight {
    Vec3f position;
    Vec3f intensity;
    int id;
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

struct HitRecord{
    Vec3f hit_point;
    Vec3f normal;
    fl t;
    int material_id;
    };


// Represents a sphere object.
struct Sphere {
    fl radius;
    int id;
    int material_id;
    int center_vertex_id;
};

// Represents a single triangle.
struct Triangle {
    Vec3i indices; // Vertex indices
    Vec3f normal;
    int id;
    int material_id;
};
struct Face{
    Vec3i indices; // Vertex indices
    //will be used for flat shading.
    Vec3f normal;
    
};

// Represents a mesh object composed of multiple triangles (faces).
struct Mesh {
    std::vector<Face> faces; // Each Vec3i contains the vertex indices for one triangle face
    int id;
    int material_id;
};

// Represents an infinite plane object.
struct Plane {
    Vec3f normal;
    int id;
    int material_id;
    int point_vertex_id; // ID of a vertex that lies on the plane
};

// The main struct that holds all scene data, parsed from the JSON file.
struct Scene {
    // Scene elements
    std::vector<Camera> cameras;
    std::vector<PointLight> point_lights;
    std::vector<Material> materials;
    std::vector<Vertex> vertex_data;
    std::vector<Sphere> spheres;
    std::vector<Triangle> triangles;
    std::vector<Mesh> meshes;
    std::vector<Plane> planes;


    VertexData vertex_data__;
    PointLightData point_light_data__;
    SphereData sphere_data__;
    TriangleData triangle_data__;

    // General scene settings
    Vec3f background_color;
    Vec3f ambient_light;
    fl shadow_ray_epsilon;
    fl intersection_test_epsilon;
    int max_recursion_depth;
};

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
    std::vector<int>  sphere_id;
    std::vector<int> sphere_mat_id;
    std::vector<int> sphere_center_vertex_id;
    std::vector<float> sphere_radius;
};

struct TriangleData{
    std::vector<int> vo_ind;
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




