#pragma once

#include <vector>
#include <string>

#include <map>



//global state
inline bool G_SMOOTH_SHADING_ENABLED = false;
inline int solo_triangle_counter = 0;


typedef float fl;

typedef enum PrimType{
    TRIANGLE,
    SPHERE
}pt;

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
struct Mat4f{
    fl m11=0,m12=0,m13=0,m14=0,
       m21=0,m22=0,m23=0,m24=0,
       m31=0,m32=0,m33=0,m34=0,
       m41=0,m42=0,m43=0,m44=0;
    };

struct AABB{
    fl bbox_xmin;
    fl bbox_ymin;
    fl bbox_zmin;

    fl bbox_xmax;
    fl bbox_ymax;
    fl bbox_zmax;
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
    fl roughness;
    int id;
};

struct ColorBlock{
    unsigned char rgb[24]={};
};

struct ColorBlockFl{
    fl rgb[24]={};
};



struct SoARayQueue{
    std::vector<fl> o_x;
    std::vector<fl> o_y;
    std::vector<fl> o_z;
    std::vector<fl> d_x;
    std::vector<fl> d_y;
    std::vector<fl> d_z;
    std::vector<fl> tp_r;
    std::vector<fl> tp_g;
    std::vector<fl> tp_b;
    std::vector<fl> area_light_u;
    std::vector<fl> area_light_v;
    std::vector<int> depth;
    std::vector<int> pixel_index; // Index to track which pixel (0-7) this ray belongs to
    std::vector<fl> time;
    
    void push(fl ox, fl oy, fl oz, fl dx, fl dy, fl dz, fl tr, fl tg, fl tb, int d, int idx,fl t) {
        o_x.push_back(ox); o_y.push_back(oy); o_z.push_back(oz);
        d_x.push_back(dx); d_y.push_back(dy); d_z.push_back(dz);
        tp_r.push_back(tr); tp_g.push_back(tg); tp_b.push_back(tb);
        depth.push_back(d);
        pixel_index.push_back(idx);time.push_back(t);
    }

    bool pop(fl& ox, fl& oy, fl& oz, fl& dx, fl& dy, fl& dz, fl& tr, fl& tg, fl& tb, int& d, int& idx) {
        if (is_empty()) return false;
        ox = o_x.back(); o_x.pop_back();
        oy = o_y.back(); o_y.pop_back();
        oz = o_z.back(); o_z.pop_back();
        dx = d_x.back(); d_x.pop_back();
        dy = d_y.back(); d_y.pop_back();
        dz = d_z.back(); d_z.pop_back();
        tr = tp_r.back(); tp_r.pop_back();
        tg = tp_g.back(); tp_g.pop_back();
        tb = tp_b.back(); tp_b.pop_back();
        d = depth.back(); depth.pop_back();
        idx = pixel_index.back(); pixel_index.pop_back();
        return true;
    }

    bool is_empty() {
        return o_x.empty();
    }
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
    int num_samples;
    float focus_distance;
    float aperture_size;
    
    Mat4f transform;      // Forward transformation
    Mat4f inv_transform;  // Inverse transformation
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
    int depth[8] = {0, 0, 0, 0, 0, 0, 0, 0};  // Bounce depth for each ray
    int pixel_index[8] = {0, 1, 2, 3, 4, 5, 6, 7}; // Pixel index for each ray (0-7)
    fl time[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
};

struct HitRecord{
    Vec3f hit_point;
    Vec3f normal;
    fl t;
    int material_id;
};



//This one will hold the necessary information for mesh indexing, which will be helpful to keep the data oriented way of storing vertices along with the blas tlas structures 
//No need for SoA since will be accessed randomly//small size
//FOR REAL BASE MESHES




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
    
    std::vector<Mat4f> pl_transform;      // Forward transformation
    std::vector<Mat4f> pl_inv_transform;  // Inverse transformation
};

struct AreaLightData{
    std::vector<fl> al_pos_x;
    std::vector<fl> al_pos_y;
    std::vector<fl> al_pos_z;

    std::vector<fl> al_norm_x;
    std::vector<fl> al_norm_y;
    std::vector<fl> al_norm_z;

    std::vector<fl> size;
    std::vector<fl> area;
    std::vector<fl> al_intensity_r;
    std::vector<fl> al_intensity_g;
    std::vector<fl> al_intensity_b;
    

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

    std::vector<fl> tri_centro_x;

    std::vector<fl> tri_centro_y;

    std::vector<fl> tri_centro_z;

    std::vector<int> triangle_id;
    //should be -1 for mesh triangles
    std::vector<int> triangle_material_id;
};

struct PlaneData{
    std::vector<fl> plane_norm_x;
    std::vector<fl> plane_norm_y;
    std::vector<fl> plane_norm_z;

    std::vector<int> plane_point_vertex_id;
    
    std::vector<int> plane_id;

    std::vector<int> plane_material_id;
    
    std::vector<Mat4f> plane_transform;      // Forward transformation
    std::vector<Mat4f> inv_plane_transform;  // Inverse transformation
};
struct TransformationData{
    
    std::map<int, Mat4f> translations;
    std::map<int, Mat4f> scalings;
    std::map<int, Mat4f> rotations;
    std::map<int, Mat4f> composites;
    
};
struct MeshInfo{
    Mat4f transformation_matrix;
    
    Mat4f inverse_transformation_matrix;
    Vec3f mb_vector = {0.f,0.f,0.f};
    
    int id;
    int base_mesh_id;
    int base_triangle_index = -1;
    int triangle_count = 0;
    int material_id;
    //will be transcended by the first mesh
    int bvh_index=-1;
    bool reset_transform=false;
    bool motion_blur = false;
};

struct PrimitiveData{
    std::vector<PrimType> primitive_type;
    std::vector<int> prim_index;
};
// The main struct that holds all scene data, parsed from the JSON file.
struct Scene {
    // SoA (Structure of Arrays) for data-oriented design
    VertexData vertex_data__;
    PointLightData point_light_data__;
    AreaLightData area_light_data__;
    SphereData sphere_data__;
    TriangleData triangle_data__;
    PlaneData plane_data__;
    TransformationData transformation_data__;
    PrimitiveData prim_data__;
    
    
    
    std::vector<MeshInfo> mesh_data;
    
    
    std::vector<Camera> cameras;
    std::vector<Material> materials;
    
    // General scene settings
    Vec3f background_color;
    Vec3f ambient_light;
    fl shadow_ray_epsilon;
    fl intersection_test_epsilon;
    int max_recursion_depth;
    int primitive_count;
};


struct BVHNode{
    AABB bbox;
    int left_child_index;
    int first_prim,prim_count;
};

struct BVH{
    BVHNode* bvh_node_list;
    std::vector<int> prim_indices;
    int next_node=1;
};


struct BLAS{
    AABB bbox;
    Mat4f inv_transform;
    Mat4f transformation_matrix;
    Vec3f mb_vector;
    int material_id;
    //this also can denote the sphere index as we have a single blas struct for both mesh instance and the spheres
    int bvh_index;
    bool motion_blur = false;
};

struct TLASNode{
    AABB bbox;
    int left_child_index;
    int blas_id_0=-1,blas_id_1=-1;
};
//ONLY ONE WILL EXIST
//NO NEED FOR THE PRIMITIVE LIST AS THIS COVERS ALL OF THE BLASSES;
struct TLAS{
    TLASNode* tlas_root;
    int next_node;


};

//Structure for tlas blas: blas array. one tlas tree. many tlas nodes.


