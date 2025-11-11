#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <functional> // Required for std::function
#include <tuple>
// GEMINI created code parts
#include <map>
#include <cmath>
// GEMINI created code parts

#include "json.hpp" 
#include "types.hpp" 
#include "utils.hpp"
#include "cmath"
#include "miniply.h"
#include "bvh_tlas.hpp"

using json = nlohmann::json;

// --- Helper Functions (unchanged) ---
static Vec3f parseVec3f(const std::string& str) {
    Vec3f vec;
    std::stringstream ss(str);
    ss >> vec.x >> vec.y >> vec.z;
    return vec;
}

// Parses a string like "x y z w" into a Vec4f
static Vec4f parseVec4f(const std::string& str) {
    Vec4f vec;
    std::stringstream ss(str);
    ss >> vec.x >> vec.y >> vec.z >> vec.w;
    return vec;
}

// Parses a string like "x y" into a Vec2i
static Vec2i parseVec2i(const std::string& str) {
    Vec2i vec;
    std::stringstream ss(str);
    ss >> vec.x >> vec.y;
    return vec;
}
static std::string get_directory(const std::string& path) {
    size_t found = path.find_last_of("/\\");
    if (found != std::string::npos) {
        return path.substr(0, found); // e.g., "inputs/scene.json" -> "inputs"
    }
    return "."; // No directory part, use current directory
}

static Mat4f get_transform_by_id(const Scene& scene, char type, int id) {
    const auto& tf = scene.transformation_data__;
    switch (type) {
        case 't':
            if (tf.translations.count(id)) return tf.translations.at(id);
            break;
        case 's':
            if (tf.scalings.count(id)) return tf.scalings.at(id);
            break;
        case 'r':
            if (tf.rotations.count(id)) return tf.rotations.at(id);
            break;
        case 'c': // Composite
            if (tf.composites.count(id)) return tf.composites.at(id);
            break;
    }
    return create_identity_matrix();
}

static Mat4f build_composite_transform(const Scene& scene, const std::string& transform_str) {
    Mat4f composite_matrix = create_identity_matrix();
    std::stringstream ss(transform_str);
    std::string tf_id_str;
    
    std::vector<std::pair<char, int>> transforms;
    while (ss >> tf_id_str) {
        if (tf_id_str.empty()) continue;
        char type = tf_id_str[0];
        int id = std::stoi(tf_id_str.substr(1));
        transforms.push_back({type, id});
    }
    
    for (const auto& tf_pair : transforms) {
        Mat4f M_next = get_transform_by_id(scene, tf_pair.first, tf_pair.second);
        composite_matrix = M_next * composite_matrix;
    }
    
    return composite_matrix;
}
// GEMINI created code parts
class Parser {
public:
    static Scene parseScene(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + filepath);
        }
        json data;
        file >> data;

        Scene scene;
        const auto& sceneData = data["Scene"];

        // --- Parse General Scene Settings ---
        scene.background_color = parseVec3f(sceneData.contains("BackgroundColor") ? 
            sceneData.at("BackgroundColor").get<std::string>() : "0 0 0");
        scene.shadow_ray_epsilon = std::stof(sceneData.contains("ShadowRayEpsilon") ? 
            sceneData.at("ShadowRayEpsilon").get<std::string>() : "1e-3");
        scene.intersection_test_epsilon = std::stof(sceneData.contains("IntersectionTestEpsilon") ? 
            sceneData.at("IntersectionTestEpsilon").get<std::string>() : "1e-9");
        
        // **FIXED**: Correctly parse MaxRecursionDepth from a string.
        scene.max_recursion_depth = std::stoi(sceneData.contains("MaxRecursionDepth") ? 
            sceneData.at("MaxRecursionDepth").get<std::string>() : "0");
        
        // Handle Lights section safely
        if (sceneData.contains("Lights")) {
            const auto& lightsData = sceneData["Lights"];
            scene.ambient_light = parseVec3f(lightsData.contains("AmbientLight") ? 
                lightsData.at("AmbientLight").get<std::string>() : "0 0 0");
        } else {
            scene.ambient_light = {0, 0, 0};
        }
        
        // --- Parse Vertex Data ---
        if (sceneData.contains("VertexData")) {
            const auto& vertexData = sceneData.at("VertexData");
            std::string vertexDataStr;
            
            // Handle both formats: direct string or object with _data field
            if (vertexData.is_string()) {
                vertexDataStr = vertexData.get<std::string>();
            } else if (vertexData.is_object() && vertexData.contains("_data")) {
                vertexDataStr = vertexData.at("_data").get<std::string>();
            }
            
            if (!vertexDataStr.empty()) {
                std::stringstream vertexStream(vertexDataStr);
                fl x, y, z;
                while (vertexStream >> x >> y >> z) {
                    // Directly populate SoA structure
                    scene.vertex_data__.v_pos_x.push_back(x);
                    scene.vertex_data__.v_pos_y.push_back(y);
                    scene.vertex_data__.v_pos_z.push_back(z);
                    
                    // Initialize normals to zero (will be accumulated later)
                    scene.vertex_data__.v_nor_x.push_back(0.0f);
                    scene.vertex_data__.v_nor_y.push_back(0.0f);
                    scene.vertex_data__.v_nor_z.push_back(0.0f);
                }
            }
        }

        
        auto processOneOrMany = [](const json& node, std::function<void(const json&)> parserFunc) {
            if (node.is_array()) {
                for (const auto& item : node) {
                    parserFunc(item);
                }
            } else if (node.is_object()) {
                parserFunc(node);
            }
        };

        // --- Parse Materials (Generalized) ---
        auto parseMaterial = [&](const json& mat) {
            Material material;
            material.id = std::stoi(mat.at("_id").get<std::string>());
            material.type = mat.contains("_type") ? mat.at("_type").get<std::string>() : "default";
            material.ambient_reflectance = parseVec3f(mat.at("AmbientReflectance").get<std::string>());
            material.diffuse_reflectance = parseVec3f(mat.at("DiffuseReflectance").get<std::string>());
            material.specular_reflectance = parseVec3f(mat.at("SpecularReflectance").get<std::string>());
            material.phong_exponent = mat.contains("PhongExponent") ? 
                std::stof(mat.at("PhongExponent").get<std::string>()) : 1.0f;
            if (mat.contains("MirrorReflectance")) material.mirror_reflectance = parseVec3f(mat.at("MirrorReflectance").get<std::string>());
            if (mat.contains("AbsorptionCoefficient")) material.absorption_coefficient = parseVec3f(mat.at("AbsorptionCoefficient").get<std::string>());
            if (mat.contains("RefractionIndex")) material.refraction_index = std::stof(mat.at("RefractionIndex").get<std::string>());
            if (mat.contains("AbsorptionIndex")) material.absorption_index = std::stof(mat.at("AbsorptionIndex").get<std::string>());
            scene.materials.push_back(material);
        };
        if (sceneData.contains("Materials") && sceneData["Materials"].contains("Material")) {
            processOneOrMany(sceneData.at("Materials").at("Material"), parseMaterial);
        }
       
        auto parseTranslation = [&](const json& trans) {
            int id = std::stoi(trans.at("_id").get<std::string>());
            Vec3f t = parseVec3f(trans.at("_data").get<std::string>());
            scene.transformation_data__.translations[id] = create_translation_matrix(t);
        };

        auto parseRotation = [&](const json& rot) {
            int id = std::stoi(rot.at("_id").get<std::string>());
            std::string data_str = rot.at("_data").get<std::string>();
            std::stringstream ss(data_str);
            fl angle;
            Vec3f axis;
            ss >> angle >> axis.x >> axis.y >> axis.z;
            scene.transformation_data__.rotations[id] = create_rotation_matrix(angle, axis);
        };

        auto parseScaling = [&](const json& scale) {
            int id = std::stoi(scale.at("_id").get<std::string>());
            Vec3f s = parseVec3f(scale.at("_data").get<std::string>());
            scene.transformation_data__.scalings[id] = create_scaling_matrix(s);
        };

        auto parseComposite = [&](const json& comp) {
            int id = std::stoi(comp.at("_id").get<std::string>());
            std::string data_str = comp.at("_data").get<std::string>();
            std::stringstream ss(data_str);
            Mat4f mat;
            ss >> mat.m11 >> mat.m12 >> mat.m13 >> mat.m14
                >> mat.m21 >> mat.m22 >> mat.m23 >> mat.m24
                >> mat.m31 >> mat.m32 >> mat.m33 >> mat.m34
                >> mat.m41 >> mat.m42 >> mat.m43 >> mat.m44;
            scene.transformation_data__.composites[id] = mat;
        };

        if(sceneData.contains("Transformations")){
            const auto& transforms = sceneData.at("Transformations");
            if(transforms.contains("Translation")){
                processOneOrMany(transforms.at("Translation"),parseTranslation);
            }
                if(transforms.contains("Rotation")){
                processOneOrMany(transforms.at("Rotation"),parseRotation);
            }
                if(transforms.contains("Scaling")){
                processOneOrMany(transforms.at("Scaling"),parseScaling);
            }
            if(transforms.contains("Composite")){
                processOneOrMany(transforms.at("Composite"),parseComposite);
            }
        }
       
        // --- Parse Lights (Generalized) ---
        auto parsePointLight = [&](const json& light) {
            int id = std::stoi(light.at("_id").get<std::string>());
            
            // Parse position directly into scalars
            fl pos_x, pos_y, pos_z;
            std::stringstream pos_stream(light.at("Position").get<std::string>());
            pos_stream >> pos_x >> pos_y >> pos_z;
            
            // Parse intensity directly into scalars
            fl int_r, int_g, int_b;
            std::stringstream int_stream(light.at("Intensity").get<std::string>());
            int_stream >> int_r >> int_g >> int_b;
            
            // Directly populate SoA structure
            scene.point_light_data__.pl_id.push_back(id);
            scene.point_light_data__.pl_pos_x.push_back(pos_x);
            scene.point_light_data__.pl_pos_y.push_back(pos_y);
            scene.point_light_data__.pl_pos_z.push_back(pos_z);
            scene.point_light_data__.pl_intensity_r.push_back(int_r);
            scene.point_light_data__.pl_intensity_g.push_back(int_g);
            scene.point_light_data__.pl_intensity_b.push_back(int_b);
        };
        if (sceneData.contains("Lights") && sceneData["Lights"].contains("PointLight")) {
            processOneOrMany(sceneData.at("Lights").at("PointLight"), parsePointLight);
        }

        
        //Calculate Normals 1- for faces(per vertex, per face, )
       
            

        
        auto parseCamera = [&](const json& cam) {
            Camera camera;
            camera.id = std::stoi(cam.at("_id").get<std::string>());
            camera.image_name = cam.at("ImageName").get<std::string>();
            camera.near_distance = std::stof(cam.at("NearDistance").get<std::string>());
            camera.image_resolution = parseVec2i(cam.at("ImageResolution").get<std::string>());
            camera.position = parseVec3f(cam.at("Position").get<std::string>());
            
            // --- Gaze Calculation ---
            if (cam.contains("Gaze")) {
                camera.gaze = parseVec3f(cam.at("Gaze").get<std::string>());
            } else if (cam.contains("GazePoint")) {
                // Calculate gaze direction from position to gaze point
                Vec3f gazePoint = parseVec3f(cam.at("GazePoint").get<std::string>());
                camera.gaze = {
                    gazePoint.x - camera.position.x,
                    gazePoint.y - camera.position.y,
                    gazePoint.z - camera.position.z
                };
            } else {
                // Default gaze direction
                camera.gaze = {0.0f, 0.0f, -1.0f};
            }
            // Normalize Gaze (Gaze / ||Gaze||)
            normalize_on_place(camera.gaze);
            
            // --- Up Vector Correction (as per homework spec) ---
            // This ensures that the final up vector ('v') is perpendicular to the gaze vector ('w').
            camera.up = parseVec3f(cam.at("Up").get<std::string>());
            
            Vec3f w = -camera.gaze; // w = -Gaze/||Gaze||
            Vec3f u = cross(camera.up, w); // u = v0 x w
            normalize_on_place(u);
            
            camera.up = cross(w, u); // v = w x u (this is the new, corrected up vector)
            // camera.up is now normalized because w and u are orthonormal.

            // --- Near Plane Calculation ---
            std::string cam_type = cam.contains("_type") ? cam.at("_type").get<std::string>() : "default";
            
            if (cam_type == "lookAt") {
                std::cout<<"LOOKAt"<<std::endl;
                // Calculate near plane from FovY and Aspect Ratio
                fl fovY_degrees = std::stof(cam.at("FovY").get<std::string>());
                fl fovY_rad = fovY_degrees * (M_PI / 180.0f);
                
                fl top = camera.near_distance * std::tan(fovY_rad / 2.0f);
                fl bottom = -top;
                
                fl aspect_ratio = (fl)camera.image_resolution.x / (fl)camera.image_resolution.y;
                fl right = top * aspect_ratio;
                fl left = -right;
                
                // Vec4f { left, right, bottom, top }
                camera.near_plane = {left, right, bottom, top};
                
            } else if (cam.contains("NearPlane")) {
                // Handle non-lookAt cameras that specify a plane
                camera.near_plane = parseVec4f(cam.at("NearPlane").get<std::string>());
            } else {
                // Default near plane if not specified (and not lookAt)
                camera.near_plane = {-1.0f, 1.0f, -1.0f, 1.0f};
            }

            scene.cameras.push_back(camera);
        };
        if (sceneData.contains("Cameras") && sceneData["Cameras"].contains("Camera")) {
            processOneOrMany(sceneData.at("Cameras").at("Camera"), parseCamera);
        }
       
        std::map<int, MeshInfo> mesh_info_map; 
       

        // --- Parse Objects (Generalized) ---
        if (sceneData.contains("Objects")) {
            const auto& objects = sceneData.at("Objects");

            if (objects.contains("Sphere")) {
                auto parseSphere = [&](const json& obj) {
                    BLAS blas;
                    int id = std::stoi(obj.at("_id").get<std::string>());
                    int material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    int center_vertex_id = std::stoi(obj.at("Center").get<std::string>()) - 1;
                    fl radius = std::stof(obj.at("Radius").get<std::string>());
                    
                    int sphere_index = scene.sphere_data__.sphere_id.size();
                    
                    scene.prim_data__.prim_index.push_back(sphere_index);
                    scene.prim_data__.primitive_type.push_back(SPHERE);
                    // Directly populate SoA structure
                    scene.sphere_data__.sphere_id.push_back(id);
                    scene.sphere_data__.sphere_mat_id.push_back(material_id);
                    scene.sphere_data__.sphere_center_vertex_id.push_back(center_vertex_id);
                    scene.sphere_data__.sphere_radius_sq.push_back(radius*radius);
                    
                    // Parse transformations if they exist
                    Mat4f transformation_matrix = create_identity_matrix();
                    if (obj.contains("Transformations")) {
                        transformation_matrix = build_composite_transform(scene, obj.at("Transformations").get<std::string>());
                    }
                    
                    // Create BLAS for this sphere
                    // Get sphere center position (untransformed)
                    fl center_x = scene.vertex_data__.v_pos_x[center_vertex_id];
                    fl center_y = scene.vertex_data__.v_pos_y[center_vertex_id];
                    fl center_z = scene.vertex_data__.v_pos_z[center_vertex_id];
                    
                    
                    
                    // Create AABB for sphere in local space with scaled radius
                    AABB local_bbox;
                    local_bbox.bbox_xmin = center_x - radius;
                    local_bbox.bbox_ymin = center_y - radius;
                    local_bbox.bbox_zmin = center_z - radius;
                    local_bbox.bbox_xmax = center_x + radius;
                    local_bbox.bbox_ymax = center_y + radius;
                    local_bbox.bbox_zmax = center_z + radius;
                    
                    // Transform the AABB to world space
                    blas.bbox = transform_aabb(local_bbox, transformation_matrix);
                    
                    // Store transforms for sphere
                    blas.transformation_matrix = transformation_matrix;
                    blas.inv_transform = mat_inv(transformation_matrix);
                    blas.material_id = material_id;
                    blas.bvh_index = 10000 + sphere_index; // Magic number for sphere identification
                    
                    blas_list.push_back(blas);
                };
                processOneOrMany(objects.at("Sphere"), parseSphere);
            }
            
            if (objects.contains("Triangle")) {
                auto parseTriangle = [&](const json& obj) {
                    int id = std::stoi(obj.at("_id").get<std::string>());
                    int material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    
                    // Parse indices directly
                    int v0, v1, v2;
                    std::stringstream ind_stream(obj.at("Indices").get<std::string>());
                    ind_stream >> v0 >> v1 >> v2;
                    v0--; v1--; v2--;  // Convert to 0-based indexing
                    
                    // Get vertex positions from SoA
                    fl v0_x = scene.vertex_data__.v_pos_x[v0];
                    fl v0_y = scene.vertex_data__.v_pos_y[v0];
                    fl v0_z = scene.vertex_data__.v_pos_z[v0];
                    
                    fl v1_x = scene.vertex_data__.v_pos_x[v1];
                    fl v1_y = scene.vertex_data__.v_pos_y[v1];
                    fl v1_z = scene.vertex_data__.v_pos_z[v1];
                    
                    fl v2_x = scene.vertex_data__.v_pos_x[v2];
                    fl v2_y = scene.vertex_data__.v_pos_y[v2];
                    fl v2_z = scene.vertex_data__.v_pos_z[v2];
                    
                    // Calculate edges using scalar operations
                    fl edge1_x, edge1_y, edge1_z;
                    fl edge2_x, edge2_y, edge2_z;
                    subtract_scalar(v2_x, v2_y, v2_z, v1_x, v1_y, v1_z, edge1_x, edge1_y, edge1_z);
                    subtract_scalar(v0_x, v0_y, v0_z, v1_x, v1_y, v1_z, edge2_x, edge2_y, edge2_z);
                    
                    // Calculate and normalize normal
                    fl norm_x, norm_y, norm_z;
                    cross_scalar(edge1_x, edge1_y, edge1_z, edge2_x, edge2_y, edge2_z, norm_x, norm_y, norm_z);
                    normalize_scalar(norm_x, norm_y, norm_z);
                    
                    fl centro_x = (v0_x + v1_x+ v2_x)/3;
                    fl centro_y = (v0_y + v1_y+ v2_y)/3;
                    fl centro_z = (v0_z + v1_z+ v2_z)/3;
                
                    scene.prim_data__.prim_index.push_back(scene.triangle_data__.v0_ind.size());
                    scene.prim_data__.primitive_type.push_back(TRIANGLE);

                    // Directly populate SoA structure
                    scene.triangle_data__.v0_ind.push_back(v0);
                    scene.triangle_data__.v1_ind.push_back(v1);
                    scene.triangle_data__.v2_ind.push_back(v2);
                    scene.triangle_data__.tri_norm_x.push_back(norm_x);
                    scene.triangle_data__.tri_norm_y.push_back(norm_y);
                    scene.triangle_data__.tri_norm_z.push_back(norm_z);
                    scene.triangle_data__.triangle_id.push_back(id);
                    scene.triangle_data__.triangle_material_id.push_back(material_id);

                    scene.triangle_data__.tri_centro_y.push_back(centro_y);

                    scene.triangle_data__.tri_centro_z.push_back(centro_z);

                    scene.triangle_data__.tri_centro_x.push_back(centro_x);

                };
                processOneOrMany(objects.at("Triangle"), parseTriangle);
            }
            

            //HW1: Parse Mesh - push all faces directly into TriangleData
            //HW2: Parse Mesh - push the data exactly as before, create a mesh info instance which will hold info required for instancing such as id and base triangle index.
            if (objects.contains("Mesh")) {
                auto parseMesh = [&](const json& obj) {
                    BVH bvh;
                    BLAS blas;
                    MeshInfo mesh_info;
                    mesh_info.bvh_index = bvh_list.size();
                    mesh_info.id = std::stoi(obj.at("_id").get<std::string>());
                    mesh_info.material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    mesh_info.base_mesh_id = -1; // This is a base mesh
                    mesh_info.reset_transform = false;
                    
                    mesh_info.base_triangle_index = scene.triangle_data__.v0_ind.size();
                    int start_triangle_count = mesh_info.base_triangle_index;

                    if (obj.contains("Transformations")) {
                        mesh_info.transformation_matrix = build_composite_transform(scene, obj.at("Transformations").get<std::string>());
                    } else {
                        mesh_info.transformation_matrix = create_identity_matrix();
                    }
                    
                    mesh_info.inverse_transformation_matrix = mat_inv(mesh_info.transformation_matrix);
                 
                    
                    int mesh_id = mesh_info.id;
                    int material_id = mesh_info.material_id;
                   
                    
                    if (obj.contains("_shadingMode") && obj.at("_shadingMode").get<std::string>() == "smooth") {
                        G_SMOOTH_SHADING_ENABLED = true;
                    }
                    
                    const auto& facesData = obj.at("Faces");
                    
                    // --- NEW LOGIC: Handle _plyFile ---
                    if (facesData.is_object() && facesData.contains("_plyFile")) {
                        std::string ply_relative_path = facesData.at("_plyFile").get<std::string>();

                        // --- MODIFICATION: Resolve path relative to JSON file ---
                        std::string json_dir = get_directory(filepath); // filepath is from parseScene
                        std::string ply_filename = json_dir + "/" + ply_relative_path;
                        
                        // 1. Get the current vertex count as our offset
                        size_t vertex_offset = scene.vertex_data__.v_pos_x.size();

                        miniply::PLYReader reader(ply_filename.c_str());
                        if (!reader.valid()) {
                            throw std::runtime_error("Could not open or parse PLY file: " + ply_filename);
                        }

                        // --- 2. Load Vertices from PLY ---
                        std::vector<float> ply_verts;
                        uint32_t num_ply_verts = 0;

                        while (reader.has_element()) {
                            if (reader.element_is(miniply::kPLYVertexElement)) {
                                if (!reader.load_element()) {
                                    throw std::runtime_error("Could not load vertex element from " + ply_filename);
                                }
                                
                                num_ply_verts = reader.num_rows();
                                if (num_ply_verts == 0) break; // No vertices

                                uint32_t pos_props[3];
                                if (!reader.find_pos(pos_props)) {
                                     throw std::runtime_error("No 'x', 'y', 'z' properties in " + ply_filename);
                                }

                                // Extract vertex data into a temporary buffer
                                ply_verts.resize(num_ply_verts * 3);
                                reader.extract_properties(pos_props, 3, miniply::PLYPropertyType::Float, ply_verts.data());
                                
                                // Reserve space in scene data for efficiency
                                scene.vertex_data__.v_pos_x.reserve(vertex_offset + num_ply_verts);
                                scene.vertex_data__.v_pos_y.reserve(vertex_offset + num_ply_verts);
                                scene.vertex_data__.v_pos_z.reserve(vertex_offset + num_ply_verts);
                                scene.vertex_data__.v_nor_x.reserve(vertex_offset + num_ply_verts);
                                scene.vertex_data__.v_nor_y.reserve(vertex_offset + num_ply_verts);
                                scene.vertex_data__.v_nor_z.reserve(vertex_offset + num_ply_verts);

                                // Append new vertices and initialize their normals to zero
                                for (size_t i = 0; i < num_ply_verts; ++i) {
                                    scene.vertex_data__.v_pos_x.push_back(ply_verts[i * 3 + 0]);
                                    scene.vertex_data__.v_pos_y.push_back(ply_verts[i * 3 + 1]);
                                    scene.vertex_data__.v_pos_z.push_back(ply_verts[i * 3 + 2]);
                                    scene.vertex_data__.v_nor_x.push_back(0.0f);
                                    scene.vertex_data__.v_nor_y.push_back(0.0f);
                                    scene.vertex_data__.v_nor_z.push_back(0.0f);
                                }
                                break; // Found vertex element, stop loop
                            }
                            reader.next_element();
                        }
                        
                        if (num_ply_verts == 0) {
                            std::cerr << "Warning: PLY file " << ply_filename << " has no vertices. Skipping mesh." << std::endl;
                            return; // Use 'return' since we are in a lambda
                        }

                        // --- 3. Load Faces from PLY ---
                        // We must reset the reader to find the face element
                        miniply::PLYReader face_reader(ply_filename.c_str());
                        while (face_reader.has_element()) {
                            if (face_reader.element_is(miniply::kPLYFaceElement)) {
                                if (!face_reader.load_element()) {
                                    throw std::runtime_error("Could not load face element from " + ply_filename);
                                }

                                uint32_t face_props[1];
                                // Find "vertex_indices" or "vertex_index"
                                if (!face_reader.find_indices(face_props)) {
                                    throw std::runtime_error("No 'vertex_indices' property in " + ply_filename);
                                }
                                
                                uint32_t prop_idx = face_props[0];
                                uint32_t num_tris = face_reader.num_triangles(prop_idx);
                                if (num_tris == 0) break; // No triangles

                                std::vector<int> face_indices(num_tris * 3);
                                
                                // This function handles triangulation of quads/polygons
                                face_reader.extract_triangles(prop_idx, ply_verts.data(), num_ply_verts, 
                                                              miniply::PLYPropertyType::Int, face_indices.data());

                                // --- 4. Create Triangles & Accumulate Normals ---
                                int triangle_counter = 0;
                                for (size_t i = 0; i < face_indices.size(); i += 3) {
                                    // Apply offset. PLY indices are 0-based.
                                    int index_v0 = face_indices[i + 0] + vertex_offset;
                                    int index_v1 = face_indices[i + 1] + vertex_offset;
                                    int index_v2 = face_indices[i + 2] + vertex_offset;

                                    // --- RE-USED LOGIC from your _data parser ---
                                    fl v0_x = scene.vertex_data__.v_pos_x[index_v0];
                                    fl v0_y = scene.vertex_data__.v_pos_y[index_v0];
                                    fl v0_z = scene.vertex_data__.v_pos_z[index_v0];
                                    
                                    fl v1_x = scene.vertex_data__.v_pos_x[index_v1];
                                    fl v1_y = scene.vertex_data__.v_pos_y[index_v1];
                                    fl v1_z = scene.vertex_data__.v_pos_z[index_v1];
                                    
                                    fl v2_x = scene.vertex_data__.v_pos_x[index_v2];
                                    fl v2_y = scene.vertex_data__.v_pos_y[index_v2];
                                    fl v2_z = scene.vertex_data__.v_pos_z[index_v2];
                                    
                                    fl edge1_x, edge1_y, edge1_z;
                                    fl edge2_x, edge2_y, edge2_z;
                                    subtract_scalar(v1_x, v1_y, v1_z, v0_x, v0_y, v0_z, edge1_x, edge1_y, edge1_z);
                                    subtract_scalar(v2_x, v2_y, v2_z, v0_x, v0_y, v0_z, edge2_x, edge2_y, edge2_z);
                                    
                                    fl face_norm_x, face_norm_y, face_norm_z;
                                    cross_scalar(edge1_x, edge1_y, edge1_z, edge2_x, edge2_y, edge2_z, 
                                               face_norm_x, face_norm_y, face_norm_z);
                                    
                                    fl norm_x = face_norm_x, norm_y = face_norm_y, norm_z = face_norm_z;
                                    normalize_scalar(norm_x, norm_y, norm_z);
                                    
                                    fl centro_x = (v0_x + v1_x+ v2_x)/3;
                                    fl centro_y = (v0_y + v1_y+ v2_y)/3;
                                    fl centro_z = (v0_z + v1_z+ v2_z)/3;
                                
                                    scene.prim_data__.prim_index.push_back(scene.triangle_data__.v0_ind.size());
                                    scene.prim_data__.primitive_type.push_back(TRIANGLE);
                                    
                                    bvh.prim_indices.push_back(scene.triangle_data__.v0_ind.size());
                                    scene.triangle_data__.v0_ind.push_back(index_v0);
                                    scene.triangle_data__.v1_ind.push_back(index_v1);
                                    scene.triangle_data__.v2_ind.push_back(index_v2);
                                    scene.triangle_data__.tri_norm_x.push_back(norm_x);
                                    scene.triangle_data__.tri_norm_y.push_back(norm_y);
                                    scene.triangle_data__.tri_norm_z.push_back(norm_z);
                                    scene.triangle_data__.triangle_id.push_back(mesh_id * 1000000 + triangle_counter);
                                    scene.triangle_data__.triangle_material_id.push_back(material_id);
                                    
                                    scene.triangle_data__.tri_centro_y.push_back(centro_y);

                                    scene.triangle_data__.tri_centro_z.push_back(centro_z);

                                    scene.triangle_data__.tri_centro_x.push_back(centro_x);
                                    
                                    // Accumulate area-weighted normals
                                    scene.vertex_data__.v_nor_x[index_v0] += face_norm_x;
                                    scene.vertex_data__.v_nor_y[index_v0] += face_norm_y;
                                    scene.vertex_data__.v_nor_z[index_v0] += face_norm_z;
                                    
                                    scene.vertex_data__.v_nor_x[index_v1] += face_norm_x;
                                    scene.vertex_data__.v_nor_y[index_v1] += face_norm_y;
                                    scene.vertex_data__.v_nor_z[index_v1] += face_norm_z;
                                    
                                    scene.vertex_data__.v_nor_x[index_v2] += face_norm_x;
                                    scene.vertex_data__.v_nor_y[index_v2] += face_norm_y;
                                    scene.vertex_data__.v_nor_z[index_v2] += face_norm_z;
                                    
                                    triangle_counter++;
                                }
                                break; // Found face element, stop loop
                            }
                            face_reader.next_element();
                        }
                    }
                    // --- ORIGINAL LOGIC: Handle _data or string ---
                    else {
                        std::string facesDataStr;
                        if (facesData.is_string()) {
                            facesDataStr = facesData.get<std::string>();
                        } else if (facesData.is_object() && facesData.contains("_data")) {
                            facesDataStr = facesData.at("_data").get<std::string>();
                        }
                        
                        if (!facesDataStr.empty()) {
                            std::stringstream faceStream(facesDataStr);
                            int v0, v1, v2;
                            int triangle_counter = 0;
                            while (faceStream >> v0 >> v1 >> v2) {
                                // Convert from 1-based (JSON) to 0-based indexing
                                int index_v0 = v0 - 1, index_v1 = v1 - 1, index_v2 = v2 - 1;
                                
                                fl v0_x = scene.vertex_data__.v_pos_x[index_v0];
                                fl v0_y = scene.vertex_data__.v_pos_y[index_v0];
                                fl v0_z = scene.vertex_data__.v_pos_z[index_v0];
                                
                                fl v1_x = scene.vertex_data__.v_pos_x[index_v1];
                                fl v1_y = scene.vertex_data__.v_pos_y[index_v1];
                                fl v1_z = scene.vertex_data__.v_pos_z[index_v1];
                                
                                fl v2_x = scene.vertex_data__.v_pos_x[index_v2];
                                fl v2_y = scene.vertex_data__.v_pos_y[index_v2];
                                fl v2_z = scene.vertex_data__.v_pos_z[index_v2];
                                
                                fl edge1_x, edge1_y, edge1_z;
                                fl edge2_x, edge2_y, edge2_z;
                                subtract_scalar(v1_x, v1_y, v1_z, v0_x, v0_y, v0_z, edge1_x, edge1_y, edge1_z);
                                subtract_scalar(v2_x, v2_y, v2_z, v0_x, v0_y, v0_z, edge2_x, edge2_y, edge2_z);
                                
                                fl face_norm_x, face_norm_y, face_norm_z;
                                cross_scalar(edge1_x, edge1_y, edge1_z, edge2_x, edge2_y, edge2_z, 
                                           face_norm_x, face_norm_y, face_norm_z);
                                
                                fl norm_x = face_norm_x, norm_y = face_norm_y, norm_z = face_norm_z;
                                normalize_scalar(norm_x, norm_y, norm_z);

                                fl centro_x = (v0_x + v1_x+ v2_x)/3;
                                fl centro_y = (v0_y + v1_y+ v2_y)/3;
                                fl centro_z = (v0_z + v1_z+ v2_z)/3;
                                
                                scene.prim_data__.prim_index.push_back(scene.triangle_data__.v0_ind.size());
                                scene.prim_data__.primitive_type.push_back(TRIANGLE);
                                
                                bvh.prim_indices.push_back(scene.triangle_data__.v0_ind.size());
                                scene.triangle_data__.v0_ind.push_back(index_v0);
                                scene.triangle_data__.v1_ind.push_back(index_v1);
                                scene.triangle_data__.v2_ind.push_back(index_v2);
                                scene.triangle_data__.tri_norm_x.push_back(norm_x);
                                scene.triangle_data__.tri_norm_y.push_back(norm_y);
                                scene.triangle_data__.tri_norm_z.push_back(norm_z);
                                
                                scene.triangle_data__.tri_centro_y.push_back(centro_y);

                                scene.triangle_data__.tri_centro_z.push_back(centro_z);

                                scene.triangle_data__.tri_centro_x.push_back(centro_x);
                                
                                scene.triangle_data__.triangle_id.push_back(mesh_id * 1000000 + triangle_counter);
                                scene.triangle_data__.triangle_material_id.push_back(material_id);
                                
                                // Accumulate area-weighted normals
                                scene.vertex_data__.v_nor_x[index_v0] += face_norm_x;
                                scene.vertex_data__.v_nor_y[index_v0] += face_norm_y;
                                scene.vertex_data__.v_nor_z[index_v0] += face_norm_z;
                                
                                scene.vertex_data__.v_nor_x[index_v1] += face_norm_x;
                                scene.vertex_data__.v_nor_y[index_v1] += face_norm_y;
                                scene.vertex_data__.v_nor_z[index_v1] += face_norm_z;
                                
                                scene.vertex_data__.v_nor_x[index_v2] += face_norm_x;
                                scene.vertex_data__.v_nor_y[index_v2] += face_norm_y;
                                scene.vertex_data__.v_nor_z[index_v2] += face_norm_z;
                                
                                triangle_counter++;
                            }
                        }
                    }
                   
                    int end_triangle_count = scene.triangle_data__.v0_ind.size();
                    mesh_info.triangle_count = end_triangle_count - start_triangle_count;
                    
                    scene.mesh_data.push_back(mesh_info);
                    mesh_info_map[mesh_info.id] = mesh_info;
                    construct_BVH(mesh_info,bvh,scene);
                    bvh_list.push_back(bvh);
                    construct_BLAS(mesh_info,blas,scene);
                    blas_list.push_back(blas);
                };
                processOneOrMany(objects.at("Mesh"), parseMesh);
            }
           
            if (objects.contains("MeshInstance")) {
                auto parseMeshInstance = [&](const json& obj) {
                    BLAS blas;
                    MeshInfo instance_info;
                    instance_info.id = std::stoi(obj.at("_id").get<std::string>());
                    instance_info.material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    instance_info.base_mesh_id = std::stoi(obj.at("_baseMeshId").get<std::string>());
                    instance_info.reset_transform = obj.contains("_resetTransform") ? 
                        (obj.at("_resetTransform").get<std::string>() == "true") : false;

                    if (mesh_info_map.find(instance_info.base_mesh_id) == mesh_info_map.end()) {
                        throw std::runtime_error("Could not find base mesh with id: " + std::to_string(instance_info.base_mesh_id));
                    }
                    const MeshInfo& base_mesh_info = mesh_info_map.at(instance_info.base_mesh_id);

                    instance_info.base_triangle_index = base_mesh_info.base_triangle_index;
                    instance_info.triangle_count = base_mesh_info.triangle_count;
                    instance_info.bvh_index = base_mesh_info.bvh_index;

                    Mat4f M_local = create_identity_matrix();
                    if (obj.contains("Transformations")) {
                        M_local = build_composite_transform(scene, obj.at("Transformations").get<std::string>());
                    }

                    if (instance_info.reset_transform) {
                        instance_info.transformation_matrix = M_local;
                    } else {
                        instance_info.transformation_matrix = M_local * base_mesh_info.transformation_matrix;
                    }
                    
                    instance_info.inverse_transformation_matrix = mat_inv(instance_info.transformation_matrix);
                    

                    scene.mesh_data.push_back(instance_info);
                    mesh_info_map[instance_info.id] = instance_info;
                    construct_BLAS(instance_info,blas,scene);
                    blas_list.push_back(blas);
                    
                };
                processOneOrMany(objects.at("MeshInstance"), parseMeshInstance);
            }
           
            
            if (objects.contains("Plane")) {
                auto parsePlane = [&](const json& obj) {
                    // Parse all data into local variables first
                    int id = std::stoi(obj.at("_id").get<std::string>());
                    int material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    int point_vertex_id = std::stoi(obj.at("Point").get<std::string>()) - 1;
                    Vec3f normal = parseVec3f(obj.at("Normal").get<std::string>());

                    // Directly populate the SoA structure (scene.plane_data__)
                    scene.plane_data__.plane_id.push_back(id);
                    scene.plane_data__.plane_material_id.push_back(material_id);
                    scene.plane_data__.plane_point_vertex_id.push_back(point_vertex_id);
                    
                    scene.plane_data__.plane_norm_x.push_back(normal.x);
                    scene.plane_data__.plane_norm_y.push_back(normal.y);
                    scene.plane_data__.plane_norm_z.push_back(normal.z);
                };
                processOneOrMany(objects.at("Plane"), parsePlane);
            }
        }
        
        // Normalize all accumulated vertex normals in SoA
        size_t vertex_count = scene.vertex_data__.v_pos_x.size();
        for(size_t i = 0; i < vertex_count; ++i){
            normalize_scalar(scene.vertex_data__.v_nor_x[i], 
                           scene.vertex_data__.v_nor_y[i], 
                           scene.vertex_data__.v_nor_z[i]);
        }
        
       
        scene.transformation_data__.translations.clear();
        scene.transformation_data__.scalings.clear();
        scene.transformation_data__.rotations.clear();
        scene.transformation_data__.composites.clear();
       
        scene.primitive_count = scene.triangle_data__.v0_ind.size()+scene.sphere_data__.sphere_center_vertex_id.size();
        return scene;
    }
};