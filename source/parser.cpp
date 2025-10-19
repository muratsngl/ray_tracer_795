#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <functional> // Required for std::function

#include "json.hpp" 
#include "types.hpp" 
#include "utils.hpp"
#include "cmath"

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

// Parses a string like "v0 v1 v2" into a Vec3i, converting from 1-based to 0-based indices
static Vec3i parseIndices(const std::string& str) {
    Vec3i vec;
    std::stringstream ss(str);
    ss >> vec.v0 >> vec.v1 >> vec.v2;
    // The file format uses 1-based indexing, so we convert to 0-based.
    vec.v0--; 
    vec.v1--; 
    vec.v2--;
    return vec;
}


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
            sceneData.at("IntersectionTestEpsilon").get<std::string>() : "1e-6");
        
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
                    scene.vertex_data__.v_pos_x.push_back(x);

                    scene.vertex_data__.v_pos_y.push_back(y);

                    scene.vertex_data__.v_pos_z.push_back(z);
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

        // --- Parse Lights (Generalized) ---
        auto parsePointLight = [&](const json& light) {
            PointLight pl;
            pl.id = std::stoi(light.at("_id").get<std::string>());
            pl.position = parseVec3f(light.at("Position").get<std::string>());
            pl.intensity = parseVec3f(light.at("Intensity").get<std::string>());
            scene.point_light_data__.pl_id.push_back(pl.id);
            scene.point_light_data__.pl_pos_x.push_back(pl.position.x);

            scene.point_light_data__.pl_pos_y.push_back(pl.position.y);

            scene.point_light_data__.pl_pos_z.push_back(pl.position.z);
            
            scene.point_light_data__.pl_intensity_r.push_back(pl.intensity.x);
            
            scene.point_light_data__.pl_intensity_g.push_back(pl.intensity.y);
            
            scene.point_light_data__.pl_intensity_b.push_back(pl.intensity.z);
        };
        if (sceneData.contains("Lights") && sceneData["Lights"].contains("PointLight")) {
            processOneOrMany(sceneData.at("Lights").at("PointLight"), parsePointLight);
        }

        
        //Calculate Normals 1- for faces(per vertex, per face, )
       
            

        
        // --- Parse Cameras (Generalized) ---
            auto parseCamera = [&](const json& cam) {
            Camera camera;
            camera.id = std::stoi(cam.at("_id").get<std::string>());
            camera.image_name = cam.at("ImageName").get<std::string>();
            camera.near_distance = std::stof(cam.at("NearDistance").get<std::string>());
            camera.image_resolution = parseVec2i(cam.at("ImageResolution").get<std::string>());
            
            // Handle different camera formats
            if (cam.contains("NearPlane")) {
                camera.near_plane = parseVec4f(cam.at("NearPlane").get<std::string>());
            } else {
                // Default near plane if not specified
                camera.near_plane = {-1.0f, 1.0f, -1.0f, 1.0f};
            }
            
            camera.position = parseVec3f(cam.at("Position").get<std::string>());
            
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
                // Normalize the gaze vector
                normalize_on_place(camera.gaze);
                
            } else {
                // Default gaze direction
                camera.gaze = {0.0f, 0.0f, -1.0f};
            }
            
            camera.up = parseVec3f(cam.at("Up").get<std::string>());
            scene.cameras.push_back(camera);
        };
        if (sceneData.contains("Cameras") && sceneData["Cameras"].contains("Camera")) {
            processOneOrMany(sceneData.at("Cameras").at("Camera"), parseCamera);
        }

        // --- Parse Objects (Generalized) ---
        if (sceneData.contains("Objects")) {
            const auto& objects = sceneData.at("Objects");

            if (objects.contains("Sphere")) {
                auto parseSphere = [&](const json& obj) {
                    Sphere sphere;
                    sphere.id = std::stoi(obj.at("_id").get<std::string>());
                    sphere.material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    sphere.center_vertex_id = std::stoi(obj.at("Center").get<std::string>()) - 1;
                    sphere.radius = std::stof(obj.at("Radius").get<std::string>());
                    scene.sphere_data__.sphere_radius.push_back(sphere.radius);
                    scene.sphere_data__.sphere_center_vertex_id.push_back(sphere.center_vertex_id);
                    scene.sphere_data__.sphere_id.push_back(sphere.id);
                    scene.sphere_data__.sphere_mat_id.push_back(sphere.material_id);
                    
                    
                };
                processOneOrMany(objects.at("Sphere"), parseSphere);
            }
            
            if (objects.contains("Triangle")) {
                auto parseTriangle = [&](const json& obj) {
                    Triangle triangle;
                    triangle.id = std::stoi(obj.at("_id").get<std::string>());
                    triangle.material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    triangle.indices = parseIndices(obj.at("Indices").get<std::string>());
                    triangle.normal = normalize(cross(scene.vertex_data[triangle.indices.v2].position-scene.vertex_data[triangle.indices.v1].position,
                        scene.vertex_data[triangle.indices.v0].position-scene.vertex_data[triangle.indices.v1].position));
                    scene.triangle_data__.tri_norm_x.push_back(triangle.normal.x);    
                    scene.triangle_data__.tri_norm_y.push_back(triangle.normal.y);    
                    scene.triangle_data__.tri_norm_z.push_back(triangle.normal.z);
                    scene.triangle_data__.vo_ind.push_back(triangle.indices.v0);    
                    scene.triangle_data__.v1_ind.push_back(triangle.indices.v1);
                    scene.triangle_data__.v2_ind.push_back(triangle.indices.v2);
                    scene.triangle_data__.triangle_id.push_back(triangle.id);
                    scene.triangle_data__.triangle_material_id.push_back(triangle.material_id);
                };
                processOneOrMany(objects.at("Triangle"), parseTriangle);
            }
            
            // Add generalized loops for Mesh and Plane as well
            if (objects.contains("Mesh")) {
                auto parseMesh = [&](const json& obj) {
                    Mesh mesh;
                    mesh.id = std::stoi(obj.at("_id").get<std::string>());
                    mesh.material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    
                    std::string facesDataStr;
                    const auto& facesData = obj.at("Faces");
                    
                    // Handle both formats: direct string or object with _data field
                    if (facesData.is_string()) {
                        facesDataStr = facesData.get<std::string>();
                    } else if (facesData.is_object() && facesData.contains("_data")) {
                        facesDataStr = facesData.at("_data").get<std::string>();
                    }
                    
                    if (!facesDataStr.empty()) {
                        std::stringstream faceStream(facesDataStr);
                        int v0, v1, v2;
                        while (faceStream >> v0 >> v1 >> v2) {
                            // Convert from 1-based to 0-based indexing
                            int index_v0 = v0 - 1, index_v1 = v1 - 1, index_v2 = v2 - 1;
                            Vertex& v0_ = scene.vertex_data[index_v0];
                            Vertex& v1_ = scene.vertex_data[index_v1];
                            Vertex& v2_ = scene.vertex_data[index_v2];
                            
                            // Calculate face normal and area
                            Vec3f normal = cross(v1_.position - v0_.position, v2_.position - v0_.position);
                            mesh.faces.push_back({
                                {index_v0, index_v1, index_v2},
                                {normalize(normal)}
                            });

                            
                            
                            // Accumulate area-weighted normals for each vertex
                            v0_.normal += normal;
                            v1_.normal += normal;
                            v2_.normal += normal;
                        }
                    }
                    scene.meshes.push_back(mesh);
                };
                processOneOrMany(objects.at("Mesh"), parseMesh);
            }
            
            if (objects.contains("Plane")) {
                auto parsePlane = [&](const json& obj) {
                    Plane plane;
                    plane.id = std::stoi(obj.at("_id").get<std::string>());
                    plane.material_id = std::stoi(obj.at("Material").get<std::string>()) - 1;
                    plane.point_vertex_id = std::stoi(obj.at("Point").get<std::string>()) - 1;
                    plane.normal = parseVec3f(obj.at("Normal").get<std::string>());
                    scene.planes.push_back(plane);
                };
                processOneOrMany(objects.at("Plane"), parsePlane);
            }
        }
        for(auto& v:scene.vertex_data){
            normalize_on_place(v.normal);
        }
        return scene;
    }
};