#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "json.hpp"
#include "types.hpp"

using json = nlohmann::json;

// Include the parser
#include "source/parser.cpp"

bool compare_float(float a, float b, float epsilon = 1e-5f) {
    return std::abs(a - b) < epsilon;
}

void test_vertex_data(const Scene& scene, const json& data) {
    std::cout << "\n=== Testing Vertex Data ===" << std::endl;
    
    const auto& sceneData = data["Scene"];
    if (!sceneData.contains("VertexData")) {
        std::cout << "✓ No vertex data in file, skipping" << std::endl;
        return;
    }
    
    const auto& vertexData = sceneData.at("VertexData");
    std::string vertexDataStr;
    
    if (vertexData.is_string()) {
        vertexDataStr = vertexData.get<std::string>();
    } else if (vertexData.is_object() && vertexData.contains("_data")) {
        vertexDataStr = vertexData.at("_data").get<std::string>();
    }
    
    // Parse expected vertices from file
    std::vector<float> expected_x, expected_y, expected_z;
    std::stringstream vertexStream(vertexDataStr);
    float x, y, z;
    while (vertexStream >> x >> y >> z) {
        expected_x.push_back(x);
        expected_y.push_back(y);
        expected_z.push_back(z);
    }
    
    // Compare with parsed data
    if (scene.vertex_data__.v_pos_x.size() != expected_x.size()) {
        std::cout << "✗ FAIL: Vertex count mismatch! Expected: " << expected_x.size() 
                  << ", Got: " << scene.vertex_data__.v_pos_x.size() << std::endl;
        return;
    }
    
    std::cout << "✓ Vertex count matches: " << expected_x.size() << std::endl;
    
    // Check first few vertices
    bool all_match = true;
    for (size_t i = 0; i < std::min(size_t(5), expected_x.size()); ++i) {
        if (!compare_float(scene.vertex_data__.v_pos_x[i], expected_x[i]) ||
            !compare_float(scene.vertex_data__.v_pos_y[i], expected_y[i]) ||
            !compare_float(scene.vertex_data__.v_pos_z[i], expected_z[i])) {
            std::cout << "✗ Vertex " << i << " mismatch!" << std::endl;
            std::cout << "  Expected: (" << expected_x[i] << ", " << expected_y[i] << ", " << expected_z[i] << ")" << std::endl;
            std::cout << "  Got:      (" << scene.vertex_data__.v_pos_x[i] << ", " 
                      << scene.vertex_data__.v_pos_y[i] << ", " << scene.vertex_data__.v_pos_z[i] << ")" << std::endl;
            all_match = false;
        }
    }
    
    if (all_match) {
        std::cout << "✓ First 5 vertices match correctly" << std::endl;
    }
}

void test_point_lights(const Scene& scene, const json& data) {
    std::cout << "\n=== Testing Point Lights ===" << std::endl;
    
    const auto& sceneData = data["Scene"];
    if (!sceneData.contains("Lights") || !sceneData["Lights"].contains("PointLight")) {
        std::cout << "✓ No point lights in file, skipping" << std::endl;
        return;
    }
    
    const auto& lights = sceneData["Lights"]["PointLight"];
    
    auto parseLights = [&](const json& light_array) {
        size_t expected_count = 0;
        if (light_array.is_array()) {
            expected_count = light_array.size();
        } else if (light_array.is_object()) {
            expected_count = 1;
        }
        
        if (scene.point_light_data__.pl_id.size() != expected_count) {
            std::cout << "✗ FAIL: Light count mismatch! Expected: " << expected_count 
                      << ", Got: " << scene.point_light_data__.pl_id.size() << std::endl;
            return;
        }
        
        std::cout << "✓ Point light count matches: " << expected_count << std::endl;
        
        // Check first light
        if (expected_count > 0) {
            const json& first_light = light_array.is_array() ? light_array[0] : light_array;
            
            float pos_x, pos_y, pos_z;
            std::stringstream pos_stream(first_light.at("Position").get<std::string>());
            pos_stream >> pos_x >> pos_y >> pos_z;
            
            float int_r, int_g, int_b;
            std::stringstream int_stream(first_light.at("Intensity").get<std::string>());
            int_stream >> int_r >> int_g >> int_b;
            
            bool match = compare_float(scene.point_light_data__.pl_pos_x[0], pos_x) &&
                         compare_float(scene.point_light_data__.pl_pos_y[0], pos_y) &&
                         compare_float(scene.point_light_data__.pl_pos_z[0], pos_z) &&
                         compare_float(scene.point_light_data__.pl_intensity_r[0], int_r) &&
                         compare_float(scene.point_light_data__.pl_intensity_g[0], int_g) &&
                         compare_float(scene.point_light_data__.pl_intensity_b[0], int_b);
            
            if (match) {
                std::cout << "✓ First point light data matches" << std::endl;
            } else {
                std::cout << "✗ First point light data mismatch!" << std::endl;
            }
        }
    };
    
    parseLights(lights);
}

void test_spheres(const Scene& scene, const json& data) {
    std::cout << "\n=== Testing Spheres ===" << std::endl;
    
    const auto& sceneData = data["Scene"];
    if (!sceneData.contains("Objects") || !sceneData["Objects"].contains("Sphere")) {
        std::cout << "✓ No spheres in file, skipping" << std::endl;
        return;
    }
    
    const auto& spheres = sceneData["Objects"]["Sphere"];
    
    size_t expected_count = 0;
    if (spheres.is_array()) {
        expected_count = spheres.size();
    } else if (spheres.is_object()) {
        expected_count = 1;
    }
    
    if (scene.sphere_data__.sphere_id.size() != expected_count) {
        std::cout << "✗ FAIL: Sphere count mismatch! Expected: " << expected_count 
                  << ", Got: " << scene.sphere_data__.sphere_id.size() << std::endl;
        return;
    }
    
    std::cout << "✓ Sphere count matches: " << expected_count << std::endl;
    
    // Check first sphere
    if (expected_count > 0) {
        const json& first_sphere = spheres.is_array() ? spheres[0] : spheres;
        
        int expected_id = std::stoi(first_sphere.at("_id").get<std::string>());
        float expected_radius = std::stof(first_sphere.at("Radius").get<std::string>());
        int expected_mat_id = std::stoi(first_sphere.at("Material").get<std::string>()) - 1;
        int expected_center = std::stoi(first_sphere.at("Center").get<std::string>()) - 1;
        
        bool match = (scene.sphere_data__.sphere_id[0] == expected_id) &&
                     compare_float(scene.sphere_data__.sphere_radius[0], expected_radius) &&
                     (scene.sphere_data__.sphere_mat_id[0] == expected_mat_id) &&
                     (scene.sphere_data__.sphere_center_vertex_id[0] == expected_center);
        
        if (match) {
            std::cout << "✓ First sphere data matches" << std::endl;
        } else {
            std::cout << "✗ First sphere data mismatch!" << std::endl;
            std::cout << "  Expected: id=" << expected_id << " radius=" << expected_radius 
                      << " mat=" << expected_mat_id << " center=" << expected_center << std::endl;
            std::cout << "  Got:      id=" << scene.sphere_data__.sphere_id[0] 
                      << " radius=" << scene.sphere_data__.sphere_radius[0]
                      << " mat=" << scene.sphere_data__.sphere_mat_id[0] 
                      << " center=" << scene.sphere_data__.sphere_center_vertex_id[0] << std::endl;
        }
    }
}

void test_triangles(const Scene& scene, const json& data) {
    std::cout << "\n=== Testing Triangles ===" << std::endl;
    
    const auto& sceneData = data["Scene"];
    size_t standalone_count = 0;
    
    if (!sceneData.contains("Objects") || !sceneData["Objects"].contains("Triangle")) {
        std::cout << "✓ No standalone triangles in file" << std::endl;
    } else {
        const auto& triangles = sceneData["Objects"]["Triangle"];
        standalone_count = triangles.is_array() ? triangles.size() : 1;
        std::cout << "  Standalone triangles in file: " << standalone_count << std::endl;
        
        // Check first standalone triangle indices
        const json& first_tri = triangles.is_array() ? triangles[0] : triangles;
        std::string indices_str = first_tri.at("Indices").get<std::string>();
        int v0_file, v1_file, v2_file;
        std::stringstream ind_stream(indices_str);
        ind_stream >> v0_file >> v1_file >> v2_file;
        
        // JSON uses 1-based, parser should convert to 0-based
        int expected_v0 = v0_file - 1;
        int expected_v1 = v1_file - 1;
        int expected_v2 = v2_file - 1;
        
        if (scene.triangle_data__.vo_ind[0] == expected_v0 &&
            scene.triangle_data__.v1_ind[0] == expected_v1 &&
            scene.triangle_data__.v2_ind[0] == expected_v2) {
            std::cout << "✓ Triangle index conversion (1-based → 0-based) correct!" << std::endl;
            std::cout << "  File: (" << v0_file << ", " << v1_file << ", " << v2_file << ") → ";
            std::cout << "Parsed: (" << expected_v0 << ", " << expected_v1 << ", " << expected_v2 << ")" << std::endl;
        } else {
            std::cout << "✗ Triangle index conversion FAILED!" << std::endl;
            std::cout << "  Expected: (" << expected_v0 << ", " << expected_v1 << ", " << expected_v2 << ")" << std::endl;
            std::cout << "  Got:      (" << scene.triangle_data__.vo_ind[0] << ", " 
                      << scene.triangle_data__.v1_ind[0] << ", " << scene.triangle_data__.v2_ind[0] << ")" << std::endl;
        }
    }
    
    // Count mesh faces
    size_t mesh_face_count = 0;
    if (sceneData.contains("Objects") && sceneData["Objects"].contains("Mesh")) {
        const auto& meshes = sceneData["Objects"]["Mesh"];
        
        auto countFaces = [&](const json& mesh) {
            const auto& facesData = mesh.at("Faces");
            std::string facesDataStr;
            
            if (facesData.is_string()) {
                facesDataStr = facesData.get<std::string>();
            } else if (facesData.is_object() && facesData.contains("_data")) {
                facesDataStr = facesData.at("_data").get<std::string>();
            }
            
            std::stringstream faceStream(facesDataStr);
            int v0, v1, v2;
            while (faceStream >> v0 >> v1 >> v2) {
                mesh_face_count++;
            }
        };
        
        if (meshes.is_array()) {
            for (const auto& mesh : meshes) {
                countFaces(mesh);
            }
        } else {
            countFaces(meshes);
        }
        
        std::cout << "  Mesh faces in file: " << mesh_face_count << std::endl;
        
        // Check first mesh face indexing
        const json& first_mesh = meshes.is_array() ? meshes[0] : meshes;
        const auto& facesData = first_mesh.at("Faces");
        std::string facesDataStr;
        
        if (facesData.is_string()) {
            facesDataStr = facesData.get<std::string>();
        } else if (facesData.is_object() && facesData.contains("_data")) {
            facesDataStr = facesData.at("_data").get<std::string>();
        }
        
        std::stringstream faceStream(facesDataStr);
        int v0_file, v1_file, v2_file;
        faceStream >> v0_file >> v1_file >> v2_file;
        
        // Mesh face should be at index standalone_count (after standalone triangles)
        size_t mesh_tri_idx = standalone_count;
        int expected_v0 = v0_file - 1;
        int expected_v1 = v1_file - 1;
        int expected_v2 = v2_file - 1;
        
        if (mesh_tri_idx < scene.triangle_data__.vo_ind.size() &&
            scene.triangle_data__.vo_ind[mesh_tri_idx] == expected_v0 &&
            scene.triangle_data__.v1_ind[mesh_tri_idx] == expected_v1 &&
            scene.triangle_data__.v2_ind[mesh_tri_idx] == expected_v2) {
            std::cout << "✓ Mesh face index conversion (1-based → 0-based) correct!" << std::endl;
            std::cout << "  File: (" << v0_file << ", " << v1_file << ", " << v2_file << ") → ";
            std::cout << "Parsed: (" << expected_v0 << ", " << expected_v1 << ", " << expected_v2 << ")" << std::endl;
        } else if (mesh_tri_idx < scene.triangle_data__.vo_ind.size()) {
            std::cout << "✗ Mesh face index conversion FAILED!" << std::endl;
            std::cout << "  Expected: (" << expected_v0 << ", " << expected_v1 << ", " << expected_v2 << ")" << std::endl;
            std::cout << "  Got:      (" << scene.triangle_data__.vo_ind[mesh_tri_idx] << ", " 
                      << scene.triangle_data__.v1_ind[mesh_tri_idx] << ", " 
                      << scene.triangle_data__.v2_ind[mesh_tri_idx] << ")" << std::endl;
        }
    }
    
    std::cout << "  Total triangles in parsed data: " << scene.triangle_data__.vo_ind.size() << std::endl;
    
    if (scene.triangle_data__.vo_ind.size() > 0) {
        std::cout << "✓ Triangle data populated" << std::endl;
    }
}

void test_sphere_indexing(const Scene& scene, const json& data) {
    std::cout << "\n=== Testing Sphere Center Indexing ===" << std::endl;
    
    const auto& sceneData = data["Scene"];
    if (!sceneData.contains("Objects") || !sceneData["Objects"].contains("Sphere")) {
        std::cout << "✓ No spheres to test" << std::endl;
        return;
    }
    
    const auto& spheres = sceneData["Objects"]["Sphere"];
    const json& first_sphere = spheres.is_array() ? spheres[0] : spheres;
    
    int center_file = std::stoi(first_sphere.at("Center").get<std::string>());
    int expected_center = center_file - 1;  // Convert to 0-based
    
    if (scene.sphere_data__.sphere_center_vertex_id[0] == expected_center) {
        std::cout << "✓ Sphere center index conversion (1-based → 0-based) correct!" << std::endl;
        std::cout << "  File: " << center_file << " → Parsed: " << expected_center << std::endl;
    } else {
        std::cout << "✗ Sphere center index conversion FAILED!" << std::endl;
        std::cout << "  Expected: " << expected_center << std::endl;
        std::cout << "  Got:      " << scene.sphere_data__.sphere_center_vertex_id[0] << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::string scene_file = "inputs/simple.json";
    
    if (argc > 1) {
        scene_file = argv[1];
    }
    
    std::cout << "\n╔════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║           DATA INTEGRITY TEST                             ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════╝" << std::endl;
    std::cout << "\nTesting: " << scene_file << std::endl;
    
    try {
        // Load JSON directly
        std::ifstream file(scene_file);
        if (!file.is_open()) {
            std::cerr << "✗ Cannot open file: " << scene_file << std::endl;
            return 1;
        }
        json data;
        file >> data;
        
        // Parse using our parser
        Scene scene = Parser::parseScene(scene_file);
        
        // Run tests
        test_vertex_data(scene, data);
        test_point_lights(scene, data);
        test_spheres(scene, data);
        test_sphere_indexing(scene, data);
        test_triangles(scene, data);
        
        std::cout << "\n╔════════════════════════════════════════════════════════════╗" << std::endl;
        std::cout << "║              TESTS COMPLETED                              ║" << std::endl;
        std::cout << "╚════════════════════════════════════════════════════════════╝\n" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ ERROR: " << e.what() << "\n" << std::endl;
        return 1;
    }
}
