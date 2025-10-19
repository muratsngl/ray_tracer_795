#include <iostream>
#include "parser.cpp"
#include "rt_functions.hpp"

int main(int argc, char* argv[]) {
    Scene scene;
    // PARSE SCENE POPULATE STRUCTS
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <scene_file.json>" << std::endl;
        return 1;
    }
    
    try {
        scene = Parser::parseScene(argv[1]);
    } catch (const nlohmann::json::exception& e) {
        std::cerr << "JSON Error: " << e.what() << std::endl;
        std::cerr << "Error ID: " << e.id << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    //MAIN LOOP
    //INITIALIZE CAM AND RAYS
    for(const auto& cam:scene.cameras){
        Vec3f w = normalize(-cam.gaze);
        Vec3f u = normalize(cross(cam.up,w));
        Vec3f v = normalize(cross(w,u));
        
        //image plane mid and top left as ray references
        Vec3f plane_mid = (cam.near_distance*-w) + cam.position;
        Vec3f plane_top_left = plane_mid + cam.near_plane.x*u + cam.near_plane.w*v;

        Vec3f delta_su = (cam.near_plane.y-cam.near_plane.x)/cam.image_resolution.x*u;
        Vec3f delta_sv = (cam.near_plane.w-cam.near_plane.z)/cam.image_resolution.y*v;

        Vec3f dir = plane_top_left + (0.5f * delta_su) - (0.5f * delta_sv) - cam.position;
        
        for(int j = 0; j<cam.image_resolution.y;j+=2){
            for(int i = 0;i<cam.image_resolution.x;i+=4){
               Vec3f d_base = dir + (float)i * delta_su - (float)j * delta_sv;
                        //j
                        Vec3f d0 = d_base;//i
                        Vec3f d1 = d_base + delta_su;//i+1
                        Vec3f d2 = d_base + 2.0f * delta_su;//i+2
                        Vec3f d3 = d_base + 3.0f * delta_su;//i+3
                        //j+1
                        Vec3f d4 = d0 - delta_sv;
                        Vec3f d5 = d1 - delta_sv;
                        Vec3f d6 = d2 - delta_sv;
                        Vec3f d7 = d3 - delta_sv;

                        
                        //INTERSECTION && COLORING
                        RP8 ray_package{
                        // Origins (all 8 rays start at the camera position)
                            .o_x{cam.position.x, cam.position.x, cam.position.x, cam.position.x,
                                cam.position.x, cam.position.x, cam.position.x, cam.position.x},
                            .o_y{cam.position.y, cam.position.y, cam.position.y, cam.position.y,
                                cam.position.y, cam.position.y, cam.position.y, cam.position.y},
                            .o_z{cam.position.z, cam.position.z, cam.position.z, cam.position.z,
                                cam.position.z, cam.position.z, cam.position.z, cam.position.z},
                            .d_x{d0.x, d1.x, d2.x, d3.x, d4.x, d5.x, d6.x, d7.x},
                            .d_y{d0.y, d1.y, d2.y, d3.y, d4.y, d5.y, d6.y, d7.y},
                            .d_z{d0.z, d1.z, d2.z, d3.z, d4.z, d5.z, d6.z, d7.z}
                        };
                                
                        
                
                };
            }
        }


        
    


 

    //OUTPUT IMAGE
    


    return 0;
}

// CORE PHILOSOPHY: USE AS MANY CPU PARALELIZATION AS POSSIBLE. ESPECIALLY FOR BULK OPERATIONS LIKE INTERSECTION AND COLORING.