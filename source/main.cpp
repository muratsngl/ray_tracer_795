#include <iostream>
#include <chrono>
#include "parser.cpp"
#include "rt_functions.hpp"
#include "shader.hpp"
#include <omp.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"





int main(int argc, char* argv[]) {
    Scene scene;
    // PARSE SCENE POPULATE STRUCTS
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <scene_file.json>" << std::endl;
        return 1;
    }
    
    // Measure parsing time
    auto parse_start = std::chrono::high_resolution_clock::now();
    
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
    
    auto parse_end = std::chrono::high_resolution_clock::now();
    auto parse_duration = std::chrono::duration_cast<std::chrono::milliseconds>(parse_end - parse_start);
    std::cout << "Parsing time: " << parse_duration.count() << " ms" << std::endl;


    // Measure rendering time
    auto render_start = std::chrono::high_resolution_clock::now();

    //MAIN LOOP
    //INITIALIZE CAM AND RAYS
    for(const auto& cam:scene.cameras){
        
        // Use vector instead of VLA to avoid stack overflow with large images
        std::vector<unsigned char> image(cam.image_resolution.x * cam.image_resolution.y * 3);
        Vec3f w = normalize(-cam.gaze);
        Vec3f u = normalize(cross(cam.up,w));
        Vec3f v = normalize(cross(w,u));
        std::string image_name = cam.image_name;
        //image plane mid and top left as ray references
        Vec3f plane_mid = (cam.near_distance*-w) + cam.position;
        Vec3f plane_top_left = plane_mid + cam.near_plane.x*u + cam.near_plane.w*v;

        Vec3f delta_su = (cam.near_plane.y-cam.near_plane.x)/cam.image_resolution.x*u;
        Vec3f delta_sv = (cam.near_plane.w-cam.near_plane.z)/cam.image_resolution.y*v;

        Vec3f dir = plane_top_left + (0.5f * delta_su) - (0.5f * delta_sv) - cam.position;
        #pragma omp parallel for schedule(dynamic)
        for(int j = 0; j<cam.image_resolution.y;j+=2){
            for(int i = 0;i<cam.image_resolution.x;i+=4){
                ColorBlock cb;
                Vec3f d_base = dir + (fl)i * delta_su - (fl)j * delta_sv;
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

                        int current_mem_ptr = 3 * (j * cam.image_resolution.x + i);
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
                                
                        return_closest_hit(ray_package,scene,scene.vertex_data__);
                        
                        shade(ray_package,scene,cb);
                        
                        
                       
                        //implement write caching here(after bvh implementation) 
                        //shade_trad(ray_package,scene,cb);
                        memcpy(&image.data()[current_mem_ptr],cb.rgb,12*sizeof(unsigned char));
                        current_mem_ptr += cam.image_resolution.x*3;
                        memcpy(&image.data()[current_mem_ptr],&cb.rgb[12],12*sizeof(unsigned char));
                       
                       
                        
                        
                        
                        
                };
            
        }
            stbi_write_png(image_name.c_str(),cam.image_resolution.x,cam.image_resolution.y,3,image.data(),cam.image_resolution.x*3);
        }

    auto render_end = std::chrono::high_resolution_clock::now();
    auto render_duration = std::chrono::duration_cast<std::chrono::milliseconds>(render_end - render_start);
    std::cout << "Rendering time: " << render_duration.count() << " ms" << std::endl;
    std::cout << "Total time: " << (parse_duration.count() + render_duration.count()) << " ms" << std::endl;

        
    


 



    return 0;
}
