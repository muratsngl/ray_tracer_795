#pragma once
#include <vector>
#include <random>
#include "types.hpp"

//use the mersenne twister ahmet hoca provided - now thread-safe
inline fl generate_rand_sample(){
    thread_local std::mt19937 generator(std::random_device{}());
    thread_local std::uniform_real_distribution<> distribution(0.0, 1.0);
    return distribution(generator);
}




//sample num samples amount of random numbers and then apply them to the pixels.
inline void jitter_sample(int n_sqrt,std::vector<RP8>& sample_packs,const RP8& reference_ray_pack,const Vec3f& delta_u, const Vec3f& delta_s,const Camera& cam,Vec3f camu,Vec3f camv){
    //for each sample
    fl step_size = 1.0f/n_sqrt;
    for(int i =0;i<n_sqrt;i++){
        for(int j =0;j<n_sqrt;j++){
            //for each pixel inside the rp8
            RP8 new_pack;
            for(int k = 0;k<8;k++){
                fl random_u = generate_rand_sample()*step_size+i*step_size;
                fl random_v = generate_rand_sample()*step_size+j*step_size;
                //if there is dof
                Vec3f random_displacement = random_u*delta_u-random_v*delta_s;
                if(cam.aperture_size>0.f){
                    //scaling with 0.5f as we are going to randomly sample from the cam center hence 
                    fl random_aperture_u=(generate_rand_sample()*2.f-1.0f)*0.5f;
                    fl random_aperture_v=(generate_rand_sample()*2.f-1.0f)*0.5f;
                    Vec3f a = cam.position+cam.aperture_size*random_aperture_v*camu+cam.aperture_size*random_aperture_u*camv;
                    Vec3f s_prim{reference_ray_pack.d_x[k] + random_displacement.x,reference_ray_pack.d_y[k] + random_displacement.y,reference_ray_pack.d_z[k] + random_displacement.z};
                    Vec3f dir = s_prim;
                    fl t = cam.focus_distance/dot(dir,cam.gaze);
                    Vec3f p = cam.position+t*dir;
                    Vec3f d = p-a;
                    new_pack.o_x[k] = a.x;
                    new_pack.o_y[k] = a.y;
                    new_pack.o_z[k] = a.z;
                    new_pack.d_x[k] = d.x;
                    new_pack.d_y[k] = d.y;
                    new_pack.d_z[k] = d.z;
                }
                //if there is no dof
                else{
                new_pack.o_x[k] = reference_ray_pack.o_x[k];
                new_pack.o_y[k] = reference_ray_pack.o_y[k];
                new_pack.o_z[k] = reference_ray_pack.o_z[k];
                new_pack.d_x[k] = reference_ray_pack.d_x[k] + random_displacement.x;
                new_pack.d_y[k] = reference_ray_pack.d_y[k] + random_displacement.y;
                new_pack.d_z[k] = reference_ray_pack.d_z[k] + random_displacement.z;
                }
                 
            }
            //std::cout<<"jittered_yah"<<std::endl;
            sample_packs.push_back(new_pack);
        }
    
    }

}
//later in a mysterious fashion, this will turn into a beautiful gaussian filter.
inline ColorBlockFl box_filter(ColorBlockFl cb_pack,int n){
    for(int i = 0;i<24;i++){
        cb_pack.rgb[i] = cb_pack.rgb[i]/n;
    }
    return cb_pack;
}