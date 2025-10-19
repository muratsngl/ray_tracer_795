



#pragma once
#ifndef RTFUNC
#define RTFUNC


#include <xsimd/xsimd.hpp>
#include "utils.hpp"
#include "types.hpp"

namespace xs = xsimd;

using batch_type = xs::batch<float,xs::avx2>;

void inline intersect_planes(const RP8& ray_pack, const PlaneData& planes,const VertexData& vertices){

}

void inline intersect_spheres(const RP8& ray_pack, const SphereData& spheres,const VertexData& vertices,const int Sphere_ID){
    auto dir_x = xs::load(&ray_pack.d_x[0]);
    auto dir_y = xs::load(&ray_pack.d_y[0]);
    auto dir_z = xs::load(&ray_pack.d_z[0]);

    auto orig_x = xs::load(&ray_pack.o_x);
    auto orig_y = xs::load(&ray_pack.o_y);
    auto orig_z = xs::load(&ray_pack.o_z);
    
    auto sphere_center_x = xs::broadcast(vertices.v_pos_x[Sphere_ID]);
    auto sphere_center_y = xs::broadcast(vertices.v_pos_y[Sphere_ID]);
    auto sphere_center_z = xs::broadcast(vertices.v_pos_z[Sphere_ID]);
    
    


    for(int i = 0;i<spheres.sphere_id.size();i++){

    }
}

void inline intersect_triangles(const RP8& ray_pack, const TriangleData& triangles,const VertexData& vertices){
    
}

void return_closest_hit(const RP8& ray_pack, const Scene& scene,const VertexData& vertices);







#endif