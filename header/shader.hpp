#pragma once
#include "types.hpp"
#include "utils.hpp"
#include "rt_functions.hpp"

// Convert ColorBlockFl (fl) to ColorBlock (unsigned char) with clamping
void inline convert_to_unsigned_char(const ColorBlockFl& color_block_int, ColorBlock& color_block){
    for(int i = 0; i < 24; i++){
        // Scale to [0, 255] range, then clamp to [0, 255]
        fl value = color_block_int.rgb[i];
        if(value < 0.0f) value = 0.0f;
        if(value > 255.0f) value = 255.0f;
        color_block.rgb[i] = static_cast<unsigned char>(value);
    }
}


void inline no_shades(const RP8& ray_pack, const Scene& scene, ColorBlock& color_block){
    // ColorBlock layout: flat array of 24 ints (8 rays × 3 RGB channels)
    // Index = ray_index * 3 + channel (0=R, 1=G, 2=B)
    
    // Ray 0 - RGB (indices 0, 1, 2)
    if(ray_pack.mat_id[0] >= 0) {
        color_block.rgb[0] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[0]].diffuse_reflectance.x * 255);
        color_block.rgb[1] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[0]].diffuse_reflectance.y * 255);
        color_block.rgb[2] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[0]].diffuse_reflectance.z * 255);
    } else {
        color_block.rgb[0] = static_cast<unsigned char>(scene.background_color.x * 255);
        color_block.rgb[1] = static_cast<unsigned char>(scene.background_color.y * 255);
        color_block.rgb[2] = static_cast<unsigned char>(scene.background_color.z * 255);
    }
    
    // Ray 1 - RGB (indices 3, 4, 5)
    if(ray_pack.mat_id[1] >= 0) {
        color_block.rgb[3] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[1]].diffuse_reflectance.x * 255);
        color_block.rgb[4] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[1]].diffuse_reflectance.y * 255);
        color_block.rgb[5] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[1]].diffuse_reflectance.z * 255);
    } else {
        color_block.rgb[3] = static_cast<unsigned char>(scene.background_color.x * 255);
        color_block.rgb[4] = static_cast<unsigned char>(scene.background_color.y * 255);
        color_block.rgb[5] = static_cast<unsigned char>(scene.background_color.z * 255);
    }
    
    // Ray 2 - RGB (indices 6, 7, 8)
    if(ray_pack.mat_id[2] >= 0) {
        color_block.rgb[6] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[2]].diffuse_reflectance.x * 255);
        color_block.rgb[7] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[2]].diffuse_reflectance.y * 255);
        color_block.rgb[8] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[2]].diffuse_reflectance.z * 255);
    } else {
        color_block.rgb[6] = static_cast<unsigned char>(scene.background_color.x * 255);
        color_block.rgb[7] = static_cast<unsigned char>(scene.background_color.y * 255);
        color_block.rgb[8] = static_cast<unsigned char>(scene.background_color.z * 255);
    }
    
    // Ray 3 - RGB (indices 9, 10, 11)
    if(ray_pack.mat_id[3] >= 0) {
        color_block.rgb[9] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[3]].diffuse_reflectance.x * 255);
        color_block.rgb[10] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[3]].diffuse_reflectance.y * 255);
        color_block.rgb[11] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[3]].diffuse_reflectance.z * 255);
    } else {
        color_block.rgb[9] = static_cast<unsigned char>(scene.background_color.x * 255);
        color_block.rgb[10] = static_cast<unsigned char>(scene.background_color.y * 255);
        color_block.rgb[11] = static_cast<unsigned char>(scene.background_color.z * 255);
    }
    
    // Ray 4 - RGB (indices 12, 13, 14)
    if(ray_pack.mat_id[4] >= 0) {
        color_block.rgb[12] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[4]].diffuse_reflectance.x * 255);
        color_block.rgb[13] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[4]].diffuse_reflectance.y * 255);
        color_block.rgb[14] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[4]].diffuse_reflectance.z * 255);
    } else {
        color_block.rgb[12] = static_cast<unsigned char>(scene.background_color.x * 255);
        color_block.rgb[13] = static_cast<unsigned char>(scene.background_color.y * 255);
        color_block.rgb[14] = static_cast<unsigned char>(scene.background_color.z * 255);
    }
    
    // Ray 5 - RGB (indices 15, 16, 17)
    if(ray_pack.mat_id[5] >= 0) {
        color_block.rgb[15] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[5]].diffuse_reflectance.x * 255);
        color_block.rgb[16] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[5]].diffuse_reflectance.y * 255);
        color_block.rgb[17] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[5]].diffuse_reflectance.z * 255);
    } else {
        color_block.rgb[15] = static_cast<unsigned char>(scene.background_color.x * 255);
        color_block.rgb[16] = static_cast<unsigned char>(scene.background_color.y * 255);
        color_block.rgb[17] = static_cast<unsigned char>(scene.background_color.z * 255);
    }
    
    // Ray 6 - RGB (indices 18, 19, 20)
    if(ray_pack.mat_id[6] >= 0) {
        color_block.rgb[18] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[6]].diffuse_reflectance.x * 255);
        color_block.rgb[19] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[6]].diffuse_reflectance.y * 255);
        color_block.rgb[20] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[6]].diffuse_reflectance.z * 255);
    } else {
        color_block.rgb[18] = static_cast<unsigned char>(scene.background_color.x * 255);
        color_block.rgb[19] = static_cast<unsigned char>(scene.background_color.y * 255);
        color_block.rgb[20] = static_cast<unsigned char>(scene.background_color.z * 255);
    }
    
    // Ray 7 - RGB (indices 21, 22, 23)
    if(ray_pack.mat_id[7] >= 0) {
        color_block.rgb[21] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[7]].diffuse_reflectance.x * 255);
        color_block.rgb[22] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[7]].diffuse_reflectance.y * 255);
        color_block.rgb[23] = static_cast<unsigned char>(scene.materials[ray_pack.mat_id[7]].diffuse_reflectance.z * 255);
    } else {
        color_block.rgb[21] = static_cast<unsigned char>(scene.background_color.x * 255);
        color_block.rgb[22] = static_cast<unsigned char>(scene.background_color.y * 255);
        color_block.rgb[23] = static_cast<unsigned char>(scene.background_color.z * 255);
    }
}

void inline no_shades_white(const RP8& ray_pack, const Scene& scene, ColorBlock& color_block){
    // ColorBlock layout: flat array of 24 ints (8 rays × 3 RGB channels)
    // Index = ray_index * 3 + channel (0=R, 1=G, 2=B)
    // White for hits, black for misses
    
    // Ray 0 - RGB (indices 0, 1, 2)
    if(ray_pack.mat_id[0] >= 0) {
        color_block.rgb[0] = 255;
        color_block.rgb[1] = 255;
        color_block.rgb[2] = 255;
    } else {
        color_block.rgb[0] = 0;
        color_block.rgb[1] = 0;
        color_block.rgb[2] = 0;
    }
    
    // Ray 1 - RGB (indices 3, 4, 5)
    if(ray_pack.mat_id[1] >= 0) {
        color_block.rgb[3] = 255;
        color_block.rgb[4] = 255;
        color_block.rgb[5] = 255;
    } else {
        color_block.rgb[3] = 0;
        color_block.rgb[4] = 0;
        color_block.rgb[5] = 0;
    }
    
    // Ray 2 - RGB (indices 6, 7, 8)
    if(ray_pack.mat_id[2] >= 0) {
        color_block.rgb[6] = 255;
        color_block.rgb[7] = 255;
        color_block.rgb[8] = 255;
    } else {
        color_block.rgb[6] = 0;
        color_block.rgb[7] = 0;
        color_block.rgb[8] = 0;
    }
    
    // Ray 3 - RGB (indices 9, 10, 11)
    if(ray_pack.mat_id[3] >= 0) {
        color_block.rgb[9] = 255;
        color_block.rgb[10] = 255;
        color_block.rgb[11] = 255;
    } else {
        color_block.rgb[9] = 0;
        color_block.rgb[10] = 0;
        color_block.rgb[11] = 0;
    }
    
    // Ray 4 - RGB (indices 12, 13, 14)
    if(ray_pack.mat_id[4] >= 0) {
        color_block.rgb[12] = 255;
        color_block.rgb[13] = 255;
        color_block.rgb[14] = 255;
    } else {
        color_block.rgb[12] = 0;
        color_block.rgb[13] = 0;
        color_block.rgb[14] = 0;
    }
    
    // Ray 5 - RGB (indices 15, 16, 17)
    if(ray_pack.mat_id[5] >= 0) {
        color_block.rgb[15] = 255;
        color_block.rgb[16] = 255;
        color_block.rgb[17] = 255;
    } else {
        color_block.rgb[15] = 0;
        color_block.rgb[16] = 0;
        color_block.rgb[17] = 0;
    }
    
    // Ray 6 - RGB (indices 18, 19, 20)
    if(ray_pack.mat_id[6] >= 0) {
        color_block.rgb[18] = 255;
        color_block.rgb[19] = 255;
        color_block.rgb[20] = 255;
    } else {
        color_block.rgb[18] = 0;
        color_block.rgb[19] = 0;
        color_block.rgb[20] = 0;
    }
    
    // Ray 7 - RGB (indices 21, 22, 23)
    if(ray_pack.mat_id[7] >= 0) {
        color_block.rgb[21] = 255;
        color_block.rgb[22] = 255;
        color_block.rgb[23] = 255;
    } else {
        color_block.rgb[21] = 0;
        color_block.rgb[22] = 0;
        color_block.rgb[23] = 0;
    }
}

void inline diffuse(const RP8& ray_pack, const Scene& scene, ColorBlockFl& color_block, xsimd::batch_bool<fl>& in_light, int light_index,
                    const f_batch& hit_x, const f_batch& hit_y, const f_batch& hit_z,
                    const f_batch& hit_norm_x, const f_batch& hit_norm_y, const f_batch& hit_norm_z,
                    const f_batch& light_pos_x, const f_batch& light_pos_y, const f_batch& light_pos_z,
                    const f_batch& light_dir_x, const f_batch& light_dir_y, const f_batch& light_dir_z,
                    const f_batch& ray_dir_x, const f_batch& ray_dir_y, const f_batch& ray_dir_z,
                    const f_batch& distance){
    
    // Step 1: Calculate N·L (diffuse intensity) using SIMD
    f_batch NdotL = dot_simd(hit_norm_x, hit_norm_y, hit_norm_z, 
                             light_dir_x, light_dir_y, light_dir_z);
    
    // Step 2: Clamp to [0, infinity] - Lambert's law
    f_batch zero_batch = xs::broadcast(0.0f);
    f_batch diffuse_intensity = xs::max(zero_batch, NdotL);
    
    // Step 3: Load material diffuse colors for all 8 rays (unfolded loop)
    alignas(32) fl mat_r[8], mat_g[8], mat_b[8];
    
    // Ray 0
    mat_r[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].diffuse_reflectance.x : 0.0f;
    mat_g[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].diffuse_reflectance.y : 0.0f;
    mat_b[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].diffuse_reflectance.z : 0.0f;
    
    // Ray 1
    mat_r[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].diffuse_reflectance.x : 0.0f;
    mat_g[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].diffuse_reflectance.y : 0.0f;
    mat_b[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].diffuse_reflectance.z : 0.0f;
    
    // Ray 2
    mat_r[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].diffuse_reflectance.x : 0.0f;
    mat_g[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].diffuse_reflectance.y : 0.0f;
    mat_b[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].diffuse_reflectance.z : 0.0f;
    
    // Ray 3
    mat_r[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].diffuse_reflectance.x : 0.0f;
    mat_g[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].diffuse_reflectance.y : 0.0f;
    mat_b[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].diffuse_reflectance.z : 0.0f;
    
    // Ray 4
    mat_r[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].diffuse_reflectance.x : 0.0f;
    mat_g[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].diffuse_reflectance.y : 0.0f;
    mat_b[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].diffuse_reflectance.z : 0.0f;
    
    // Ray 5
    mat_r[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].diffuse_reflectance.x : 0.0f;
    mat_g[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].diffuse_reflectance.y : 0.0f;
    mat_b[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].diffuse_reflectance.z : 0.0f;
    
    // Ray 6
    mat_r[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].diffuse_reflectance.x : 0.0f;
    mat_g[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].diffuse_reflectance.y : 0.0f;
    mat_b[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].diffuse_reflectance.z : 0.0f;
    
    // Ray 7
    mat_r[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].diffuse_reflectance.x : 0.0f;
    mat_g[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].diffuse_reflectance.y : 0.0f;
    mat_b[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].diffuse_reflectance.z : 0.0f;
    
    // Load into SIMD batches
    f_batch material_r = xs::load_aligned(mat_r);
    f_batch material_g = xs::load_aligned(mat_g);
    f_batch material_b = xs::load_aligned(mat_b);
    
    // Step 4: Load light intensity as broadcast (same for all 8 rays)
    f_batch light_r = xs::broadcast(scene.point_light_data__.pl_intensity_r[light_index]);
    f_batch light_g = xs::broadcast(scene.point_light_data__.pl_intensity_g[light_index]);
    f_batch light_b = xs::broadcast(scene.point_light_data__.pl_intensity_b[light_index]);
    
    // Step 5: Apply distance attenuation (1 / distance^2)
    f_batch distance_squared = distance * distance;
    f_batch attenuation = xs::broadcast(1.0f) / distance_squared;
    
    // Step 6: Calculate diffuse color with attenuation (Lambert: I = kd * light * max(0, N·L) / d^2)
    f_batch color_r = material_r * light_r * diffuse_intensity * attenuation;
    f_batch color_g = material_g * light_g * diffuse_intensity * attenuation;
    f_batch color_b = material_b * light_b * diffuse_intensity * attenuation;
    
    // Step 7: Apply shadow mask (only add color where not in shadow)
    color_r = xs::select(in_light, color_r, zero_batch);
    color_g = xs::select(in_light, color_g, zero_batch);
    color_b = xs::select(in_light, color_b, zero_batch);
    
    // Step 8: Store results back to color_block as FL (no conversion yet)
    alignas(32) fl final_r[8], final_g[8], final_b[8];
    xs::store_aligned(final_r, color_r);
    xs::store_aligned(final_g, color_g);
    xs::store_aligned(final_b, color_b);
    
    // Accumulate into ColorBlockFl (keep as fl, no *255)
    color_block.rgb[0]  += final_r[0];
    color_block.rgb[1]  += final_g[0];
    color_block.rgb[2]  += final_b[0];
    
    color_block.rgb[3]  += final_r[1];
    color_block.rgb[4]  += final_g[1];
    color_block.rgb[5]  += final_b[1];
    
    color_block.rgb[6]  += final_r[2];
    color_block.rgb[7]  += final_g[2];
    color_block.rgb[8]  += final_b[2];
    
    color_block.rgb[9]  += final_r[3];
    color_block.rgb[10] += final_g[3];
    color_block.rgb[11] += final_b[3];
    
    color_block.rgb[12] += final_r[4];
    color_block.rgb[13] += final_g[4];
    color_block.rgb[14] += final_b[4];
    
    color_block.rgb[15] += final_r[5];
    color_block.rgb[16] += final_g[5];
    color_block.rgb[17] += final_b[5];
    
    color_block.rgb[18] += final_r[6];
    color_block.rgb[19] += final_g[6];
    color_block.rgb[20] += final_b[6];
    
    color_block.rgb[21] += final_r[7];
    color_block.rgb[22] += final_g[7];
    color_block.rgb[23] += final_b[7];
}

void inline ambient(const RP8& ray_pack, const Scene& scene, ColorBlockFl& color_block){
    
    // Step 1: Load material ambient colors for all 8 rays (unfolded loop)
    alignas(32) fl mat_r[8], mat_g[8], mat_b[8];
    
    // Ray 0
    mat_r[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].ambient_reflectance.x : 0.0f;
    mat_g[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].ambient_reflectance.y : 0.0f;
    mat_b[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].ambient_reflectance.z : 0.0f;
    
    // Ray 1
    mat_r[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].ambient_reflectance.x : 0.0f;
    mat_g[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].ambient_reflectance.y : 0.0f;
    mat_b[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].ambient_reflectance.z : 0.0f;
    
    // Ray 2
    mat_r[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].ambient_reflectance.x : 0.0f;
    mat_g[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].ambient_reflectance.y : 0.0f;
    mat_b[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].ambient_reflectance.z : 0.0f;
    
    // Ray 3
    mat_r[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].ambient_reflectance.x : 0.0f;
    mat_g[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].ambient_reflectance.y : 0.0f;
    mat_b[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].ambient_reflectance.z : 0.0f;
    
    // Ray 4
    mat_r[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].ambient_reflectance.x : 0.0f;
    mat_g[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].ambient_reflectance.y : 0.0f;
    mat_b[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].ambient_reflectance.z : 0.0f;
    
    // Ray 5
    mat_r[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].ambient_reflectance.x : 0.0f;
    mat_g[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].ambient_reflectance.y : 0.0f;
    mat_b[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].ambient_reflectance.z : 0.0f;
    
    // Ray 6
    mat_r[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].ambient_reflectance.x : 0.0f;
    mat_g[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].ambient_reflectance.y : 0.0f;
    mat_b[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].ambient_reflectance.z : 0.0f;
    
    // Ray 7
    mat_r[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].ambient_reflectance.x : 0.0f;
    mat_g[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].ambient_reflectance.y : 0.0f;
    mat_b[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].ambient_reflectance.z : 0.0f;
    
    // Load into SIMD batches
    f_batch material_r = xs::load_aligned(mat_r);
    f_batch material_g = xs::load_aligned(mat_g);
    f_batch material_b = xs::load_aligned(mat_b);
    
    // Step 2: Load ambient light intensity as broadcast (same for all 8 rays)
    f_batch ambient_r = xs::broadcast(scene.ambient_light.x);
    f_batch ambient_g = xs::broadcast(scene.ambient_light.y);
    f_batch ambient_b = xs::broadcast(scene.ambient_light.z);
    
    // Step 3: Calculate ambient color (I = ka * ambient_light)
    f_batch color_r = material_r * ambient_r;
    f_batch color_g = material_g * ambient_g;
    f_batch color_b = material_b * ambient_b;
    
    // Step 4: Store results back to color_block as FL (no conversion yet)
    alignas(32) fl final_r[8], final_g[8], final_b[8];
    xs::store_aligned(final_r, color_r);
    xs::store_aligned(final_g, color_g);
    xs::store_aligned(final_b, color_b);
    
    // Accumulate into ColorBlockFl (keep as fl)
    color_block.rgb[0]  += final_r[0];
    color_block.rgb[1]  += final_g[0];
    color_block.rgb[2]  += final_b[0];
    
    color_block.rgb[3]  += final_r[1];
    color_block.rgb[4]  += final_g[1];
    color_block.rgb[5]  += final_b[1];
    
    color_block.rgb[6]  += final_r[2];
    color_block.rgb[7]  += final_g[2];
    color_block.rgb[8]  += final_b[2];
    
    color_block.rgb[9]  += final_r[3];
    color_block.rgb[10] += final_g[3];
    color_block.rgb[11] += final_b[3];
    
    color_block.rgb[12] += final_r[4];
    color_block.rgb[13] += final_g[4];
    color_block.rgb[14] += final_b[4];
    
    color_block.rgb[15] += final_r[5];
    color_block.rgb[16] += final_g[5];
    color_block.rgb[17] += final_b[5];
    
    color_block.rgb[18] += final_r[6];
    color_block.rgb[19] += final_g[6];
    color_block.rgb[20] += final_b[6];
    
    color_block.rgb[21] += final_r[7];
    color_block.rgb[22] += final_g[7];
    color_block.rgb[23] += final_b[7];
}

void inline specular(const RP8& ray_pack, const Scene& scene, ColorBlockFl& color_block, xsimd::batch_bool<fl>& in_light, int light_index,
                     const f_batch& hit_x, const f_batch& hit_y, const f_batch& hit_z,
                     const f_batch& hit_norm_x, const f_batch& hit_norm_y, const f_batch& hit_norm_z,
                     const f_batch& light_pos_x, const f_batch& light_pos_y, const f_batch& light_pos_z,
                     const f_batch& light_dir_x, const f_batch& light_dir_y, const f_batch& light_dir_z,
                     const f_batch& ray_dir_x, const f_batch& ray_dir_y, const f_batch& ray_dir_z,
                     const f_batch& distance){
    
    // --- BLINN-PHONG IMPLEMENTATION ---

    // Step 1: Get View vector V (from hit point to camera)
    f_batch view_x = -ray_dir_x;
    f_batch view_y = -ray_dir_y;
    f_batch view_z = -ray_dir_z;
    
    // Step 2: Calculate Halfway vector H = L + V
    // (L is light_dir_x, etc.)
    f_batch half_x = light_dir_x + view_x;
    f_batch half_y = light_dir_y + view_y;
    f_batch half_z = light_dir_z + view_z;
    // Normalize the halfway vector
    normalize_simd_overwrite(half_x, half_y, half_z);
    
    // Step 3: Calculate N·H (dot product of normal and halfway vector)
    f_batch NdotH = dot_simd(hit_norm_x, hit_norm_y, hit_norm_z, 
                             half_x, half_y, half_z);

    // Step 4: Clamp to [0, infinity]
    f_batch zero_batch = xs::broadcast(0.0f);
    f_batch specular_base = xs::max(zero_batch, NdotH);

    // --- END BLINN-PHONG ---


    // Step 5: Load material specular colors and phong exponents (unfolded loop)
    alignas(32) fl mat_r[8], mat_g[8], mat_b[8], phong_exp[8];
    
    // Ray 0
    mat_r[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].specular_reflectance.x : 0.0f;
    mat_g[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].specular_reflectance.y : 0.0f;
    mat_b[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].specular_reflectance.z : 0.0f;
    phong_exp[0] = (ray_pack.mat_id[0] >= 0) ? scene.materials[ray_pack.mat_id[0]].phong_exponent : 1.0f;
    
    // Ray 1
    mat_r[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].specular_reflectance.x : 0.0f;
    mat_g[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].specular_reflectance.y : 0.0f;
    mat_b[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].specular_reflectance.z : 0.0f;
    phong_exp[1] = (ray_pack.mat_id[1] >= 0) ? scene.materials[ray_pack.mat_id[1]].phong_exponent : 1.0f;
    
    // Ray 2
    mat_r[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].specular_reflectance.x : 0.0f;
    mat_g[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].specular_reflectance.y : 0.0f;
    mat_b[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].specular_reflectance.z : 0.0f;
    phong_exp[2] = (ray_pack.mat_id[2] >= 0) ? scene.materials[ray_pack.mat_id[2]].phong_exponent : 1.0f;
    
    // Ray 3
    mat_r[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].specular_reflectance.x : 0.0f;
    mat_g[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].specular_reflectance.y : 0.0f;
    mat_b[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].specular_reflectance.z : 0.0f;
    phong_exp[3] = (ray_pack.mat_id[3] >= 0) ? scene.materials[ray_pack.mat_id[3]].phong_exponent : 1.0f;
    
    // Ray 4
    mat_r[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].specular_reflectance.x : 0.0f;
    mat_g[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].specular_reflectance.y : 0.0f;
    mat_b[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].specular_reflectance.z : 0.0f;
    phong_exp[4] = (ray_pack.mat_id[4] >= 0) ? scene.materials[ray_pack.mat_id[4]].phong_exponent : 1.0f;
    
    // Ray 5
    mat_r[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].specular_reflectance.x : 0.0f;
    mat_g[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].specular_reflectance.y : 0.0f;
    mat_b[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].specular_reflectance.z : 0.0f;
    phong_exp[5] = (ray_pack.mat_id[5] >= 0) ? scene.materials[ray_pack.mat_id[5]].phong_exponent : 1.0f;
    
    // Ray 6
    mat_r[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].specular_reflectance.x : 0.0f;
    mat_g[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].specular_reflectance.y : 0.0f;
    mat_b[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].specular_reflectance.z : 0.0f;
    phong_exp[6] = (ray_pack.mat_id[6] >= 0) ? scene.materials[ray_pack.mat_id[6]].phong_exponent : 1.0f;
    
    // Ray 7
    mat_r[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].specular_reflectance.x : 0.0f;
    mat_g[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].specular_reflectance.y : 0.0f;
    mat_b[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].specular_reflectance.z : 0.0f;
    phong_exp[7] = (ray_pack.mat_id[7] >= 0) ? scene.materials[ray_pack.mat_id[7]].phong_exponent : 1.0f;
    
    // Load into SIMD batches
    f_batch material_r = xs::load_aligned(mat_r);
    f_batch material_g = xs::load_aligned(mat_g);
    f_batch material_b = xs::load_aligned(mat_b);
    f_batch phong_exponent = xs::load_aligned(phong_exp);
    
    // Step 6: Calculate specular_base^phong_exponent using SIMD pow
    f_batch specular_intensity = xs::pow(specular_base, phong_exponent);
    
    // Step 7: Load light intensity as broadcast (same for all 8 rays)
    f_batch light_r = xs::broadcast(scene.point_light_data__.pl_intensity_r[light_index]);
    f_batch light_g = xs::broadcast(scene.point_light_data__.pl_intensity_g[light_index]);
    f_batch light_b = xs::broadcast(scene.point_light_data__.pl_intensity_b[light_index]);
    
    // Step 8: Apply distance attenuation (1 / distance^2)
    f_batch distance_squared = distance * distance;
    f_batch attenuation = xs::broadcast(1.0f) / distance_squared;
    
    // Step 9: Calculate specular color with attenuation (Blinn-Phong: I = ks * light * (N·H)^n / d^2)
    f_batch color_r = material_r * light_r * specular_intensity * attenuation;
    f_batch color_g = material_g * light_g * specular_intensity * attenuation;
    f_batch color_b = material_b * light_b * specular_intensity * attenuation;
    
    // Step 10: Apply shadow mask (only add color where not in shadow)
    color_r = xs::select(in_light, color_r, zero_batch);
    color_g = xs::select(in_light, color_g, zero_batch);
    color_b = xs::select(in_light, color_b, zero_batch);
    
    // Step 11: Store results back to color_block as FL (no conversion yet)
    alignas(32) fl final_r[8], final_g[8], final_b[8];
    xs::store_aligned(final_r, color_r);
    xs::store_aligned(final_g, color_g);
    xs::store_aligned(final_b, color_b);
    
    // Accumulate into ColorBlockFl (keep as fl)
    color_block.rgb[0]  += final_r[0];
    color_block.rgb[1]  += final_g[0];
    color_block.rgb[2]  += final_b[0];
    
    color_block.rgb[3]  += final_r[1];
    color_block.rgb[4]  += final_g[1];
    color_block.rgb[5]  += final_b[1];
    
    color_block.rgb[6]  += final_r[2];
    color_block.rgb[7]  += final_g[2];
    color_block.rgb[8]  += final_b[2];
    
    color_block.rgb[9]  += final_r[3];
    color_block.rgb[10] += final_g[3];
    color_block.rgb[11] += final_b[3];
    
    color_block.rgb[12] += final_r[4];
    color_block.rgb[13] += final_g[4];
    color_block.rgb[14] += final_b[4];
    
    color_block.rgb[15] += final_r[5];
    color_block.rgb[16] += final_g[5];
    color_block.rgb[17] += final_b[5];
    
    color_block.rgb[18] += final_r[6];
    color_block.rgb[19] += final_g[6];
    color_block.rgb[20] += final_b[6];
    
    color_block.rgb[21] += final_r[7];
    color_block.rgb[22] += final_g[7];
    color_block.rgb[23] += final_b[7];
}

// ------------------- FIXED FUNCTION -------------------
void inline shadow(RP8& ray_pack, const Scene& scene, xsimd::batch_bool<fl>& in_light, int light_index, const f_batch& light_distance){
    // Test for shadow ray intersection.
    // return_closest_hit normalizes the ray direction in ray_pack
    // internally, so t_mins will be the actual world-space distance.
    return_closest_hit(ray_pack, scene, scene.vertex_data__);
    
    // Load t_min values (this is the distance to the first hit)
    f_batch t_mins = xs::load(ray_pack.t_min);
    
    // A point is in light if the first intersection (t_mins)
    // is *farther away* than the light source (light_distance).
    // This means no object was hit between the surface and the light.
    in_light = t_mins > light_distance;
}
// ----------------- END FIXED FUNCTION -----------------

void inline mirror();


void inline shade(const RP8& ray_pack,const Scene& scene, ColorBlock& color_block){
    
    // Use ColorBlockFl for accumulation (fl values) - explicitly zero-initialize
    ColorBlockFl color_block_int = {};
    
    // Deconstruct incoming ray direction into fl batches
    f_batch ray_dir_x = xs::load(ray_pack.d_x);
    f_batch ray_dir_y = xs::load(ray_pack.d_y);
    f_batch ray_dir_z = xs::load(ray_pack.d_z);

    normalize_simd_overwrite(ray_dir_x,ray_dir_y,ray_dir_z);


    f_batch hit_x = xs::load(ray_pack.hit_pos_x); 
    f_batch hit_y = xs::load(ray_pack.hit_pos_y); 
    f_batch hit_z = xs::load(ray_pack.hit_pos_z);
    
    // Precalculated hit normal
    f_batch hit_norm_x = xs::load(ray_pack.hit_norm_x);
    f_batch hit_norm_y = xs::load(ray_pack.hit_norm_y);
    f_batch hit_norm_z = xs::load(ray_pack.hit_norm_z);
    
    // Handle background color for rays that missed all geometry
    for(int ray_idx = 0; ray_idx < 8; ray_idx++) {
        if(ray_pack.mat_id[ray_idx] < 0) {
            // Ray missed - set to background color
            color_block_int.rgb[ray_idx * 3 + 0] = scene.background_color.x;
            color_block_int.rgb[ray_idx * 3 + 1] = scene.background_color.y;
            color_block_int.rgb[ray_idx * 3 + 2] = scene.background_color.z;
        }
    }
    
    // Add ambient lighting first (independent of lights) - only for rays that hit geometry
    ambient(ray_pack, scene, color_block_int);
    
    for(int i = 0;i<scene.point_light_data__.pl_id.size();i++){
        b_batch in_light;
        
        // Precalculated light position
        f_batch light_pos_x = xs::broadcast(scene.point_light_data__.pl_pos_x[i]);
        f_batch light_pos_y = xs::broadcast(scene.point_light_data__.pl_pos_y[i]);
        f_batch light_pos_z = xs::broadcast(scene.point_light_data__.pl_pos_z[i]);
        
        // Calculate object to light vector
        f_batch object_to_light_x = light_pos_x - hit_x;
        f_batch object_to_light_y = light_pos_y - hit_y;
        f_batch object_to_light_z = light_pos_z - hit_z;
        f_batch distance = length_simd(object_to_light_x, object_to_light_y, object_to_light_z);
        
        // Create shadow rays with epsilon offset to avoid self-intersection
        RP8 shadow_rays;
        
        // Shadow ray direction is unnormalized (object to light)
        // return_closest_hit will normalize this vector internally
        xs::store(shadow_rays.d_x, object_to_light_x);
        xs::store(shadow_rays.d_y, object_to_light_y);
        xs::store(shadow_rays.d_z, object_to_light_z);
        
        // Offset shadow ray origin by epsilon along the normal to avoid shadow acne
        f_batch epsilon = xs::broadcast(scene.shadow_ray_epsilon);
        f_batch shadow_origin_x = hit_x + hit_norm_x * epsilon;
        f_batch shadow_origin_y = hit_y + hit_norm_y * epsilon;
        f_batch shadow_origin_z = hit_z + hit_norm_z * epsilon;
        
        xs::store(shadow_rays.o_x, shadow_origin_x);
        xs::store(shadow_rays.o_y, shadow_origin_y);
        xs::store(shadow_rays.o_z, shadow_origin_z);
        
        // Test shadow rays - pass the light distance for comparison
        shadow(shadow_rays, scene, in_light, i, distance);
        
        // If all rays are in shadow for this light, skip to next light
        if(xs::none(in_light)) continue;
        
        // Normalize to get light direction for shading calculations
        f_batch light_dir_x = object_to_light_x;
        f_batch light_dir_y = object_to_light_y;
        f_batch light_dir_z = object_to_light_z;
        normalize_simd_overwrite(light_dir_x, light_dir_y, light_dir_z);


        diffuse(ray_pack, scene, color_block_int, in_light, i,
                hit_x, hit_y, hit_z,
                hit_norm_x, hit_norm_y, hit_norm_z,
                light_pos_x, light_pos_y, light_pos_z,
                light_dir_x, light_dir_y, light_dir_z,
                ray_dir_x, ray_dir_y, ray_dir_z,
                distance);
        
        specular(ray_pack, scene, color_block_int, in_light, i,
                 hit_x, hit_y, hit_z,
                 hit_norm_x, hit_norm_y, hit_norm_z,
                 light_pos_x, light_pos_y, light_pos_z,
                 light_dir_x, light_dir_y, light_dir_z,
                 ray_dir_x, ray_dir_y, ray_dir_z,
                 distance);
        
        //mirror();
    }
    
    // Convert from fl to unsigned char with clamping at the very end
    convert_to_unsigned_char(color_block_int, color_block);
};