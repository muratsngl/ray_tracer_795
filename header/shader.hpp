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
    return_any_hit_shadow(ray_pack, scene, in_light,light_index,light_distance);
}
// ----------------- END FIXED FUNCTION -----------------

void inline mirror();


void inline shade_trad(const RP8& ray_pack,const Scene& scene, ColorBlock& color_block){

    
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



// --- Local SoA Ray Queue (Your "ray_data" struct) ---
// This will hold the refracted rays *locally* for one 2x4 pixel block.
struct SoARayQueue {
    std::vector<fl> o_x, o_y, o_z;
    std::vector<fl> d_x, d_y, d_z;
    std::vector<fl> tp_r, tp_g, tp_b;

    void push(fl ox, fl oy, fl oz, fl dx, fl dy, fl dz, fl tr, fl tg, fl tb) {
        o_x.push_back(ox); o_y.push_back(oy); o_z.push_back(oz);
        d_x.push_back(dx); d_y.push_back(dy); d_z.push_back(dz);
        tp_r.push_back(tr); tp_g.push_back(tg); tp_b.push_back(tb);
    }

    bool pop(fl& ox, fl& oy, fl& oz, fl& dx, fl& dy, fl& dz, fl& tr, fl& tg, fl& tb) {
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
        return true;
    }

    bool is_empty() {
        return o_x.empty();
    }
};

// --- Masked Local Shading Functions (NEW) ---
// These are the new local shaders. They now accept a MASK.

void inline ambient_masked(const RP8& ray_pack, const Scene& scene, ColorBlockFl& color_block, const b_batch& active_mask){
    
    // Step 1: Load material ambient colors (only for active rays)
    alignas(32) fl mat_r[8] = {0.f}, mat_g[8] = {0.f}, mat_b[8] = {0.f};
    
    for(int i = 0; i < 8; ++i) {
        if (active_mask.get(i)) {
            int mat_id = ray_pack.mat_id[i];
            if (mat_id >= 0) {
                mat_r[i] = scene.materials[mat_id].ambient_reflectance.x;
                mat_g[i] = scene.materials[mat_id].ambient_reflectance.y;
                mat_b[i] = scene.materials[mat_id].ambient_reflectance.z;
            }
        }
    }
    
    f_batch material_r = xs::load_aligned(mat_r);
    f_batch material_g = xs::load_aligned(mat_g);
    f_batch material_b = xs::load_aligned(mat_b);
    
    // Step 2: Load ambient light intensity
    f_batch ambient_r = xs::broadcast(scene.ambient_light.x);
    f_batch ambient_g = xs::broadcast(scene.ambient_light.y);
    f_batch ambient_b = xs::broadcast(scene.ambient_light.z);
    
    // Step 3: Calculate ambient color
    f_batch color_r = material_r * ambient_r;
    f_batch color_g = material_g * ambient_g;
    f_batch color_b = material_b * ambient_b;
    
    // Step 4: Accumulate into ColorBlockFl (scalar loop for masked add)
    alignas(32) fl final_r[8], final_g[8], final_b[8];
    xs::store_aligned(final_r, color_r);
    xs::store_aligned(final_g, color_g);
    xs::store_aligned(final_b, color_b);
    
    for(int i = 0; i < 8; ++i) {
        if (active_mask.get(i)) {
            color_block.rgb[i*3 + 0] += final_r[i];
            color_block.rgb[i*3 + 1] += final_g[i];
            color_block.rgb[i*3 + 2] += final_b[i];
        }
    }
}

void inline diffuse_masked(const RP8& ray_pack, const Scene& scene, ColorBlockFl& color_block, const b_batch& active_mask, int light_index,
                    const f_batch& light_dir_x, const f_batch& light_dir_y, const f_batch& light_dir_z,
                    const f_batch& distance){
    
    f_batch zero_batch = xs::broadcast(0.0f);
    if(xs::none(active_mask)) return;

    // Step 1: Calculate N·L
    f_batch NdotL = dot_simd(xs::load(ray_pack.hit_norm_x), xs::load(ray_pack.hit_norm_y), xs::load(ray_pack.hit_norm_z), 
                             light_dir_x, light_dir_y, light_dir_z);
    
    f_batch diffuse_intensity = xs::max(zero_batch, NdotL);
    
    // Step 2: Load material diffuse colors
    alignas(32) fl mat_r[8] = {0.f}, mat_g[8] = {0.f}, mat_b[8] = {0.f};
    for(int i = 0; i < 8; ++i) {
        if (active_mask.get(i)) {
            int mat_id = ray_pack.mat_id[i];
            if (mat_id >= 0) {
                mat_r[i] = scene.materials[mat_id].diffuse_reflectance.x;
                mat_g[i] = scene.materials[mat_id].diffuse_reflectance.y;
                mat_b[i] = scene.materials[mat_id].diffuse_reflectance.z;
            }
        }
    }
    f_batch material_r = xs::load_aligned(mat_r);
    f_batch material_g = xs::load_aligned(mat_g);
    f_batch material_b = xs::load_aligned(mat_b);
    
    // Step 3: Load light intensity
    f_batch light_r = xs::broadcast(scene.point_light_data__.pl_intensity_r[light_index]);
    f_batch light_g = xs::broadcast(scene.point_light_data__.pl_intensity_g[light_index]);
    f_batch light_b = xs::broadcast(scene.point_light_data__.pl_intensity_b[light_index]);
    
    // Step 4: Apply attenuation
    f_batch distance_squared = distance * distance;
    f_batch attenuation = xs::broadcast(1.0f) / distance_squared;
    
    // Step 5: Calculate diffuse color
    f_batch color_r = material_r * light_r * diffuse_intensity * attenuation;
    f_batch color_g = material_g * light_g * diffuse_intensity * attenuation;
    f_batch color_b = material_b * light_b * diffuse_intensity * attenuation;

    // Step 6: Accumulate into ColorBlockFl (scalar loop for masked add)
    alignas(32) fl final_r[8], final_g[8], final_b[8];
    xs::store_aligned(final_r, color_r);
    xs::store_aligned(final_g, color_g);
    xs::store_aligned(final_b, color_b);
    
    for(int i = 0; i < 8; ++i) {
        if (active_mask.get(i)) {
            color_block.rgb[i*3 + 0] += final_r[i];
            color_block.rgb[i*3 + 1] += final_g[i];
            color_block.rgb[i*3 + 2] += final_b[i];
        }
    }
}

void inline specular_masked(const RP8& ray_pack, const Scene& scene, ColorBlockFl& color_block, const b_batch& active_mask, int light_index,
                     const f_batch& light_dir_x, const f_batch& light_dir_y, const f_batch& light_dir_z,
                     const f_batch& ray_dir_x, const f_batch& ray_dir_y, const f_batch& ray_dir_z,
                     const f_batch& distance){
    
    f_batch zero_batch = xs::broadcast(0.0f);
    if(xs::none(active_mask)) return;

    // --- BLINN-PHONG ---
    f_batch view_x = -ray_dir_x;
    f_batch view_y = -ray_dir_y;
    f_batch view_z = -ray_dir_z;
    
    f_batch half_x = light_dir_x + view_x;
    f_batch half_y = light_dir_y + view_y;
    f_batch half_z = light_dir_z + view_z;
    normalize_simd_overwrite(half_x, half_y, half_z);
    
    f_batch NdotH = dot_simd(xs::load(ray_pack.hit_norm_x), xs::load(ray_pack.hit_norm_y), xs::load(ray_pack.hit_norm_z), 
                             half_x, half_y, half_z);
    f_batch specular_base = xs::max(zero_batch, NdotH);

    // --- Load material properties ---
    alignas(32) fl mat_r[8] = {0.f}, mat_g[8] = {0.f}, mat_b[8] = {0.f}, phong_exp[8] = {1.f};
    for(int i = 0; i < 8; ++i) {
        if (active_mask.get(i)) {
            int mat_id = ray_pack.mat_id[i];
            if (mat_id >= 0) {
                mat_r[i] = scene.materials[mat_id].specular_reflectance.x;
                mat_g[i] = scene.materials[mat_id].specular_reflectance.y;
                mat_b[i] = scene.materials[mat_id].specular_reflectance.z;
                phong_exp[i] = scene.materials[mat_id].phong_exponent;
            }
        }
    }
    f_batch material_r = xs::load_aligned(mat_r);
    f_batch material_g = xs::load_aligned(mat_g);
    f_batch material_b = xs::load_aligned(mat_b);
    f_batch phong_exponent = xs::load_aligned(phong_exp);
    
    // --- Calculate Color ---
    f_batch specular_intensity = xs::pow(specular_base, phong_exponent);
    
    f_batch light_r = xs::broadcast(scene.point_light_data__.pl_intensity_r[light_index]);
    f_batch light_g = xs::broadcast(scene.point_light_data__.pl_intensity_g[light_index]);
    f_batch light_b = xs::broadcast(scene.point_light_data__.pl_intensity_b[light_index]);
    
    f_batch distance_squared = distance * distance;
    f_batch attenuation = xs::broadcast(1.0f) / distance_squared;
    
    f_batch color_r = material_r * light_r * specular_intensity * attenuation;
    f_batch color_g = material_g * light_g * specular_intensity * attenuation;
    f_batch color_b = material_b * light_b * specular_intensity * attenuation;
    
    // *** FIX: Apply active_mask to zero out inactive/shadowed rays ***
    color_r = xs::select(active_mask, color_r, zero_batch);
    color_g = xs::select(active_mask, color_g, zero_batch);
    color_b = xs::select(active_mask, color_b, zero_batch);
    
    // --- Accumulate ---
    alignas(32) fl final_r[8], final_g[8], final_b[8];
    xs::store_aligned(final_r, color_r);
    xs::store_aligned(final_g, color_g);
    xs::store_aligned(final_b, color_b);
    
    for(int i = 0; i < 8; ++i) {
        if (active_mask.get(i)) {
            color_block.rgb[i*3 + 0] += final_r[i];
            color_block.rgb[i*3 + 1] += final_g[i];
            color_block.rgb[i*3 + 2] += final_b[i];
        }
    }
}

// Helper to add masked color to the final buffer
void inline accumulate_color(ColorBlockFl& color_block, const b_batch& mask,
                             const f_batch& r, const f_batch& g, const f_batch& b) {
    
    alignas(32) fl add_r[8], add_g[8], add_b[8];
    xs::store_aligned(add_r, r);
    xs::store_aligned(add_g, g);
    xs::store_aligned(add_b, b);

    for (int i = 0; i < 8; i++) {
        if (mask.get(i)) {
            color_block.rgb[i*3 + 0] += add_r[i];
            color_block.rgb[i*3 + 1] += add_g[i];
            color_block.rgb[i*3 + 2] += add_b[i];
        }
    }
}

// Helper to do local shading (ambient + all lights)
void inline shade_local_masked(RP8& ray_pack, const Scene& scene, const b_batch& local_mask,
                               ColorBlockFl& local_color_buffer,
                               const f_batch& ray_dir_x, const f_batch& ray_dir_y, const f_batch& ray_dir_z) {
    
    if (xs::none(local_mask)) return;

    // --- 1. Add Ambient Light ---
    ambient_masked(ray_pack, scene, local_color_buffer, local_mask);

    // --- 2. Add Light Contributions ---
    f_batch hit_x = xs::load(ray_pack.hit_pos_x); 
    f_batch hit_y = xs::load(ray_pack.hit_pos_y); 
    f_batch hit_z = xs::load(ray_pack.hit_pos_z);
    f_batch hit_norm_x = xs::load(ray_pack.hit_norm_x);
    f_batch hit_norm_y = xs::load(ray_pack.hit_norm_y);
    f_batch hit_norm_z = xs::load(ray_pack.hit_norm_z);

    for(int i = 0; i < scene.point_light_data__.pl_id.size(); i++){
        
        f_batch light_pos_x = xs::broadcast(scene.point_light_data__.pl_pos_x[i]);
        f_batch light_pos_y = xs::broadcast(scene.point_light_data__.pl_pos_y[i]);
        f_batch light_pos_z = xs::broadcast(scene.point_light_data__.pl_pos_z[i]);
        
        f_batch object_to_light_x = light_pos_x - hit_x;
        f_batch object_to_light_y = light_pos_y - hit_y;
        f_batch object_to_light_z = light_pos_z - hit_z;
        f_batch distance = length_simd(object_to_light_x, object_to_light_y, object_to_light_z);
        
        // --- 3. Shadow Test ---
        RP8 shadow_rays;
        f_batch epsilon = xs::broadcast(scene.shadow_ray_epsilon);
        xs::store(shadow_rays.o_x, hit_x + hit_norm_x * epsilon);
        xs::store(shadow_rays.o_y, hit_y + hit_norm_y * epsilon);
        xs::store(shadow_rays.o_z, hit_z + hit_norm_z * epsilon);
        xs::store(shadow_rays.d_x, object_to_light_x);
        xs::store(shadow_rays.d_y, object_to_light_y);
        xs::store(shadow_rays.d_z, object_to_light_z);

        b_batch in_light;
        return_any_hit_shadow_masked(shadow_rays, scene, in_light, i, distance, local_mask);
        
        b_batch final_shade_mask = in_light & local_mask;
        if(xs::none(final_shade_mask)) continue;

        // --- 4. Add Diffuse & Specular ---
        f_batch light_dir_x = object_to_light_x;
        f_batch light_dir_y = object_to_light_y;
        f_batch light_dir_z = object_to_light_z;
        normalize_simd_overwrite(light_dir_x, light_dir_y, light_dir_z);

        diffuse_masked(ray_pack, scene, local_color_buffer, final_shade_mask, i,
                       light_dir_x, light_dir_y, light_dir_z, distance);
        
        specular_masked(ray_pack, scene, local_color_buffer, final_shade_mask, i,
                        light_dir_x, light_dir_y, light_dir_z,
                        ray_dir_x, ray_dir_y, ray_dir_z, distance);
    }
}


// --- THE NEW RECURSIVE SHADER ---
void inline shade_recursive(RP8& ray_pack, const Scene& scene, ColorBlock& color_block){
    
    // Final pixel colors are accumulated here
    ColorBlockFl final_color_buffer = {};
    
    // This is the "throughput" (color filter) for each ray
    alignas(32) fl tp_r[8] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
    alignas(32) fl tp_g[8] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
    alignas(32) fl tp_b[8] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
    
    f_batch zero_batch = xs::broadcast(0.f);
    f_batch one_batch = xs::broadcast(1.f);
    

    // This is your "hashmap" - a local queue for this 8-pixel block
    SoARayQueue refracted_queue;
    
    // This is your "indices" array, mapping RP8 lane to ColorBlock pixel (0-7)
    // For this architecture, it's just 0-7 and never changes.
    // We don't even need the array, we just use the loop index 'i'.

    // The "bitmask" of which rays are still active and bouncing
    b_batch active_mask = true;

    // This is the main "bounce" loop
    for(int depth = 0; depth <= scene.max_recursion_depth; ++depth) {
        
        // 1. BACK-FILLING (Your design)
        // If some lanes are dead, fill them from the queue
        if (!xs::all(active_mask) && !refracted_queue.is_empty()) {
            alignas(32) fl active_fl[8];
            xs::store(active_fl, xs::select(active_mask, one_batch, zero_batch));

            for(int i = 0; i < 8; ++i) {
                if(active_fl[i] == 0.0f) { // This lane is dead, fill it
                    if (refracted_queue.is_empty()) break; // Queue is empty, stop filling
                    
                    // Pop ray from queue
                    ray_pack.o_x[i] = refracted_queue.o_x.back(); refracted_queue.o_x.pop_back();
                    ray_pack.o_y[i] = refracted_queue.o_y.back(); refracted_queue.o_y.pop_back();
                    ray_pack.o_z[i] = refracted_queue.o_z.back(); refracted_queue.o_z.pop_back();
                    ray_pack.d_x[i] = refracted_queue.d_x.back(); refracted_queue.d_x.pop_back();
                    ray_pack.d_y[i] = refracted_queue.d_y.back(); refracted_queue.d_y.pop_back();
                    ray_pack.d_z[i] = refracted_queue.d_z.back(); refracted_queue.d_z.pop_back();
                    tp_r[i] = refracted_queue.tp_r.back(); refracted_queue.tp_r.pop_back();
                    tp_g[i] = refracted_queue.tp_g.back(); refracted_queue.tp_g.pop_back();
                    tp_b[i] = refracted_queue.tp_b.back(); refracted_queue.tp_b.pop_back();
                    // We must re-activate this lane
                    active_fl[i] = 1.0f;
                }
            }
            active_mask = xs::load_aligned(active_fl) > 0.0f;
        }

        // If all 8 lanes are dead AND the queue is empty, we are done.
        if (xs::none(active_mask) && refracted_queue.is_empty()) {
            break;
        }

        // 2. TRACE
        // Reset t_min and mat_id for active rays
        f_batch t_inf = xs::broadcast(__builtin_inff());
        i_batch id_none = xs::broadcast(-1);
        xs::store_aligned(ray_pack.t_min, t_inf);
        xs::store_aligned(ray_pack.mat_id, id_none);
        
        return_closest_hit_masked(ray_pack, scene, scene.vertex_data__, active_mask);

        // 3. SCALAR DISPATCH LOOP (Your design)
        f_batch throughput_r = xs::load_aligned(tp_r);
        f_batch throughput_g = xs::load_aligned(tp_g);
        f_batch throughput_b = xs::load_aligned(tp_b);
        f_batch current_throughput_r = xs::load_aligned(tp_r); // Load current throughput
        f_batch current_throughput_g = xs::load_aligned(tp_g);
        f_batch current_throughput_b = xs::load_aligned(tp_b);
        f_batch zero_batch = xs::broadcast(0.0f);

        // Load all hit data now for SIMD functions
        f_batch ray_dir_x = xs::load(ray_pack.d_x);
        f_batch ray_dir_y = xs::load(ray_pack.d_y);
        f_batch ray_dir_z = xs::load(ray_pack.d_z);
        // *** FIX: Normalize ray directions for view vector calculation ***
        normalize_simd_overwrite(ray_dir_x, ray_dir_y, ray_dir_z);
        
        f_batch hit_norm_x = xs::load(ray_pack.hit_norm_x);
        f_batch hit_norm_y = xs::load(ray_pack.hit_norm_y);
        f_batch hit_norm_z = xs::load(ray_pack.hit_norm_z);
        f_batch hit_pos_x = xs::load(ray_pack.hit_pos_x);
        f_batch hit_pos_y = xs::load(ray_pack.hit_pos_y);
        f_batch hit_pos_z = xs::load(ray_pack.hit_pos_z);
        
        // These masks will track what to do *after* the scalar loop
        b_batch next_active_mask = b_batch(false);
        b_batch local_shade_mask = b_batch(false);
        
        // Use boolean arrays for easier modification
        alignas(32) bool next_active[8] = {false};
        alignas(32) bool local_shade[8] = {false};

        alignas(32) fl next_tp_r[8];
        alignas(32) fl next_tp_g[8];
        alignas(32) fl next_tp_b[8];
        // Initialize next_tp with current tp (rays that don't hit mirrors keep their tp)
        xs::store_aligned(next_tp_r, current_throughput_r);
        xs::store_aligned(next_tp_g, current_throughput_g);
        xs::store_aligned(next_tp_b, current_throughput_b);
        
        for(int i = 0; i < 8; ++i) {
            if (!active_mask.get(i)) continue; // Skip inactive rays

            int mat_id = ray_pack.mat_id[i];

            if (mat_id < 0) {
                // --- Ray Missed (Hit Background) ---
                final_color_buffer.rgb[i*3 + 0] += scene.background_color.x * tp_r[i];
                final_color_buffer.rgb[i*3 + 1] += scene.background_color.y * tp_g[i];
                final_color_buffer.rgb[i*3 + 2] += scene.background_color.z * tp_b[i];
                // Ray is now inactive
            
            } else {
                const Material& mat = scene.materials[mat_id];

                if (mat.type == "mirror") {
                    // --- Ray Hit Mirror ---
                    // Don't shade. Update throughput, reflect ray, keep active.
                    next_tp_r[i] = tp_r[i] * mat.mirror_reflectance.x;
                    next_tp_g[i] = tp_g[i] * mat.mirror_reflectance.y;
                    next_tp_b[i] = tp_b[i] * mat.mirror_reflectance.z;

                    // Reflect: wr = -wo + 2n(n.wo)
                    fl cos_i = dot_scalar(ray_dir_x.get(i), ray_dir_y.get(i), ray_dir_z.get(i),
                                          hit_norm_x.get(i), hit_norm_y.get(i), hit_norm_z.get(i));
                    
                    ray_pack.d_x[i] = ray_dir_x.get(i) - 2.0f * hit_norm_x.get(i) * cos_i;
                    ray_pack.d_y[i] = ray_dir_y.get(i) - 2.0f * hit_norm_y.get(i) * cos_i;
                    ray_pack.d_z[i] = ray_dir_z.get(i) - 2.0f * hit_norm_z.get(i) * cos_i;
                    
                    ray_pack.o_x[i] = hit_pos_x.get(i) + hit_norm_x.get(i) * scene.shadow_ray_epsilon;
                    ray_pack.o_y[i] = hit_pos_y.get(i) + hit_norm_y.get(i) * scene.shadow_ray_epsilon;
                    ray_pack.o_z[i] = hit_pos_z.get(i) + hit_norm_z.get(i) * scene.shadow_ray_epsilon;

                   
                   
                   

                    local_shade[i]= true;
                    next_active[i] = true; // Keep this lane active
                    

                } else if (mat.type == "dielectric") {
                    // // --- Ray Hit Dielectric ---
                    // // TODO: This is where the complex logic goes.
                    // // For now, treat as a mirror.
                    // tp_r[i] *= 0.8f; // Placeholder
                    // tp_g[i] *= 0.8f;
                    // tp_b[i] *= 0.8f;

                    // fl cos_i = dot_scalar(ray_dir_x.get(i), ray_dir_y.get(i), ray_dir_z.get(i),
                    //                       hit_norm_x.get(i), hit_norm_y.get(i), hit_norm_z.get(i));
                    
                    // ray_pack.d_x[i] = ray_dir_x.get(i) - 2.0f * hit_norm_x.get(i) * cos_i;
                    // ray_pack.d_y[i] = ray_dir_y.get(i) - 2.0f * hit_norm_y.get(i) * cos_i;
                    // ray_pack.d_z[i] = ray_dir_z.get(i) - 2.0f * hit_norm_z.get(i) * cos_i;
                    
                    // ray_pack.o_x[i] = hit_pos_x.get(i) + hit_norm_x.get(i) * scene.shadow_ray_epsilon;
                    // ray_pack.o_y[i] = hit_pos_y.get(i) + hit_norm_y.get(i) * scene.shadow_ray_epsilon;
                    // ray_pack.o_z[i] = hit_pos_z.get(i) + hit_norm_z.get(i) * scene.shadow_ray_epsilon;

                    // next_active[i] = true; // Keep this lane active

                } else {
                    // --- Ray Hit Diffuse/Default ---
                    // Mark for local shading. Ray will become inactive.
                    local_shade[i] = true;
                }
            }
        } // End scalar dispatch loop
        
        // Load boolean arrays into batch_bool masks
        next_active_mask = b_batch::load_aligned(next_active);
        local_shade_mask = b_batch::load_aligned(local_shade);

        // 4. SHADE LOCAL RAYS (in SIMD)
        // We do this *after* the loop for max SIMD efficiency
        if(xs::any(local_shade_mask)) {
            ColorBlockFl local_color_buffer = {};
            
            // Run all local shading (ambient, diffuse, specular)
            shade_local_masked(ray_pack, scene, local_shade_mask, local_color_buffer,
                               ray_dir_x, ray_dir_y, ray_dir_z);
            
            // Accumulate the final, throughput-modulated color
            alignas(32) fl local_r[8], local_g[8], local_b[8];
            for(int i=0; i<8; i++) {
                local_r[i] = local_color_buffer.rgb[i*3+0] * tp_r[i];
                local_g[i] = local_color_buffer.rgb[i*3+1] * tp_g[i];
                local_b[i] = local_color_buffer.rgb[i*3+2] * tp_b[i];

                

            }
            
            accumulate_color(final_color_buffer, local_shade_mask,
                             xs::load_aligned(local_r),
                             xs::load_aligned(local_g),
                             xs::load_aligned(local_b));
        }

        // 5. PREPARE FOR NEXT LOOP
        active_mask = next_active_mask; // Only mirror/dielectric rays survive
        
        xs::store_aligned(tp_r, xs::load_aligned(next_tp_r));
        xs::store_aligned(tp_g, xs::load_aligned(next_tp_g));
        xs::store_aligned(tp_b, xs::load_aligned(next_tp_b));
    } // End bounce loop

    // Convert from fl to unsigned char with clamping at the very end
    convert_to_unsigned_char(final_color_buffer, color_block);
}

// --- LEGACY SHADE FUNCTION (for compatibility) ---
// We keep this, but the new code should call shade_recursive
void inline shade(const RP8& ray_pack,const Scene& scene, ColorBlock& color_block){
    
    // This function is now just a wrapper for the new recursive shader
    // We need a mutable copy of the ray_pack to pass to the recursive function
    RP8 mutable_ray_pack = ray_pack;
    shade_recursive(mutable_ray_pack, scene, color_block);
}
