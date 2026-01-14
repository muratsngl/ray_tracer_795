
# ray_tracer::devlog_5
## Part 1: Breaking Free from 8-bit Hell

The Dual-Buffer Architecture

My first challenge was moving away from the limited [0, 255] range. I needed to store floating-point radiance values throughout the rendering pipeline without clipping or quantization. Here's the image structure I designed:

```cpp
struct Image {
    int id;
    std::string path;
    int width, height, channels;
    std::vector<unsigned char> data;      // LDR path (existing)
    std::vector<float> data_float;        // NEW: HDR storage
    bool is_hdr = false;                  // Runtime discrimination
};
```

Why dual storage? I wanted to maintain backward compatibility with my existing LDR textures while supporting HDR environment maps. The is_hdr flag lets me dispatch to the correct buffer during texture sampling—no code duplication needed.

Tone Mapping Infrastructure

I extended my Camera structure to support multiple tone mapping operators per render. This was crucial because I wanted to generate multiple output images (different exposures, gamma values) from a single expensive HDR render:

```cpp
struct Camera {
    // ... existing fields ...
    std:: vector<Tonemap> tonemaps;    // Batch processing support
    bool is_hdr = false;              // Output format control
};

struct Tonemap {
    std::string tmo;           // "Photographic" | "Filmic" | "ACES"
    fl key;                    // 'a':  luminance key value
    fl burn_out;               // Percentile for L_white threshold
    fl saturation;             // 's': color reconstruction param
    fl gamma;                  // 'g':  gamma correction (applied as 1/g)
    std::string extension;     // Output filename suffix
};
```

The tonemaps vector allows me to process one HDR framebuffer multiple times with different parameters—effectively batch rendering without re-tracing rays.

## Part 2: Photographic Tone Mapping—Getting the Math Right

Implementing Reinhard's photographic operator was trickier than I expected. 

Luminance Conversion

I used the standard coefficients for RGB-to-luminance conversion:

```text
Y = 0.2126·R + 0.7152·G + 0.0722·B
```


The Burn-Out Logic

The burn_out parameter specifies a percentile (e.g., 0.95 = 95th percentile). Here's my approach:

Collect all luminance values from the HDR framebuffer
Sort the array (yes, expensive—but only done once per image)
Find L_white as the value at index burn_out × num_pixels
Apply the modified photographic operator:

```cpp
L_d = L(x,y) * (a / L_avg) / (1 + L(x,y) * (a / L_avg));
```

where L_white acts as a threshold for selective highlight compression.

Color Reconstruction & Gamma

After tone mapping luminance, I reconstructed color using the saturation parameter:

```cpp
C_out = pow(C_in / L_in, s) * L_out;
```

Then applied gamma correction:

```cpp
C_final = pow(C_out, 1.0f / gamma);
```

before quantizing to [0, 255]. This two-stage process preserves hue while controlling saturation in the final output.

## Part 3: Analytic Lights—Beyond Point Sources

I needed to support directional and spot lights for more realistic lighting scenarios. Both required careful SIMD implementation to maintain performance.

Directional Lights: Modeling the Sun

Directional lights simulate infinitely distant sources (like the sun) with parallel rays. The key insight: no distance-based attenuation.

```cpp
void inline directional_masked(const RP8& ray_pack, const Scene& scene, 
                               ColorBlockFl& color_block, const b_batch& active_mask, int i,
                               const f_batch& ray_dir_x, const f_batch& ray_dir_y, const f_batch& ray_dir_z,
                               const fl* diff_r, const fl* diff_g, const fl* diff_b,
                               const fl* spec_r, const fl* spec_g, const fl* spec_b, const fl* ph_exp) {
    
    // Light direction (negated—points TOWARD light)
    f_batch light_dir_x = xs::broadcast(-scene.directional_light_data__.dl_dir_x[i]);
    f_batch light_dir_y = xs::broadcast(-scene.directional_light_data__.dl_dir_y[i]);
    f_batch light_dir_z = xs::broadcast(-scene.directional_light_data__.dl_dir_z[i]);
    
    // Shadow ray with epsilon offset to prevent self-intersection
    RP8 shadow_rays;
    f_batch epsilon = xs::broadcast(scene.shadow_ray_epsilon);
    xs::store(shadow_rays.o_x, xs::load(ray_pack.hit_pos_x) + xs::load(ray_pack.hit_norm_x) * epsilon);
    xs::store(shadow_rays.o_y, xs::load(ray_pack.hit_pos_y) + xs::load(ray_pack.hit_norm_y) * epsilon);
    xs::store(shadow_rays.o_z, xs::load(ray_pack.hit_pos_z) + xs::load(ray_pack.hit_norm_z) * epsilon);
    xs::store(shadow_rays.d_x, light_dir_x);
    xs::store(shadow_rays.d_y, light_dir_y);
    xs::store(shadow_rays.d_z, light_dir_z);
    
    // Any-hit shadow test (terminates on first occlusion)
    b_batch is_occluded;
    intersect_tlas_any_hit_wrapper(shadow_rays, scene, xs::broadcast(__builtin_inff()), 
                                   is_occluded, active_mask);
    b_batch in_light = active_mask & (! is_occluded);
    if (xs::none(in_light)) return;

    // Lambert's cosine law: N·L
    f_batch norm_x = xs::load(ray_pack.hit_norm_x);
    f_batch norm_y = xs::load(ray_pack.hit_norm_y);
    f_batch norm_z = xs::load(ray_pack.hit_norm_z);
    f_batch NdotL = xs::max(f_batch(0.0f), dot_simd(norm_x, norm_y, norm_z, 
                                                     light_dir_x, light_dir_y, light_dir_z));
    
    // Irradiance E = Radiance (constant—no 1/r² term)
    f_batch rad_r = xs::broadcast(scene.directional_light_data__.dl_radiance_r[i]);
    f_batch rad_g = xs::broadcast(scene.directional_light_data__.dl_radiance_g[i]);
    f_batch rad_b = xs::broadcast(scene.directional_light_data__.dl_radiance_b[i]);

    // Diffuse:  kd · E · cos(θ)
    f_batch diffuse_r = xs::load_aligned(diff_r) * rad_r * NdotL;
    f_batch diffuse_g = xs::load_aligned(diff_g) * rad_g * NdotL;
    f_batch diffuse_b = xs::load_aligned(diff_b) * rad_b * NdotL;

    // Specular (Blinn-Phong): ks · E · (N·H)^p
    f_batch half_x = light_dir_x - ray_dir_x;
    f_batch half_y = light_dir_y - ray_dir_y;
    f_batch half_z = light_dir_z - ray_dir_z;
    normalize_simd_overwrite(half_x, half_y, half_z);
    f_batch spec_int = xs::pow(xs::max(f_batch(0.0f), 
                                       dot_simd(norm_x, norm_y, norm_z, half_x, half_y, half_z)), 
                               xs::load_aligned(ph_exp));
    
    f_batch specular_r = xs::load_aligned(spec_r) * rad_r * spec_int;
    f_batch specular_g = xs::load_aligned(spec_g) * rad_g * spec_int;
    f_batch specular_b = xs:: load_aligned(spec_b) * rad_b * spec_int;

    accumulate_color(color_block, in_light, diffuse_r + specular_r, 
                                            diffuse_g + specular_g, 
                                            diffuse_b + specular_b);
}
```

I trace shadow rays to infinity (t_max = inf) because directional lights have no position—only direction. This simplified the occlusion test significantly.

Spot Lights: Three-Zone Falloff

Spot lights were more complex. I needed to implement a smooth falloff between the inner (no-falloff) cone and outer (coverage) cone. Here's my three-zone attenuation logic:

```cpp
// α: angle between ray FROM light TO hit point and spotlight axis
f_batch cos_alpha = dot_simd(-L_dir_x, -L_dir_y, -L_dir_z, spot_dir_x, spot_dir_y, spot_dir_z);
f_batch cos_coverage_half = xs::cos(xs::broadcast(scene.spot_light_data__.sl_coverage[i] * 0.5f));
f_batch cos_falloff_half = xs::cos(xs::broadcast(scene.spot_light_data__.sl_falloff[i] * 0.5f));

// Zone 1: Inside inner cone (α < θ_f/2) → s = 1
// Zone 2: Falloff region (θ_f/2 ≤ α ≤ θ_c/2) → s = smooth interpolation
// Zone 3: Outside coverage (α > θ_c/2) → s = 0

b_batch in_no_falloff = (cos_alpha >= cos_falloff_half);
f_batch s_numerator = cos_alpha - cos_coverage_half;
f_batch s_denominator = cos_falloff_half - cos_coverage_half;
f_batch s_ratio = s_numerator / xs::max(s_denominator, xs::broadcast(1e-6f));
f_batch s_falloff = xs::pow(xs::max(f_batch(0.0f), xs::min(s_ratio, f_batch(1.0f))), 
                            f_batch(4.0f));  // Quartic smoothstep
f_batch s = xs::select(in_no_falloff, xs::broadcast(1.0f), s_falloff);
```


```cpp
s = pow((cos(alpha) - cos(theta_c/2)) / (cos(theta_f/2) - cos(theta_c/2)), 4);
```


## Part 4: Environment Lighting—The Real Challenge

Implementing environment map sampling was the most mathematically intensive part. I needed to handle spherical sampling, UV mapping, and Monte Carlo integration correctly.

Inversion Sampling: Two Approaches

I implemented both uniform and cosine-weighted hemisphere sampling via analytical CDF inversion.

Uniform Hemisphere Sampling

```cpp
inline void sample_uniform_hemisphere(fl xi1, fl xi2, fl& lx, fl& ly, fl& lz) {
    fl phi = 2.0f * M_PI * xi1;            // Azimuth:  φ = 2π·ξ₁
    fl cos_theta = xi2;                     // Direct inversion:  cos(θ) = ξ₂
    fl sin_theta = std::sqrt(1.0f - cos_theta * cos_theta);
    
    lx = sin_theta * std::cos(phi);
    ly = sin_theta * std::sin(phi);
    lz = cos_theta;                         // cos(θ) in local space
}
```

PDF: p(ω) = 1/(2π)

This generates uniform samples over the hemisphere. Simple, but not optimal for diffuse surfaces.

Cosine-Weighted Sampling

```cpp
inline void sample_cosine_hemisphere(fl xi1, fl xi2, fl& lx, fl& ly, fl& lz) {
    fl phi = 2.0f * M_PI * xi1;
    fl theta = std::acos(std::sqrt(xi2));   // KEY:  θ = arccos(√ξ₂)
    
    lx = std::sin(theta) * std::cos(phi);
    ly = std::sin(theta) * std::sin(phi);
    lz = std::cos(theta);                   // cos(θ) = √ξ₂
}
```

PDF: p(ω) = cos(θ)/π

This concentrates samples near the normal direction—perfect for importance sampling Lambertian BRDFs.

Building an Orthonormal Basis

To use these samples, I needed to transform from tangent space (where +Z = normal) to world space. This required constructing an orthonormal frame:

```cpp
inline void build_onb(fl nx, fl ny, fl nz, 
                     fl& tx, fl& ty, fl& tz,
                     fl& bx, fl& by, fl& bz) {
    // Choose auxiliary vector not parallel to N
    fl ux, uy, uz;
    if (std::abs(nx) > 0.9f) {
        ux = 0.0f; uy = 1.0f; uz = 0.0f;    // Use Y if N ≈ X-axis
    } else {
        ux = 1.0f; uy = 0.0f; uz = 0.0f;    // Use X otherwise
    }
    
    // T = normalize(N × up)
    tx = ny * uz - nz * uy;
    ty = nz * ux - nx * uz;
    tz = nx * uy - ny * ux;
    fl t_len = std::sqrt(tx*tx + ty*ty + tz*tz);
    tx /= t_len; ty /= t_len; tz /= t_len;
    
    // B = T × N (already normalized)
    bx = ty * nz - tz * ny;
    by = tz * nx - tx * nz;
    bz = tx * ny - ty * nx;
}
```

The auxiliary vector selection (X or Y axis) prevents numerical instability when the normal is nearly parallel to one axis. I check if |nx| > 0.9 to decide which axis to use.

UV Mapping: Two Projections

I implemented both latitude-longitude (equirectangular) and light probe (angular) mappings.

Latitude-Longitude Mapping

```cpp
void inline get_latlong_uv(const f_batch& dx, const f_batch& dy, const f_batch& dz, 
                           f_batch& u, f_batch& v) {
    f_batch pi_batch = xs::broadcast(static_cast<fl>(M_PI));
    f_batch two_pi = xs::broadcast(static_cast<fl>(2.0 * M_PI));
    
    // u = 0.5 + atan2(dx, -dz) / (2π)  →  Azimuth angle
    u = xs::broadcast(0.5f) + xs::atan2(dx, -dz) / two_pi;
    
    // v = acos(dy) / π  →  Polar angle
    v = xs::acos(xs::min(xs::max(dy, f_batch(-1.0f)), f_batch(1.0f))) / pi_batch;
}
```

Why clamp dy? Floating-point errors can produce |dy| > 1, causing acos to return NaN. Clamping to [-1, 1] prevents this.

Light Probe (Angular) Mapping

```cpp
void inline get_probe_uv(const f_batch& dx, const f_batch& dy, const f_batch& dz, 
                         f_batch& u, f_batch& v) {
    f_batch dx_sq_plus_dy_sq = dx * dx + dy * dy;
    f_batch sqrt_term = xs::sqrt(dx_sq_plus_dy_sq);
    f_batch acos_term = xs::acos(xs::min(xs::max(-dz, f_batch(-1.0f)), f_batch(1.0f)));
    
    f_batch pi_batch = xs::broadcast(static_cast<fl>(M_PI));
    
    // Handle singularity at pole (dx ≈ dy ≈ 0)
    b_batch at_pole = (sqrt_term < xs::broadcast(1e-6f));
    f_batch safe_sqrt = xs::select(at_pole, xs::broadcast(1.0f), sqrt_term);
    
    // r = θ / (π · √(dx² + dy²))
    f_batch r = acos_term / (pi_batch * safe_sqrt);
    
    // u = dx·r + 0.5, v = -dy·r + 0.5
    u = xs::select(at_pole, xs::broadcast(0.5f), dx * r + xs::broadcast(0.5f));
    v = xs::select(at_pole, xs::broadcast(0.5f), -dy * r + xs::broadcast(0.5f));
}
```

The Pole Problem: When the direction points directly at -Z (dx ≈ dy ≈ 0), we'd divide by zero. I detect this condition and explicitly set u = v = 0.5 (texture center) to avoid NaN propagation.

Monte Carlo Integration: Getting PDF Compensation Right

This was the subtlest part. I needed to correctly weight samples to ensure unbiased radiance estimates. The rendering equation for diffuse surfaces with environment lighting:

```text
L_o = ∫_Ω (k_d / π) · L_env(ω) · cos(θ) dω
```

Using Monte Carlo integration with N samples:

```text
L_o ≈ (1/N) Σ [ (k_d / π) · L_env(ω_i) · cos(θ_i) / p(ω_i) ]
```

The weight factor w = cos(θ) / p(ω) depends on the sampling strategy:

```cpp
fl cos_theta = lz;  // cos(θ) in local tangent space where +Z = normal

if (sampler_type == "cosine") {
    // PDF: p(ω) = cos(θ) / π
    // Weight: w = cos(θ) / (cos(θ)/π) = π
    // But (k_d/π) · π = k_d, so: 
    fl weight = 1.0f;  // Factors cancel perfectly! 
    env_contrib_r[lane] = diffuse_r[lane] * env_r.get(0) * weight;
    env_contrib_g[lane] = diffuse_g[lane] * env_g.get(0) * weight;
    env_contrib_b[lane] = diffuse_b[lane] * env_b.get(0) * weight;
    
} else {  // "uniform"
    // PDF: p(ω) = 1 / (2π)
    // Weight: w = cos(θ) / (1/(2π)) = 2π · cos(θ)
    // Combined: (k_d/π) · 2π · cos(θ) = 2 · k_d · cos(θ)
    fl weight = 2.0f * cos_theta;
    env_contrib_r[lane] = diffuse_r[lane] * env_r.get(0) * weight;
    env_contrib_g[lane] = diffuse_g[lane] * env_g.get(0) * weight;
    env_contrib_b[lane] = diffuse_b[lane] * env_b.get(0) * weight;
}
```

With cosine-weighted sampling, the PDF exactly matches the cos(θ) term in the integrand, causing perfect cancellation. This is why importance sampling is so powerful—it reduces varianc[...]

## Part 5: Putting It All Together

Here's how I integrated environment lighting into my shader:

```cpp
// --- 5. Environment Map as Light Source with Importance Sampling ---
if (! scene.env_light_data__.el_image_id.empty()) {
    std::string sampler_type = scene.env_light_data__.el_sampler.empty() ? 
                               "uniform" : scene.env_light_data__.el_sampler[0];
    
    alignas(32) fl env_contrib_r[8] = {0};
    alignas(32) fl env_contrib_g[8] = {0};
    alignas(32) fl env_contrib_b[8] = {0};
    
    for (int lane = 0; lane < 8; ++lane) {
        if (! local_mask.get(lane)) continue;
        
        // Get surface normal
        fl nx = ray_pack.hit_norm_x[lane];
        fl ny = ray_pack.hit_norm_y[lane];
        fl nz = ray_pack.hit_norm_z[lane];
        
        // Build orthonormal basis
        fl tx, ty, tz, bx, by, bz;
        build_onb(nx, ny, nz, tx, ty, tz, bx, by, bz);
        
        // Generate random numbers (thread-safe)
        fl xi1 = generate_rand_sample();
        fl xi2 = generate_rand_sample();
        
        // Sample direction in tangent space
        fl lx, ly, lz;
        if (sampler_type == "cosine") {
            sample_cosine_hemisphere(xi1, xi2, lx, ly, lz);
        } else {
            sample_uniform_hemisphere(xi1, xi2, lx, ly, lz);
        }
        
        // Transform to world space
        fl wx, wy, wz;
        local_to_world(lx, ly, lz, tx, ty, tz, bx, by, bz, nx, ny, nz, wx, wy, wz);
        
        // Sample environment map
        f_batch wx_b = xs::broadcast(wx);
        f_batch wy_b = xs::broadcast(wy);
        f_batch wz_b = xs:: broadcast(wz);
        f_batch env_r, env_g, env_b;
        sample_environment_map(scene, wx_b, wy_b, wz_b, env_r, env_g, env_b);
        
        // Apply PDF weighting (as shown above)
        // ...
    }
    
    accumulate_color(local_color_buffer, local_mask, 
                    xs::load_aligned(env_contrib_r), 
                    xs::load_aligned(env_contrib_g), 
                    xs::load_aligned(env_contrib_b));
}
```

## Part 6: Performance Measurements and Conclusion
Below are the performance measurements for this assignment. All the measurements are conducted on AMD Ryzen 5600h 6 core processor.

| Scene Preview | Scene Name | Preprocessing Time | Render Time |
| :---: | :---: | :---: | :---: |
| <img width="600" height="600" alt="cube_directional" src="https://github.com/user-attachments/assets/293dbe86-a8c1-4c76-ac00-5a6efc8282da" /> | `cube_directional` | 0 ms | 142 ms |
| <img width="600" height="600" alt="glass_sphere_env" src="https://github.com/user-attachments/assets/739f166f-18d3-4e00-be1a-704b5c796661" /> | `glass_sphere_env` | 193 ms | 352 ms|
| <img width="600" height="600" alt="head_env_light" src="https://github.com/user-attachments/assets/f1579754-9082-4ab0-9ed9-98aa1a5f1c35" /> | `head_env_light` | 466 ms | 3.99 s |
| <img width="600" height="600" alt="sphere_env_light" src="https://github.com/user-attachments/assets/fea7e751-95ef-4c2b-b993-bd338b847642" /> | `sphere_env_light` | 34 ms | 32.97 s |
| <img width="600" height="600" alt="audi-tt-glacier" src="https://github.com/user-attachments/assets/b4bcab69-a23b-410f-930b-a4acbc98b092" /> | `audi-tt-glacier` | 774 ms | 110.4 s |
| <img width="600" height="600" alt="audi-tt-pisa" src="https://github.com/user-attachments/assets/3366ce30-e9d6-4088-83d5-769169b18451"  /> | `audi-tt-pisa` | 1089 ms | 110.86 s |
