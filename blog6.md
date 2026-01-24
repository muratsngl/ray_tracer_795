# ray_tracer::devlog_6 - Advanced Path Tracing: NEE, MIS, Russian Roulette & Splitting

After implementing HDR rendering and environment lighting, I embarked on the most challenging upgrade yet: transforming my Whitted-style ray tracer into a full-fledged Monte Carlo path tracer with variance reduction techniques. This devlog covers my implementation of Next Event Estimation (NEE), Multiple Importance Sampling (MIS), Russian Roulette termination, and Adaptive Splitting.

---

## Part 1: Next Event Estimation (NEE) - Explicit Light Sampling

### The Problem: Naive Path Tracing is Noisy

In basic path tracing, I relied on random bounces to eventually hit a light source. For small or distant lights, this is incredibly inefficient—most rays miss the light entirely, producing black pixels (fireflies) that only slowly converge with more samples.

### The Solution: Sample Lights Directly

NEE (also called explicit light sampling) splits the rendering equation into **direct** and **indirect** components. At every surface hit, I explicitly cast a shadow ray toward a randomly sampled point on each light source.

**Mathematical Foundation:**

The rendering equation:
```
L_o(p, ω_o) = L_e(p, ω_o) + ∫_Ω f_r(p, ω_i, ω_o) L_i(p, ω_i) cos(θ_i) dω_i
```

For NEE, I change the integration variable from solid angle (ω) to area (A) on the light source using the Jacobian transformation:

```
dω = (cos(θ_light) / r²) dA
```

This gives the Monte Carlo estimator for direct lighting:

```
L_direct ≈ f_r · L_source · cos(θ_i) · cos(θ_light) · / (r² · p(A))
```

where:
- `r²` is the squared distance to the light
- `p(A) = 1/Area_light` is the probability density of sampling that point
- `cos(θ_light)` is the angle between the light normal and the incoming ray (absolute value used per Equation 18 in PBRT)

### Implementation: Sampling Light Spheres and Light Meshes

Here's my implementation for sampling spherical area lights with proper PDF conversion:

```cpp
// Sample LightSpheres with NEE
for (size_t light_idx = 0; light_idx < scene.light_spheres.size(); light_idx++) {
    const LightSphere& light = scene.light_spheres[light_idx];
    
    // Sample a random point on the sphere surface
    Vec3f center = {scene.vertex_data__.v_pos_x[light.center_vertex_id],
                   scene.vertex_data__.v_pos_y[light.center_vertex_id],
                   scene.vertex_data__.v_pos_z[light.center_vertex_id]};
    
    if (light.transform.m11 != 0.0f) {
        center = transform_point(center, light.transform);
    }
    
    for (int lane = 0; lane < 8; lane++) {
        if (!active_mask.get(lane)) continue;
        
        // Generate uniform samples on sphere surface
        fl xi1 = generate_rand_sample();
        fl xi2 = generate_rand_sample();
        
        Vec3f sampled_point, sphere_normal;
        sample_sphere_point(xi1, xi2, light.radius, center, sampled_point, sphere_normal);
        
        sampled_x[lane] = sampled_point.x;
        sampled_y[lane] = sampled_point.y;
        sampled_z[lane] = sampled_point.z;
        light_normal_x[lane] = sphere_normal.x;
        light_normal_y[lane] = sphere_normal.y;
        light_normal_z[lane] = sphere_normal.z;
    }
    
    // ... (continued below with visibility and PDF calculation)
}
```

**Uniform Sphere Sampling:**

```cpp
inline void sample_sphere_point(fl xi1, fl xi2, fl radius,
                                const Vec3f& center,
                                Vec3f& sampled_point, Vec3f& sphere_normal) {
    // Uniform sphere sampling
    fl phi = 2.0f * M_PI * xi1;
    fl cos_theta = 1.0f - 2.0f * xi2;  // Maps [0,1] to [-1,1]
    fl sin_theta = std::sqrt(1.0f - cos_theta * cos_theta);
    
    sphere_normal.x = sin_theta * std::cos(phi);
    sphere_normal.y = sin_theta * std::sin(phi);
    sphere_normal.z = cos_theta;
    
    sampled_point.x = center.x + sphere_normal.x * radius;
    sampled_point.y = center.y + sphere_normal.y * radius;
    sampled_point.z = center.z + sphere_normal.z * radius;
}
```

**PDF Conversion from Area to Solid Angle:**

```cpp
// Calculate cos(θ_light) = -wi · n_light (can be negative if back-facing)
f_batch cos_theta_light_raw = -(wi_x * light_norm_x + wi_y * light_norm_y + wi_z * light_norm_z);

// For geometry check: reject back-facing lights
b_batch light_facing = (cos_theta_light_raw > xs::broadcast(0.0f));

// For PDF calculation: use absolute value (per Equation 18 in PBRT)
f_batch cos_theta_light = xs::max(xs::abs(cos_theta_light_raw), xs::broadcast(1e-6f));

// Shadow test
b_batch is_occluded;
intersect_tlas_any_hit_wrapper(shadow_rays, scene, distance - epsilon, is_occluded, valid_geometry);
b_batch in_light = valid_geometry & (!is_occluded);

// PDF conversion: p(A) = 1 / (4π r²) for sphere
fl sphere_area = 4.0f * M_PI * light.radius * light.radius;
f_batch pdf_area = xs::broadcast(1.0f / sphere_area);
f_batch pdf_light = pdf_area * distance_sqr / cos_theta_light;  // Convert to solid angle measure
```

### Light Mesh Sampling with Pre-Built CDFs

For mesh lights (arbitrary triangulated geometries), I needed a way to sample proportional to triangle area. I pre-built Cumulative Distribution Functions (CDFs) during parsing to avoid runtime overhead:

```cpp
// Sample LightMeshes with area-weighted triangle selection
for (size_t light_idx = 0; light_idx < scene.light_meshes.size(); light_idx++) {
    const LightMesh& light = scene.light_meshes[light_idx];
    
    // Use pre-built CDF (computed in parser for efficiency)
    if (!light.cdf || light.cdf->triangle_areas.empty()) continue;
    const LightMeshCDF& cdf = *light.cdf;
    
    for (int lane = 0; lane < 8; lane++) {
        if (!active_mask.get(lane)) continue;
        
        // Sample triangle based on area using CDF
        fl xi_tri = generate_rand_sample();
        int tri_idx = cdf.sample_triangle(xi_tri);
        
        // Get triangle vertices
        int idx0 = light.face_indices[tri_idx * 3 + 0] + light.vertex_offset;
        int idx1 = light.face_indices[tri_idx * 3 + 1] + light.vertex_offset;
        int idx2 = light.face_indices[tri_idx * 3 + 2] + light.vertex_offset;
        
        Vec3f v0 = {scene.vertex_data__.v_pos_x[idx0], scene.vertex_data__.v_pos_y[idx0], scene.vertex_data__.v_pos_z[idx0]};
        Vec3f v1 = {scene.vertex_data__.v_pos_x[idx1], scene.vertex_data__.v_pos_y[idx1], scene.vertex_data__.v_pos_z[idx1]};
        Vec3f v2 = {scene.vertex_data__.v_pos_x[idx2], scene.vertex_data__.v_pos_y[idx2], scene.vertex_data__.v_pos_z[idx2]};
        
        // Apply transform if present
        if (light.transform.m11 != 0.0f) {
            v0 = shader_transform_point(v0, light.transform);
            v1 = shader_transform_point(v1, light.transform);
            v2 = shader_transform_point(v2, light.transform);
        }
        
        // Sample point on triangle using barycentric coordinates
        fl xi1 = generate_rand_sample();
        fl xi2 = generate_rand_sample();
        Vec3f sampled_point, triangle_normal;
        sample_triangle_point(xi1, xi2, v0, v1, v2, sampled_point, triangle_normal);
        
        sampled_x[lane] = sampled_point.x;
        sampled_y[lane] = sampled_point.y;
        sampled_z[lane] = sampled_point.z;
        tri_normal_x[lane] = triangle_normal.x;
        tri_normal_y[lane] = triangle_normal.y;
        tri_normal_z[lane] = triangle_normal.z;
        pdf_values[lane] = 1.0f / cdf.total_area;  // Uniform area PDF
    }
    
    // ... (continued with same visibility/PDF logic as spheres)
}
```

**Triangle Sampling:**

```cpp
inline void sample_triangle_point(fl xi1, fl xi2, 
                                   const Vec3f& v0, const Vec3f& v1, const Vec3f& v2,
                                   Vec3f& sampled_point, Vec3f& triangle_normal) {
    // Uniform sampling: u = 1 - sqrt(ξ₁), v = ξ₂ * sqrt(ξ₁)
    fl sqrt_xi1 = std::sqrt(xi1);
    fl u = 1.0f - sqrt_xi1;
    fl v = xi2 * sqrt_xi1;
    fl w = 1.0f - u - v;
    
    sampled_point.x = v0.x * w + v1.x * u + v2.x * v;
    sampled_point.y = v0.y * w + v1.y * u + v2.y * v;
    sampled_point.z = v0.z * w + v1.z * u + v2.z * v;
    
    // Calculate triangle normal
    Vec3f edge1 = v1 - v0;
    Vec3f edge2 = v2 - v0;
    triangle_normal = normalize(cross(edge1, edge2));
}
```

---

## Part 2: Multiple Importance Sampling (MIS) - Combining BRDF and Light Sampling

### The Problem: Which Sampling Strategy is Better?

With NEE implemented, I had two ways to sample lighting:
1. **BRDF Sampling:** Generate directions based on surface properties (cosine-weighted hemisphere)
2. **Light Sampling (NEE):** Sample points directly on lights

Neither strategy is universally better:
- **Glossy surfaces:** BRDF sampling is better (tight specular lobe)
- **Large area lights:** Light sampling is better (high probability of hitting)
- **Small distant lights + diffuse surfaces:** Light sampling wins

Using only one strategy produces high variance (noise) in scenes where the other would excel. Using both naively causes **double counting** (the same light path is added twice).

### The Solution: MIS Weighting

MIS combines both strategies by weighting their contributions based on their probability densities. The **Balance Heuristic** minimizes variance:

```
w_A(x) = p_A(x) / (p_A(x) + p_B(x))
```

where:
- `p_A(x)` is the PDF of strategy A (e.g., BRDF sampling)
- `p_B(x)` is the PDF of strategy B (e.g., Light sampling)

If one PDF is very small (e.g., sampling a tiny light via BRDF), the weight drops to near zero, suppressing the noise spike.

### Implementation: Balance and Power Heuristics

```cpp
// Balance heuristic: w = p / (p + q)
inline fl mis_weight_balance(fl pdf_this, fl pdf_other) {
    fl sum = pdf_this + pdf_other;
    if (sum < 1e-10f) return 0.0f;
    return pdf_this / sum;
}

// Power heuristic: w = p² / (p² + q²)
// Provides slightly better variance reduction in practice
inline fl mis_weight_power(fl pdf_this, fl pdf_other) {
    fl p_sqr = pdf_this * pdf_this;
    fl q_sqr = pdf_other * pdf_other;
    fl sum = p_sqr + q_sqr;
    if (sum < 1e-10f) return 0.0f;
    return p_sqr / sum;
}

// Calculate MIS weight based on camera settings
inline fl calculate_mis_weight(const Camera& cam, fl pdf_this, fl pdf_other) {
    if (cam.mis_heuristic == "power") {
        return mis_weight_power(pdf_this, pdf_other);
    } else if (cam.mis_heuristic == "01" || cam.mis_heuristic == "MIS_01") {
        // Binary MIS: randomly choose one strategy with probability proportional to PDF
        // For deterministic execution, pick strategy with higher PDF
        return 1.0f;  // Already chose this path, so weight = 1
    } else {  // "balance" or default
        return mis_weight_balance(pdf_this, pdf_other);
    }
}
```

### Applying MIS in NEE

When NEE hits a light, I need to calculate what the BRDF sampling PDF *would have been* for that direction:

```cpp
// Calculate BRDF PDF for MIS (what would the PDF be if we sampled via BRDF?)
f_batch pdf_brdf;
if (cam.importance_sampling) {
    // Cosine-weighted: PDF = cos(θ) / π
    pdf_brdf = cos_theta / xs::broadcast(static_cast<fl>(M_PI));
} else {
    // Uniform hemisphere: PDF = 1 / (2π)
    pdf_brdf = xs::broadcast(static_cast<fl>(1.0f / (2.0f * M_PI)));
}

// Calculate MIS weight
alignas(32) fl pdf_light_arr[8], pdf_brdf_arr[8], mis_weight_arr[8];
xs::store_aligned(pdf_light_arr, pdf_light);
xs::store_aligned(pdf_brdf_arr, pdf_brdf);

for (int lane = 0; lane < 8; lane++) {
    if (in_light.get(lane)) {
        mis_weight_arr[lane] = calculate_mis_weight(cam, pdf_light_arr[lane], pdf_brdf_arr[lane]);
    } else {
        mis_weight_arr[lane] = 0.0f;
    }
}
f_batch mis_weight = xs::load_aligned(mis_weight_arr);

// Monte Carlo estimator with MIS: MIS_weight * L_e * BRDF * cos(θ) / pdf_light(ω)
f_batch contrib_r = mis_weight * radiance_r * brdf_r * cos_theta / pdf_light;
f_batch contrib_g = mis_weight * radiance_g * brdf_g * cos_theta / pdf_light;
f_batch contrib_b = mis_weight * radiance_b * brdf_b * cos_theta / pdf_light;
```

### BRDF Path Hitting Lights (Indirect Hits)

When a BRDF-sampled ray accidentally hits a light, I apply the reverse logic:

```cpp
// Calculate light PDF with proper Jacobian transformation
// pdf_direction = pdf_area * (distance² / |cos(θ_light)|)
fl pdf_light = 0.0f;
int emissive_idx = ray_pack.emissive_id[i];
fl distance = ray_pack.t_min[i];
fl distance_sqr = distance * distance;

// Get light normal
fl light_norm_x = hit_norm_x.get(i);
fl light_norm_y = hit_norm_y.get(i);
fl light_norm_z = hit_norm_z.get(i);

// Calculate cos(θ_light) = |dot(light_normal, -ray_direction)|
fl cos_theta_light = std::abs(light_norm_x * (-wi_x) + 
                              light_norm_y * (-wi_y) + 
                              light_norm_z * (-wi_z));
cos_theta_light = std::max(cos_theta_light, 1e-6f);  // Avoid division by zero

if (emissive_type == 0) {  // LightSphere
    const auto& light_sphere = scene.light_spheres[emissive_idx];
    fl sphere_area = 4.0f * M_PI * light_sphere.radius * light_sphere.radius;
    pdf_light = (1.0f / sphere_area) * (distance_sqr / cos_theta_light);
} else if (emissive_type == 1) {  // LightMesh
    const auto& light_mesh = scene.light_meshes[emissive_idx];
    if (light_mesh.cdf && light_mesh.cdf->total_area > 0.0f) {
        pdf_light = (1.0f / light_mesh.cdf->total_area) * (distance_sqr / cos_theta_light);
    }
}

// Calculate MIS weight
fl mis_weight = calculate_mis_weight(cam, pdf_brdf, pdf_light);

// Add MIS-weighted radiance
fl contrib_r = hit_radiance_r * tp_r[i] * mis_weight;
fl contrib_g = hit_radiance_g * tp_g[i] * mis_weight;
fl contrib_b = hit_radiance_b * tp_b[i] * mis_weight;
```



---

## Part 3: Russian Roulette - Unbiased Path Termination

### The Problem: When to Stop Bouncing?

Naive approaches:
1. **Fixed depth cutoff:** Biased (ignores energy beyond max depth)
2. **Throughput threshold:** Also biased (ignores low but non-zero contributions)

I needed an **unbiased** termination strategy.

### The Solution: Probabilistic Termination with Compensation

Russian Roulette randomly terminates rays with probability `1 - q`, but boosts survivors by `1/q` to conserve energy on average.

**Mathematical Proof:**

Given a contribution `F` and survival probability `q`:

```
E[F] = q · (F/q) + (1-q) · 0 = F
```

The expected value is unchanged—the estimator remains **unbiased**.

### Implementation

```cpp
// Russian Roulette path termination (after minimum depth)
bool path_continues = true;
fl rr_boost = 1.0f;

if (cam.russian_roulette && ray_pack.depth[i] >= cam.min_recursion_depth) {
    // Calculate survival probability based on current throughput
    // Use max channel to avoid killing paths with colored light
    fl max_throughput = std::max(tp_r[i], std::max(tp_g[i], tp_b[i]));
    
    // Clamp to 0.99 to prevent bias (survival probability must never exceed 1.0)
    // This is CRITICAL for unbiasedness
    fl q = std::min(max_throughput, 0.99f);
    
    // Generate random number for Russian Roulette
    fl rr_sample = generate_rand_sample();
    
    if (rr_sample > q) {
        // Ray dies (terminate path)
        path_continues = false;
    } else {
        // Ray survives - boost throughput to compensate for termination probability
        // This ensures the estimator remains unbiased
        rr_boost = 1.0f / q;
    }
}

if (!path_continues) {
    // Path terminated by Russian Roulette
    next_active[i] = false;
    local_shade[i] = false;
} else {
    // Path continues - apply rr_boost to throughput later
    // ...
}
```
I clamp `q` to 0.99 to ensure `1/q` doesn't explode. Even very bright paths (throughput near 1.0) must have a small termination probability to avoid infinite bounces.

### Applying the Boost

The boost is applied to the **throughput** of surviving rays:

```cpp
// Apply Russian Roulette boost (already calculated before splitting)
split_tp_r *= rr_boost;
split_tp_g *= rr_boost;
split_tp_b *= rr_boost;
```

---














## Part 4: Results 

### Performance
Render Times for all of the 
Cornellbox_Diffuse scenes was:189125ms
Cornellbox_Glass_mirror was:427745ms
| Scene Preview | Scene Name | Preprocessing Time | Render Time |
| :---: | :---: | :---: | :---: |
|<img width="512" height="512" alt="cornell_box_importance_nee_mis_balance_phot" src="https://github.com/user-attachments/assets/458fe14c-5250-434d-b4df-f10de59453f5" />| `cornell_box_importance_nee_mis_balance_phot` | - | - |
|<img width="512" height="512" alt="cornell_box_importance_nee_mis_balance_russian_phot" src="https://github.com/user-attachments/assets/e8320685-82f8-4865-b44b-bb344297753e" /> | `cornell_box_importance_nee_mis_balance_russian_phot` | - | -|
|<img width="512" height="512" alt="cornellbox_prism_light_phot" src="https://github.com/user-attachments/assets/28267b61-9220-4fc8-ae0b-7af4be85e151" />| `cornellbox_prism_light_phot` |0|35412|
|<img width="512" height="512" alt="cornellbox_sphere_light_phot" src="https://github.com/user-attachments/assets/43412e44-3e97-4de8-8957-05e05a487206" />|`cornellbox_sphere_light_phot`|0|31439|
|<img width="512" height="512" alt="diffuse_cornell_box_importance_nee_mis_balance_russian_1600_phot" src="https://github.com/user-attachments/assets/180d177f-2c78-4add-94dd-0d5047b0c58f" />| `diffuse_cornell_box_importance_nee_mis_balance_russian_1600_phot` | - | - |
## Part 5: Missing Features & Scenes 
<img width="512" height="512" alt="cornellbox_prism_light_phot_weird_graze" src="https://github.com/user-attachments/assets/b04c1e08-61c2-4f2d-a297-ae488e56758a" />
<img width="512" height="512" alt="cornellbox_sphere_light_phot_weird_graze" src="https://github.com/user-attachments/assets/e27cb291-72b6-47de-8e4b-f2f0f1834a49" />
#### Inaccurate Shadows

Solved by changing the shadow ray offsetting from ray direction*epsilon to surface_normal * epsilon. FIXED

#### Splitting Correlation Bug

<img width="512" height="512" alt="image" src="https://github.com/user-attachments/assets/7b107157-ad1e-40d1-bb71-91f46683f204" />

There was a source of correlation which I could not find even though I have attempted numerous times. As you can see the image has a tiled appereance. Those tiles are the ray packets I use for path tracing. NOT FIXED


