# ray_tracer::devlog_7 - Volumetric Rendering: Clouds, Smoke, and the Physics of Light in Participating Media

This devlog covers Fatih & Murat's journey from integrating NanoVDB for efficient volume data access to implementing physically-based scattering with the Henyey-Greenstein phase function.


## Part 1: NanoVDB Integration - CPU-Based Volume Data Access

### The Volume Struct: First-Class Scene Objects

We wanted volumes to be treated as primitives, just like spheres and meshes. Our first step was adding a Volume struct to the scene definition:

```cpp
struct Volume {
    int id;
    std::string vdb_path;           // Path to .nvdb file
    nanovdb::GridHandle<> handle;   // NanoVDB grid handle
    nanovdb::FloatGrid* grid;       // Pointer to the grid
    Vec3f absorption;               // σ_a: absorption coefficient
    Vec3f scattering;               // σ_s: scattering coefficient  
    fl density_multiplier;          // Global density scale
    fl phase_g;                     // Henyey-Greenstein anisotropy [-1, 1]
    Mat4x4 transform;               // World-to-volume transform
    Mat4x4 inv_transform;           // Volume-to-world transform
    BoundingBox aabb;               // Axis-aligned bounding box in world space
};
```
Storing both `absorption` and `scattering` separately lets us balance Beer's Law attenuation (absorption) against in-scattering contribution (scattering) independently—critical for getting diverse lighting outputs.

### JSON Parsing: Handling OpenVDB Paths

We extended our JSON parser to recognize a new `"Volumes"` array:

```cpp
if (doc.HasMember("Volumes") && doc["Volumes"].IsArray()) {
    for (auto& vol_obj : doc["Volumes"].GetArray()) {
        Volume vol;
        vol.id = vol_obj["Id"].GetInt();
        vol.vdb_path = vol_obj["VdbPath"].GetString();
        
        // Parse extinction coefficients
        const auto& abs = vol_obj["Absorption"].GetArray();
        vol.absorption = {abs[0].GetFloat(), abs[1].GetFloat(), abs[2].GetFloat()};
        
        const auto& scat = vol_obj["Scattering"].GetArray();
        vol.scattering = {scat[0].GetFloat(), scat[1].GetFloat(), scat[2].GetFloat()};
        
        vol.density_multiplier = vol_obj["DensityMultiplier"].GetFloat();
        vol.phase_g = vol_obj["PhaseG"].GetFloat();
        
        // Parse transformation matrix (4x4)
        const auto& transf = vol_obj["Transformation"].GetArray();
        parse_mat4x4(transf, vol.transform);
        vol.inv_transform = inverse(vol.transform);
        
        scene.volumes.push_back(vol);
    }
}
```

### On-the-Fly NanoVDB Loading

During scene preprocessing, we load the `.nvdb` files and extract grid pointers:

```cpp
for (Volume& vol : scene.volumes) {
    vol.handle = nanovdb::io::readGrid(vol.vdb_path);
    vol.grid = vol.handle.grid<float>();
    
    if (!vol.grid) {
        throw std::runtime_error("Failed to load NanoVDB grid: " + vol.vdb_path);
    }
    
    // Compute world-space AABB from index space bounds
    auto index_bbox = vol.grid->indexBBox();
    Vec3f corners[8] = { /* transform all 8 corners */ };
    vol.aabb = compute_world_aabb(corners);
}
```

Why NanoVDB? It's compact and fast.

---

## Part 2: Ray Marching Implementation - Fixed Step-Size Integration

### The Rendering Equation for Volumes

The volume rendering equation describes how light changes as it travels through a participating medium:

```text
L_o(x, ω) = ∫₀ᵗ T(x, x_t) · [σ_s(x_t) · L_i(x_t, ω) + σ_a(x_t) · 0] dt
```

where:
- **T(x, x_t)** is the transmittance (Beer's Law exponential extinction)
- **σ_s** is the scattering coefficient (how much light scatters *into* the ray)
- **σ_a** is the absorption coefficient (how much light is absorbed)
- **L_i** is the incoming radiance at point x_t (from lights or environment)

The total extinction coefficient is **σ_t = σ_a + σ_s**.

### World-to-Index Space Mapping

Before marching, we transform rays from world space to the volume's index space:

```cpp
void transform_ray_to_index_space(const Vec3f& world_origin, const Vec3f& world_dir,
                                  const Volume& vol,
                                  Vec3f& index_origin, Vec3f& index_dir) {
    // Apply inverse transformation matrix
    index_origin = transform_point(world_origin, vol.inv_transform);
    index_dir = transform_direction(world_dir, vol.inv_transform);
}
```

This lets us directly query the NanoVDB grid using voxel coordinates.

### AABB Optimization: Clipping Ray Intervals

We only march through the volume's bounding box:

```cpp
bool intersect_aabb(const Vec3f& ray_o, const Vec3f& ray_d,
                    const BoundingBox& aabb,
                    fl& t_near, fl& t_far) {
    Vec3f inv_d = {1.0f / ray_d.x, 1.0f / ray_d.y, 1.0f / ray_d.z};
    
    fl t1 = (aabb.min.x - ray_o.x) * inv_d.x;
    fl t2 = (aabb.max.x - ray_o.x) * inv_d.x;
    fl t3 = (aabb.min.y - ray_o.y) * inv_d.y;
    fl t4 = (aabb.max.y - ray_o.y) * inv_d.y;
    fl t5 = (aabb.min.z - ray_o.z) * inv_d.z;
    fl t6 = (aabb.max.z - ray_o.z) * inv_d.z;
    
    t_near = std::max({std::min(t1, t2), std::min(t3, t4), std::min(t5, t6)});
    t_far = std::min({std::max(t1, t2), std::max(t3, t4), std::max(t5, t6)});
    
    return t_far >= t_near && t_far >= 0.0f;
}
```

This clips the marching interval [t_near, t_far] to only the volume's extent—no wasted samples outside.

### Fixed Step-Size Integration Loop

Here's our ray marching kernel for a single volume:

```cpp
Vec3f march_volume(const Vec3f& ray_o, const Vec3f& ray_d,
                   fl t_near, fl t_far,
                   const Volume& vol,
                   const Scene& scene) {
    const fl step_size = 0.1f;  // Fixed step (tunable)
    const int max_steps = 512;
    
    Vec3f accumulated_light = {0.0f, 0.0f, 0.0f};
    Vec3f transmittance = {1.0f, 1.0f, 1.0f};  // Start fully transparent
    
    fl t = t_near;
    int steps = 0;
    
    while (t < t_far && steps < max_steps) {
        Vec3f sample_pos_world = ray_o + ray_d * t;
        Vec3f sample_pos_index = transform_point(sample_pos_world, vol.inv_transform);
        
        // Sample density from NanoVDB grid
        fl density = sample_grid(vol.grid, sample_pos_index) * vol.density_multiplier;
        
        if (density > 1e-6f) {
            // Compute extinction
            Vec3f sigma_t = (vol.absorption + vol.scattering) * density;
            
            // Calculate transmittance decrement (Beer's Law)
            Vec3f tau = sigma_t * step_size;
            Vec3f step_transmittance = {exp(-tau.x), exp(-tau.y), exp(-tau.z)};
            
            // Accumulate in-scattering light
            Vec3f in_scatter = compute_inscattering(sample_pos_world, ray_d, vol, scene);
            Vec3f sigma_s_rgb = vol.scattering * density;
            Vec3f contribution = transmittance * (1.0f - step_transmittance) * sigma_s_rgb * in_scatter;
            accumulated_light += contribution;
            
            // Update transmittance
            transmittance *= step_transmittance;
            
            // Early termination if transmittance too low
            if (transmittance.x < 0.01f && transmittance.y < 0.01f && transmittance.z < 0.01f) {
                break;
            }
        }
        
        t += step_size;
        steps++;
    }
    
    return accumulated_light;
}
```

**Simplification:** We use a fixed step size instead of adaptive sampling. This trades some accuracy for predictable performance and simpler code.

---

## Part 3: Volumetric Shadows - Transmittance to Light Sources

### Computing In-Scattering

At each sample point, we need to know how much light reaches it from each light source:

```cpp
Vec3f compute_inscattering(const Vec3f& sample_pos, const Vec3f& view_dir,
                           const Volume& vol, const Scene& scene) {
    Vec3f total_light = {0.0f, 0.0f, 0.0f};
    
    for (const PointLight& light : scene.point_light_data__) {
        Vec3f to_light = light.position - sample_pos;
        fl distance = length(to_light);
        to_light /= distance;  // Normalize
        
        // Compute transmittance from sample point to light (volumetric shadow)
        Vec3f light_transmittance = compute_transmittance(sample_pos, to_light, distance, vol);
        
        // Apply inverse square attenuation
        Vec3f light_intensity = light.radiance / (distance * distance);
        
        // Apply phase function for anisotropic scattering
        fl phase = henyey_greenstein(-view_dir, to_light, vol.phase_g);
        
        total_light += light_transmittance * light_intensity * phase;
    }
    
    return total_light;
}
```

### Beer's Law Application: Exponential Extinction

To compute transmittance from the sample point to the light, we march *again* (nested marching):

```cpp
Vec3f compute_transmittance(const Vec3f& from, const Vec3f& direction,
                            fl max_distance, const Volume& vol) {
    const fl shadow_step = 0.2f;  // Coarser steps for performance
    Vec3f transmittance = {1.0f, 1.0f, 1.0f};
    
    fl t = 0.0f;
    while (t < max_distance) {
        Vec3f pos_world = from + direction * t;
        Vec3f pos_index = transform_point(pos_world, vol.inv_transform);
        
        fl density = sample_grid(vol.grid, pos_index) * vol.density_multiplier;
        
        if (density > 1e-6f) {
            Vec3f sigma_t = (vol.absorption + vol.scattering) * density;
            Vec3f tau = sigma_t * shadow_step;
            transmittance *= Vec3f{exp(-tau.x), exp(-tau.y), exp(-tau.z)};
            
            // Early exit if nearly opaque
            if (transmittance.x < 0.01f && transmittance.y < 0.01f && transmittance.z < 0.01f) {
                return {0.0f, 0.0f, 0.0f};
            }
        }
        
        t += shadow_step;
    }
    
    return transmittance;
}
```

**Beer's Law:**

```text
T(a → b) = exp(-∫ σ_t(x) dx)
```

For a homogeneous segment:

```text
T = exp(-σ_t · Δt)
```



---

## Part 4: Henyey-Greenstein Phase Function - Anisotropic Scattering

### The Math Behind Phase Functions

The phase function p(θ) describes how light scatters at a given angle θ between the incoming and outgoing directions. For isotropic scattering (uniform in all directions):

```text
p(θ) = 1 / (4π)
```

But clouds aren't isotropic! They exhibit **forward scattering** (light prefers to continue in the same direction). The Henyey-Greenstein phase function models this:

```text
p(θ, g) = (1 - g²) / (4π · (1 + g² - 2g·cos(θ))^(3/2))
```

where:
- **g ∈ [-1, 1]** is the anisotropy parameter
- g = 0: isotropic
- g > 0: forward scattering (typical for clouds, g ≈ 0.6 to 0.8)
- g < 0: backward scattering

### Implementation

```cpp
fl henyey_greenstein(const Vec3f& view_dir, const Vec3f& light_dir, fl g) {
    fl cos_theta = dot(view_dir, light_dir);
    fl g2 = g * g;
    fl denom = 1.0f + g2 - 2.0f * g * cos_theta;
    
    if (denom < 1e-6f) return 1.0f / (4.0f * M_PI);  // Avoid division by zero
    
    fl numerator = 1.0f - g2;
    fl denominator = 4.0f * M_PI * pow(denom, 1.5f);
    
    return numerator / denominator;
}
```



---

## Part 5: Shortcuts

### The "Shadow Density Multiplier" Trick


```cpp
// In compute_transmittance (shadow ray):
fl shadow_density = density * 0.5f;  // Reduce density for shadow rays
Vec3f sigma_t = (vol.absorption + vol.scattering) * shadow_density;
```

This is physically incorrect but visually pleasing. It simulates multiple scattering (where light bounces around inside the cloud before exiting) without the computational cost.

This creates detailed clouds by adding variance to the light transport.

### From Scalar to Vector: Colored Shadows

Initially, we computed transmittance as a scalar:

```cpp
fl transmittance = exp(-sigma_t * distance);
```

But colored volumes (like fire or stained glass fog) need **RGB transmittance**:

```cpp
Vec3f transmittance = {exp(-sigma_t.x * distance),
                       exp(-sigma_t.y * distance),
                       exp(-sigma_t.z * distance)};
```
This creates realistic colored fog effects.

---



## Part 6: Integration Wiring - Blending Volumes with Surface Geometry

### The Hybrid Rendering Loop

Volumes and surfaces coexist in our scenes, so we needed to decide the rendering order:

```cpp
for (int lane = 0; lane < 8; lane++) {
    if (!active_mask[lane]) continue;
    
    Vec3f ray_o = {ray_pack.o_x[lane], ray_pack.o_y[lane], ray_pack.o_z[lane]};
    Vec3f ray_d = {ray_pack.d_x[lane], ray_pack.d_y[lane], ray_pack.d_z[lane]};
    
    fl t_surface = ray_pack.t_min[lane];  // Distance to nearest surface
    Vec3f volume_radiance = {0.0f, 0.0f, 0.0f};
    Vec3f transmittance = {1.0f, 1.0f, 1.0f};
    
    // 1. March through ALL volumes up to the surface
    for (const Volume& vol : scene.volumes) {
        fl t_near, t_far;
        if (intersect_aabb(ray_o, ray_d, vol.aabb, t_near, t_far)) {
            t_near = std::max(t_near, 0.0f);
            t_far = std::min(t_far, t_surface);  // Clip to surface
            
            if (t_far > t_near) {
                Vec3f vol_contrib = march_volume(ray_o, ray_d, t_near, t_far, vol, scene);
                Vec3f vol_transmittance = compute_transmittance(ray_o, ray_d, t_far - t_near, vol);
                
                volume_radiance += transmittance * vol_contrib;
                transmittance *= vol_transmittance;
            }
        }
    }
    
    // 2. Add attenuated surface contribution
    Vec3f surface_radiance = {color_block_float.r[lane],
                              color_block_float.g[lane],
                              color_block_float.b[lane]};
    
    Vec3f final_radiance = volume_radiance + transmittance * surface_radiance;
    
    color_block_float.r[lane] = final_radiance.x;
    color_block_float.g[lane] = final_radiance.y;
    color_block_float.b[lane] = final_radiance.z;
}
```

**Key Insight:** Volumes are composited *in front of* surfaces using the over operator:

```text
C_final = C_volume + T_volume · C_surface
```

where T_volume is the transmittance through all volumes between the camera and the surface.

---

