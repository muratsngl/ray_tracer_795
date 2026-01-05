
# Ray_tracer::devlog_4

## Part 1: Background Textures
<img width="597" height="847" alt="image" src="https://github.com/user-attachments/assets/a2e07c83-a5da-4547-b1cf-57a34c67558e" />
<img width="607" height="116" alt="image" src="https://github.com/user-attachments/assets/84fea210-d385-46d8-97e4-dad710e40cbd" />
<img width="653" height="734" alt="image" src="https://github.com/user-attachments/assets/af1371e7-e206-4732-8611-a7dcabc8699f" />

The first picture here shows my main image sampling function. If a ray didn't intersect any geometry I would simply call sample_background function(3rd picture). Which by transforming the current ray direction vector into u v coordinates and passing the coordinates as input to my main sampling function, samples from the image and assings the pixel color to the background texture 

## Part 2: Image Based Color Sampling (diffuse/specular maps)

For this part I again used my main sampling function, in the intersection part I have calculated the UV coordinates of all geometries then passed them to the main sampling function in order to obtain the texture rgb value. Later on based on the texture type(which is again bundled with the ray hitrecords) I assign the reflectance values dictated by the textures. Below is the code part that assigns the values to the textures based on the decal mode.

<img width="252" height="319" alt="image" src="https://github.com/user-attachments/assets/b0289072-d108-40b6-bb0b-add4dfb35adb" />

## Part 3: Normal Maps

For normal maps I have done all my calculations inside the intersection routines, while in the shading part I interpreted every normal the same thus reducing the need to carry additonal information regarding rays. For this inside the intersection tests I check whether or not there exists a normal map texture attached to the object. If there is any I use the below code to sample from it. 
<img width="712" height="679" alt="image" src="https://github.com/user-attachments/assets/af185b8e-fbe5-42e1-9d55-e34c57de9a9d" />

To properly implement normal mapping, I needed a way to calculate the tangent and bitangent vectors for my triangles. Since my renderer uses batch processing for performance (using SIMD types like f_batch), I wrote a specialized function to compute the TBN (Tangent, Bitangent, Normal) basis for a batch of intersections simultaneously.
The code first calculates the determinant of the UV coordinates (Δu1​Δv2​−Δu2​Δv1​). This is crucial for mapping the 2D texture space coordinates back into 3D world space. I included an epsilon check here to handle degenerate triangles or cases where UV mapping might result in a zero determinant, preventing division-by-zero errors.
Once I have the initial tangent and bitangent vectors using the inverse determinant, I perform a Gram-Schmidt orthogonalization. This ensures that the tangent is perfectly perpendicular to the surface normal by subtracting the projection of the tangent onto the normal:
Tortho​=T−N(N⋅T)
After normalizing the tangent to unit length, I calculate the final bitangent using a cross product (N×T). This approach is more robust than using the raw bitangent from the UV calculation because it guarantees a perfectly orthonormal basis.

<img width="595" height="616" alt="image" src="https://github.com/user-attachments/assets/86886e6b-1637-463d-bce0-afd3dd06df7b" />

After calculating the basis and applying the perturbation. I select among the valid rays that the normal texture exists for. 

## Part 4: Bump Maps

Again bump mapping is applied at the intersection level for the same reasons as the normal mapping. In bump mapping I again calculated the same tbn basis then applied the perturbations with the below code.

<img width="806" height="920" alt="image" src="https://github.com/user-attachments/assets/4e37180a-77be-460d-a7c8-09a1c49bce32" />


## Part 5: Perlin Noise Textures
For the perling noise calculations I have used the original ken perlin implementation with the 6t^5 - 15t^4 + 10t^3 fading function. In addition to that turbulence is also added. You can see the implementation below. 

<img width="806" height="923" alt="Screenshot from 2025-12-24 21-50-33" src="https://github.com/user-attachments/assets/729fbb1f-4697-4778-8756-2165bc15f9b3" />













## Part 6: Results
Below are the performance measurements for this assignment. All the measurements are conducted on AMD Ryzen 5600h 6 core processor.
<video src="https://github.com/user-attachments/assets/b87d809f-d637-4c83-a3fd-55fee0d36187" autoplay muted loop playsinline width="100%"></video>

| Scene Preview | Scene Name | Preprocessing Time | Render Time |
| :---: | :---: | :---: | :---: |
|<img width="800" height="800" alt="image" src="https://github.com/user-attachments/assets/491046ba-35f9-4f55-93c2-5b26209403c9" /> | `brickwall with normal map` | 17 ms | 130 ms |
|<img width="800" height="800" alt="image" src="https://github.com/user-attachments/assets/46753c3a-e2e6-465d-9824-1d4882cb2baf" />| `bump_mapping_transformed` | 57 ms | 124 ms |
|<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/741c505d-a48c-433e-8f32-913fb44612c0" />| `my_tap_final` | 77 ms |77.5 s |
|<img width="800" height="800" alt="image" src="https://github.com/user-attachments/assets/c075d0f1-77f1-4681-a5e5-0c463346e3b4" />| `bump_mapping_transformed` | 0 ms | 73 ms |
|<img width="800" height="800" alt="image" src="https://github.com/user-attachments/assets/dfde246a-c9da-4e2c-a606-58edf1f1cba3" />| `tunnel_of_doom_010` | 0 ms | 690 ms |
|<img width="1280" height="720" alt="image" src="https://github.com/user-attachments/assets/d0d0f21f-2af8-48a5-9904-bf0dd401c44c" />| `Veach Ajar` | 322 ms | 2192 ms |

