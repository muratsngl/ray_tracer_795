# Ray_tracer::devlog_3
## Part 1: Compulsory Migration
<img width="1356" height="306" alt="image" src="https://github.com/user-attachments/assets/69808650-1d73-4973-ac0f-60b1dc71cd33" />
As blog.metu.edu.tr servers are down I have decided to move the blogposts to Github from now on.   

## Part 2: Previous Bugs
From the previous assignments I had 2 unsolved issues. One was that with Chinese Dragon scene due to limited floating point precision and the triangles are being too small I was getting incorrect/incomplete renders. The other was due to my refraction and internal absorption logic I was not observing green glass effect in the metal_mirror_planes scene.
<img width="8000" height="8000" alt="chinese_dragon" src="https://github.com/user-attachments/assets/08f8af89-a6a3-42a4-893a-82bbf80cb62a" />
<img width="800" height="480" alt="dragon_dynamic_false_shadows" src="https://github.com/user-attachments/assets/6b114746-a270-415c-8f99-e3e7aaa99116" />

Green glass is missing in the above dragon. Which was a result of backface culling as the refracted rays inside the objects, when intersected while leaving the object fails at the front face test and hence the intersection was not accepted as a result the rays were never accepted as exiting so the attenuation was never applied. I have managed to escape this by culling on the condition that the currently intersected material was not a dielectric.

<img width="408" height="308" alt="image" src="https://github.com/user-attachments/assets/13c76f65-ef83-40f5-9b07-176fb1556a9c" />

For the precision issues I have discussed with my friends regarding to their solutions. Some of them included using doubles and lower default epsilons, while the others suggested scaling the scene altogether then dividing the t value by the scale factor. Both of which are valid solutions. I haven't implemented them in my raytracer as I believe it is more appropriate to fix the scene then to implement scaling for only one scene that requests it. I will be uploading a more forgiving and visually equivalent version of it to this repository for the following attendees of this course. 

## Part 3: Multisampling

For this part I have used Mersenne Twister random number generator that Course Instructor Ahmet Oguz Akyuz suggested. While generating N*N samples for each pixel I have used jittered sampling which divides each pixel to N*N squares each of which holding one sample. Later on results of each sample are averaged to produce the final color output. At this part I have made the error to create the pixel buffer as char before averaging which resulted in overflow, result of which you can see below.

<img width="1314" height="736" alt="image" src="https://github.com/user-attachments/assets/5db2d229-5dcd-4d86-86ed-abce0662f203" />



<img width="1024" height="1024" alt="lobster_ij" src="https://github.com/user-attachments/assets/b69c7076-9e8a-48fa-946d-db38893e9594" />

<img width="1024" height="1024" alt="lobster" src="https://github.com/user-attachments/assets/dbfc002a-60cf-4158-9acc-d58bcd3d7625" />

## Part 4: Aperture and Depth of Field

As you can see in the above code the depth of field effect is achieved by sampling a point on the aperture for each sample and then applying the lens calculations and finding the direction vector of the incoming ray which the lens focuses the to the given position on the aperture.
The effect is really beautiful and adds a lot of character to the rendered images. 


<img width="800" height="800" alt="spheres_dof" src="https://github.com/user-attachments/assets/cfd0f5ca-327d-4767-a469-671ce8960a00" />

<img width="1920" height="1080" alt="focusing_dragons" src="https://github.com/user-attachments/assets/3ce7ec86-6fcd-40fd-93c4-ea5747ec6c00" />


## Part 5: Area Lights and Soft Shadows

In order to cast soft shadows point lights are not sufficient as they are singular light points and their shadows are in a sense binary. Area lights are different in that sense, with the power of random sampling they can provide soft shadows and thus improve realism. The mechanism it achieves this can be summarized as follows. We first sample a random point on the area light by constructing a orthonormal base on the area light plane using the normal and the position of the area light. Later on we calculate the distance and the cosine of the angle between the light-object ray and the area light plane normal to calculate the recieved irradiance at the sample position. Below you can find the code that implements the above functionality.

<img width="629" height="298" alt="image" src="https://github.com/user-attachments/assets/8a2a7cce-bf98-427d-bb94-ee341c3f6562" />

Sampling a random position on the area light for 8 rays in a ray pack.

<img width="847" height="240" alt="image" src="https://github.com/user-attachments/assets/b9c0ed80-5260-486e-805a-0f77d9503f5a" />

Calculating intensity for the area light which will unify the shading pipeline afterwards, meaning after this point we have a ray for which we know the intensity and hence we can make our calculations just as before we did for the point lights. 

<img width="800" height="800" alt="cornellbox_area" src="https://github.com/user-attachments/assets/49690c1c-f228-482b-8002-83107266f133" />

<img width="800" height="800" alt="chessboard_arealight" src="https://github.com/user-attachments/assets/55e5903e-4473-4314-91d6-67ab083b04bc" />

## Part 6: Time Sampling and Motion Blur

This part was the most difficult one for me to understand completely. At first I have completely ignored the transformations for the motion blur enabled objects, then I have sampled the time and then applied the same transform to both the object and the ray at the same time which resulted in no motion blur at all. At last I have forgotten to set the time variable for the shadow rays which are initialized at t=0, resulting in misfitting shadows. In addition to these I have forgotten to update the ray_tracer internal loop that kept track of refracted rays to have time which had no effect on the scenes rendered but I have fixed it nevertheless. You can see the result of the motion blur below.



<img width="750" height="600" alt="cornellbox_boxes_dynamic" src="https://github.com/user-attachments/assets/f1f92219-6077-4f4a-8a95-4865a48870c7" />

<img width="800" height="480" alt="dragon_dynamic" src="https://github.com/user-attachments/assets/afcb22c8-9324-4e95-a182-696e19e77f13" />

Dragon also features the non-vacuum medium attenuation effect.

## Part 7: Roughness and Glossy Reflections

This part was the easiest and most straightforward of all. For the mirror, conductor and dielectric materials, in order to create a brushed/patina look on the surface, we perturb the normals with help of 2 random samples, Orthonormal base of which the samples are taken and roughness parameter which scales the effect.   

<img width="643" height="166" alt="image" src="https://github.com/user-attachments/assets/61284828-967b-46f3-835a-a29acfb9644f" />

<img width="828" height="950" alt="Screenshot from 2025-12-06 17-12-52" src="https://github.com/user-attachments/assets/b9c9dd1f-0913-4f0d-880b-702d547b1803" />

<img width="800" height="800" alt="metal_glass_plates" src="https://github.com/user-attachments/assets/ba14ae25-6df9-4949-8be7-d063530b8ca9" />

<img width="800" height="800" alt="cornellbox_brushed_metal" src="https://github.com/user-attachments/assets/405cf152-4d2b-43a5-8a51-c6a88b3ce736" />


## Part 8: Performance Measurements and Conclusion
https://github.com/user-attachments/assets/6caa2412-1ad2-480a-b62b-953e5ecd4297
| Scene Preview | Scene Name | Preprocessing Time | Render Time |
| :---: | :---: | :---: | :---: |
| <img width="400" height="400" alt="chessboard_arealight_dof_glass_queen" src="https://github.com/user-attachments/assets/d7931b4c-80b7-4985-92e4-fed668e8b5ef" /> | `chessboard_arealight_dof_glass_queen` | 87 ms | 7.74 s |
| <img width="400" height="400" alt="dragon_dynamic" src="https://github.com/user-attachments/assets/32fb532f-224c-489d-a473-ddfa47207d18" /> | `dragon_dynamic` | 3.9 s | 60.8 s |
| <img width="400" height="400" alt="tap_0121" src="https://github.com/user-attachments/assets/0c5c977a-0664-4eae-90c2-4ea88492ad3c" /> | `tap_0121` | 57 ms | 3.1 s |
| <img width="400" height="400" alt="wine_glass" src="https://github.com/user-attachments/assets/a2739043-8c2f-46f2-9bb1-cba6fb70e671" /> | `wine_glass` | 7 ms | 97.9 s |
