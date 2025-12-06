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

