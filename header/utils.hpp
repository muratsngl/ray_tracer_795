//contains intersection logic & mathematical operations
#include "types.hpp"
#include "cmath"



//VEC3F OPERATIONS
inline Vec3f operator+(const Vec3f& v0, const Vec3f& v1){
    return Vec3f{v0.x+v1.x,v0.y+v1.y,v0.z+v1.z};
}
inline Vec3f operator-(const Vec3f& v0, const Vec3f& v1){
    return Vec3f{v0.x-v1.x,v0.y-v1.y,v0.z-v1.z};
}
inline Vec3f operator*(float alpha, const Vec3f& v0){
    return {v0.x*alpha,v0.y*alpha,v0.z*alpha};
}
inline Vec3f operator/( const Vec3f& v0,float alpha){
     return {v0.x/alpha,v0.y/alpha,v0.z/alpha};
}

inline void operator+=(Vec3f& v0, const Vec3f& v1){
    v0.x += v1.x;v0.y += v1.y;v0.z += v1.z;
};
inline void operator/=(Vec3f& v0,float alpha){
      v0.x /= alpha; v0.y /= alpha; v0.z /= alpha;
}
inline fl dot(const Vec3f& v0, const Vec3f& v1){
    return v0.x*v1.x+v0.y*v1.y+v0.z*v1.z;
}
inline Vec3f cross(const Vec3f& v0, const Vec3f& v1){
    return Vec3f{v0.y*v1.z-v0.z*v1.y, v0.z*v1.x-v0.x*v1.z, v0.x*v1.y-v0.y*v1.x};
}
inline fl length(const Vec3f& v0){
   return sqrt(v0.x*v0.x+v0.y*v0.y+v0.z*v0.z);
}
//adding an if check for robustness
inline Vec3f normalize(const Vec3f& v0){
    fl len = length(v0);
    return len>1e-6?v0/len:Vec3f{0.f,0.f,0.f};
}
inline void normalize_on_place(Vec3f& v0){
     fl len = length(v0);
     if(len>1e-6f)
     v0/=len;
     else 
     v0 = {0.f,0.f,0.f};
};


//MATRIX OPERATIONS
inline float determinant_3x3(const Mat3f& mat){
     return mat.m11 * (mat.m22 * mat.m33 - mat.m23 * mat.m32)
         - mat.m12 * (mat.m21 * mat.m33 - mat.m23 * mat.m31)
         + mat.m13 * (mat.m21 * mat.m32 - mat.m22 * mat.m31);
}


//GEOMETRY INTERSECTION
inline HitRecord ray_triangle_intersection(const Ray& r, const Triangle& t);
