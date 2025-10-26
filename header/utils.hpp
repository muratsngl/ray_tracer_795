#ifndef UTILS_HPP
#define UTILS_HPP
//contains intersection logic & mathematical operations
#include "types.hpp"
#include "cmath"
#include <xsimd/xsimd.hpp>


namespace xs = xsimd;

// Acronym type aliases for xsimd batches
using f_batch = xs::batch<fl>;
using i_batch = xs::batch<int>;
using b_batch = xs::batch_bool<fl>;


//VEC3F OPERATIONS (kept for backward compatibility)
inline Vec3f operator+(const Vec3f& v0, const Vec3f& v1){
    return Vec3f{v0.x+v1.x,v0.y+v1.y,v0.z+v1.z};
}
//binary
inline Vec3f operator-(const Vec3f& v0, const Vec3f& v1){
    return Vec3f{v0.x-v1.x,v0.y-v1.y,v0.z-v1.z};
}
//unary
inline Vec3f operator-(const Vec3f& v0){
return  Vec3f{-v0.x,-v0.y,-v0.z};
}
inline Vec3f operator*(fl alpha, const Vec3f& v0){
    return {v0.x*alpha,v0.y*alpha,v0.z*alpha};
}
inline Vec3f operator/( const Vec3f& v0,fl alpha){
     return {v0.x/alpha,v0.y/alpha,v0.z/alpha};
}

inline void operator+=(Vec3f& v0, const Vec3f& v1){
    v0.x += v1.x;v0.y += v1.y;v0.z += v1.z;
};
inline void operator/=(Vec3f& v0,fl alpha){
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


//SCALAR-BASED VECTOR OPERATIONS (for SoA data structures)
// Dot product: returns x0*x1 + y0*y1 + z0*z1
inline fl dot_scalar(fl x0, fl y0, fl z0, fl x1, fl y1, fl z1){
    return x0*x1 + y0*y1 + z0*z1;
}

// Cross product: outputs the cross product components
inline void cross_scalar(fl x0, fl y0, fl z0, fl x1, fl y1, fl z1, 
                         fl& out_x, fl& out_y, fl& out_z){
    out_x = y0*z1 - z0*y1;
    out_y = z0*x1 - x0*z1;
    out_z = x0*y1 - y0*x1;
}

// Vector subtraction: v0 - v1
inline void subtract_scalar(fl x0, fl y0, fl z0, fl x1, fl y1, fl z1,
                           fl& out_x, fl& out_y, fl& out_z){
    out_x = x0 - x1;
    out_y = y0 - y1;
    out_z = z0 - z1;
}

// Vector addition: v0 + v1
inline void add_scalar(fl x0, fl y0, fl z0, fl x1, fl y1, fl z1,
                      fl& out_x, fl& out_y, fl& out_z){
    out_x = x0 + x1;
    out_y = y0 + y1;
    out_z = z0 + z1;
}

// Scalar multiplication: alpha * v
inline void multiply_scalar(fl alpha, fl x, fl y, fl z,
                           fl& out_x, fl& out_y, fl& out_z){
    out_x = alpha * x;
    out_y = alpha * y;
    out_z = alpha * z;
}

// Length (magnitude) of vector
inline fl length_scalar(fl x, fl y, fl z){
    return sqrt(x*x + y*y + z*z);
}

// Normalize vector (in-place)
inline void normalize_scalar(fl& x, fl& y, fl& z){
    fl len = length_scalar(x, y, z);
    if(len > 1e-6f){
        x /= len;
        y /= len;
        z /= len;
    } else {
        x = y = z = 0.0f;
    }
}

// Normalize vector (returns new values)
inline void normalize_scalar_copy(fl x, fl y, fl z,
                                  fl& out_x, fl& out_y, fl& out_z){
    fl len = length_scalar(x, y, z);
    if(len > 1e-6f){
        out_x = x / len;
        out_y = y / len;
        out_z = z / len;
    } else {
        out_x = out_y = out_z = 0.0f;
    }
}


//MATRIX OPERATIONS
inline fl determinant_3x3(const Mat3f& mat){
     return mat.m11 * (mat.m22 * mat.m33 - mat.m23 * mat.m32)
         - mat.m12 * (mat.m21 * mat.m33 - mat.m23 * mat.m31)
         + mat.m13 * (mat.m21 * mat.m32 - mat.m22 * mat.m31);
}

inline f_batch length_simd(f_batch& vec0_x,f_batch& vec0_y,f_batch& vec0_z){
    return xs::sqrt(vec0_x * vec0_x + vec0_y*vec0_y + vec0_z*vec0_z);
}
inline void normalize_simd_overwrite(f_batch& vec0_x,f_batch& vec0_y,f_batch& vec0_z){
    f_batch accumulator = length_simd(vec0_x,vec0_y,vec0_z);
    vec0_x/=accumulator;
    vec0_y/=accumulator;
    vec0_z/=accumulator;
}
inline void cross_simd_copy(const f_batch& vec0_x,const f_batch& vec0_y,const f_batch& vec0_z,const f_batch& vec1_x,const f_batch& vec1_y,const f_batch& vec1_z,f_batch& rvec_x,f_batch& rvec_y,f_batch& rvec_z){
    rvec_x = vec0_y*vec1_z-vec1_y*vec0_z;
    rvec_y = vec0_x*vec1_z-vec1_x*vec0_z;
    rvec_z = vec0_x*vec1_y-vec1_x*vec0_y;

}

inline void normmalize_simd_copy();


inline f_batch dot_simd(const f_batch& vec0_x,const f_batch& vec0_y,const f_batch& vec0_z,const f_batch& vec1_x,const f_batch& vec1_y,const f_batch& vec1_z){
        return vec0_x*vec1_x+vec0_y*vec1_y+vec0_z*vec1_z;
}
#endif