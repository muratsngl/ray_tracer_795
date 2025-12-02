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


inline void operator+=(ColorBlockFl& c1, ColorBlock& c2){
    for(int i = 0;i<24;i++)c1.rgb[i] += c2.rgb[i]; 
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
inline Mat4f operator*(Mat4f mat1, Mat4f mat2){
    Mat4f res;
    res.m11 = mat1.m11*mat2.m11 + mat1.m12*mat2.m21 + mat1.m13*mat2.m31 + mat1.m14*mat2.m41;
    res.m12 = mat1.m11*mat2.m12 + mat1.m12*mat2.m22 + mat1.m13*mat2.m32 + mat1.m14*mat2.m42;
    res.m13 = mat1.m11*mat2.m13 + mat1.m12*mat2.m23 + mat1.m13*mat2.m33 + mat1.m14*mat2.m43;
    res.m14 = mat1.m11*mat2.m14 + mat1.m12*mat2.m24 + mat1.m13*mat2.m34 + mat1.m14*mat2.m44;
    res.m21 = mat1.m21*mat2.m11 + mat1.m22*mat2.m21 + mat1.m23*mat2.m31 + mat1.m24*mat2.m41;
    res.m22 = mat1.m21*mat2.m12 + mat1.m22*mat2.m22 + mat1.m23*mat2.m32 + mat1.m24*mat2.m42;
    res.m23 = mat1.m21*mat2.m13 + mat1.m22*mat2.m23 + mat1.m23*mat2.m33 + mat1.m24*mat2.m43;
    res.m24 = mat1.m21*mat2.m14 + mat1.m22*mat2.m24 + mat1.m23*mat2.m34 + mat1.m24*mat2.m44;
    res.m31 = mat1.m31*mat2.m11 + mat1.m32*mat2.m21 + mat1.m33*mat2.m31 + mat1.m34*mat2.m41;
    res.m32 = mat1.m31*mat2.m12 + mat1.m32*mat2.m22 + mat1.m33*mat2.m32 + mat1.m34*mat2.m42;
    res.m33 = mat1.m31*mat2.m13 + mat1.m32*mat2.m23 + mat1.m33*mat2.m33 + mat1.m34*mat2.m43;
    res.m34 = mat1.m31*mat2.m14 + mat1.m32*mat2.m24 + mat1.m33*mat2.m34 + mat1.m34*mat2.m44;
    res.m41 = mat1.m41*mat2.m11 + mat1.m42*mat2.m21 + mat1.m43*mat2.m31 + mat1.m44*mat2.m41;
    res.m42 = mat1.m41*mat2.m12 + mat1.m42*mat2.m22 + mat1.m43*mat2.m32 + mat1.m44*mat2.m42;
    res.m43 = mat1.m41*mat2.m13 + mat1.m42*mat2.m23 + mat1.m43*mat2.m33 + mat1.m44*mat2.m43;
    res.m44 = mat1.m41*mat2.m14 + mat1.m42*mat2.m24 + mat1.m43*mat2.m34 + mat1.m44*mat2.m44;
    return res;
}
inline Mat4f mat_inv(Mat4f mat1){
    fl det = mat1.m11 * (mat1.m22 * (mat1.m33 * mat1.m44 - mat1.m34 * mat1.m43) - mat1.m23 * (mat1.m32 * mat1.m44 - mat1.m34 * mat1.m42) + mat1.m24 * (mat1.m32 * mat1.m43 - mat1.m33 * mat1.m42))
          - mat1.m12 * (mat1.m21 * (mat1.m33 * mat1.m44 - mat1.m34 * mat1.m43) - mat1.m23 * (mat1.m31 * mat1.m44 - mat1.m34 * mat1.m41) + mat1.m24 * (mat1.m31 * mat1.m43 - mat1.m33 * mat1.m41))
          + mat1.m13 * (mat1.m21 * (mat1.m32 * mat1.m44 - mat1.m34 * mat1.m42) - mat1.m22 * (mat1.m31 * mat1.m44 - mat1.m34 * mat1.m41) + mat1.m24 * (mat1.m31 * mat1.m42 - mat1.m32 * mat1.m41))
          - mat1.m14 * (mat1.m21 * (mat1.m32 * mat1.m43 - mat1.m33 * mat1.m42) - mat1.m22 * (mat1.m31 * mat1.m43 - mat1.m33 * mat1.m41) + mat1.m23 * (mat1.m31 * mat1.m42 - mat1.m32 * mat1.m41));
    if (std::abs(det) < 1e-6) {
        return Mat4f{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    }
    fl inv_det = 1.0f / det;
    Mat4f inv;
    
    // M_ij = C_ji
    inv.m11 = inv_det * (mat1.m22 * (mat1.m33 * mat1.m44 - mat1.m34 * mat1.m43) - mat1.m23 * (mat1.m32 * mat1.m44 - mat1.m34 * mat1.m42) + mat1.m24 * (mat1.m32 * mat1.m43 - mat1.m33 * mat1.m42));
    inv.m21 = inv_det * - (mat1.m21 * (mat1.m33 * mat1.m44 - mat1.m34 * mat1.m43) - mat1.m23 * (mat1.m31 * mat1.m44 - mat1.m34 * mat1.m41) + mat1.m24 * (mat1.m31 * mat1.m43 - mat1.m33 * mat1.m41));
    inv.m31 = inv_det * (mat1.m21 * (mat1.m32 * mat1.m44 - mat1.m34 * mat1.m42) - mat1.m22 * (mat1.m31 * mat1.m44 - mat1.m34 * mat1.m41) + mat1.m24 * (mat1.m31 * mat1.m42 - mat1.m32 * mat1.m41));
    inv.m41 = inv_det * - (mat1.m21 * (mat1.m32 * mat1.m43 - mat1.m33 * mat1.m42) - mat1.m22 * (mat1.m31 * mat1.m43 - mat1.m33 * mat1.m41) + mat1.m23 * (mat1.m31 * mat1.m42 - mat1.m32 * mat1.m41));
    
    inv.m12 = inv_det * - (mat1.m12 * (mat1.m33 * mat1.m44 - mat1.m34 * mat1.m43) - mat1.m13 * (mat1.m32 * mat1.m44 - mat1.m34 * mat1.m42) + mat1.m14 * (mat1.m32 * mat1.m43 - mat1.m33 * mat1.m42));
    inv.m22 = inv_det * (mat1.m11 * (mat1.m33 * mat1.m44 - mat1.m34 * mat1.m43) - mat1.m13 * (mat1.m31 * mat1.m44 - mat1.m34 * mat1.m41) + mat1.m14 * (mat1.m31 * mat1.m43 - mat1.m33 * mat1.m41));
    inv.m32 = inv_det * - (mat1.m11 * (mat1.m32 * mat1.m44 - mat1.m34 * mat1.m42) - mat1.m12 * (mat1.m31 * mat1.m44 - mat1.m34 * mat1.m41) + mat1.m14 * (mat1.m31 * mat1.m42 - mat1.m32 * mat1.m41));
    inv.m42 = inv_det * (mat1.m11 * (mat1.m32 * mat1.m43 - mat1.m33 * mat1.m42) - mat1.m12 * (mat1.m31 * mat1.m43 - mat1.m33 * mat1.m41) + mat1.m13 * (mat1.m31 * mat1.m42 - mat1.m32 * mat1.m41));
    
    inv.m13 = inv_det * (mat1.m12 * (mat1.m23 * mat1.m44 - mat1.m24 * mat1.m43) - mat1.m13 * (mat1.m22 * mat1.m44 - mat1.m24 * mat1.m42) + mat1.m14 * (mat1.m22 * mat1.m43 - mat1.m23 * mat1.m42));
    inv.m23 = inv_det * - (mat1.m11 * (mat1.m23 * mat1.m44 - mat1.m24 * mat1.m43) - mat1.m13 * (mat1.m21 * mat1.m44 - mat1.m24 * mat1.m41) + mat1.m14 * (mat1.m21 * mat1.m43 - mat1.m23 * mat1.m41));
    inv.m33 = inv_det * (mat1.m11 * (mat1.m22 * mat1.m44 - mat1.m24 * mat1.m42) - mat1.m12 * (mat1.m21 * mat1.m44 - mat1.m24 * mat1.m41) + mat1.m14 * (mat1.m21 * mat1.m42 - mat1.m22 * mat1.m41));
    inv.m43 = inv_det * - (mat1.m11 * (mat1.m22 * mat1.m43 - mat1.m23 * mat1.m42) - mat1.m12 * (mat1.m21 * mat1.m43 - mat1.m23 * mat1.m41) + mat1.m13 * (mat1.m21 * mat1.m42 - mat1.m22 * mat1.m41));
    
    inv.m14 = inv_det * - (mat1.m12 * (mat1.m23 * mat1.m34 - mat1.m24 * mat1.m33) - mat1.m13 * (mat1.m22 * mat1.m34 - mat1.m24 * mat1.m32) + mat1.m14 * (mat1.m22 * mat1.m33 - mat1.m23 * mat1.m32));
    inv.m24 = inv_det * (mat1.m11 * (mat1.m23 * mat1.m34 - mat1.m24 * mat1.m33) - mat1.m13 * (mat1.m21 * mat1.m34 - mat1.m24 * mat1.m31) + mat1.m14 * (mat1.m21 * mat1.m33 - mat1.m23 * mat1.m31));
    inv.m34 = inv_det * - (mat1.m11 * (mat1.m22 * mat1.m34 - mat1.m24 * mat1.m32) - mat1.m12 * (mat1.m21 * mat1.m34 - mat1.m24 * mat1.m31) + mat1.m14 * (mat1.m21 * mat1.m32 - mat1.m22 * mat1.m31));
    inv.m44 = inv_det * (mat1.m11 * (mat1.m22 * mat1.m33 - mat1.m23 * mat1.m32) - mat1.m12 * (mat1.m21 * mat1.m33 - mat1.m23 * mat1.m31) + mat1.m13 * (mat1.m21 * mat1.m32 - mat1.m22 * mat1.m31));
    
    return inv;
}
// GEMINI created code parts
inline Mat4f mat_transpose(Mat4f mat) {
    std::swap(mat.m12, mat.m21);
    std::swap(mat.m13, mat.m31);
    std::swap(mat.m14, mat.m41);
    std::swap(mat.m23, mat.m32);
    std::swap(mat.m24, mat.m42);
    std::swap(mat.m34, mat.m43);
    return mat;
}
// GEMINI created code parts
inline Vec4f operator*(Mat4f mat1,Vec4f vec1){
    Vec4f res;
    res.x = mat1.m11 * vec1.x + mat1.m12 * vec1.y + mat1.m13 * vec1.z + mat1.m14 * vec1.w;
    res.y = mat1.m21 * vec1.x + mat1.m22 * vec1.y + mat1.m23 * vec1.z + mat1.m24 * vec1.w;
    res.z = mat1.m31 * vec1.x + mat1.m32 * vec1.y + mat1.m33 * vec1.z + mat1.m34 * vec1.w;
    res.w = mat1.m41 * vec1.x + mat1.m42 * vec1.y + mat1.m43 * vec1.z + mat1.m44 * vec1.w;
    return res;
}
// GEMINI created code parts
inline Mat4f create_identity_matrix() {
    Mat4f mat;
    mat.m11 = 1.0f; mat.m22 = 1.0f; mat.m33 = 1.0f; mat.m44 = 1.0f;
    return mat;
}

inline Mat4f create_translation_matrix(Vec3f t) {
    Mat4f mat = create_identity_matrix();
    mat.m14 = t.x;
    mat.m24 = t.y;
    mat.m34 = t.z;
    return mat;
}

inline Mat4f create_scaling_matrix(Vec3f s) {
    Mat4f mat = create_identity_matrix();
    mat.m11 = s.x;
    mat.m22 = s.y;
    mat.m33 = s.z;
    return mat;
}

inline Mat4f create_rotation_matrix(fl angle_degrees, Vec3f axis) {
    Mat4f mat = create_identity_matrix();
    fl angle_rad = angle_degrees * (M_PI / 180.0f);
    normalize_on_place(axis);
    fl c = cos(angle_rad);
    fl s = sin(angle_rad);
    fl t = 1.0f - c;
    fl x = axis.x, y = axis.y, z = axis.z;

    mat.m11 = t*x*x + c;   mat.m12 = t*x*y - s*z; mat.m13 = t*x*z + s*y;
    mat.m21 = t*x*y + s*z; mat.m22 = t*y*y + c;   mat.m23 = t*y*z - s*x;
    mat.m31 = t*x*z - s*y; mat.m32 = t*y*z + s*x; mat.m33 = t*z*z + c;

    return mat;
}
// GEMINI created code parts


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