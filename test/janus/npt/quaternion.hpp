#ifndef HPP_QUARTERNION
#define HPP_QUARTERNION

#include<iostream>
#include <cmath>

#include "matrix3.hpp"
#include "matrix4.hpp"

template <class T>
Matrix4<double> S(const T a){
  return Matrix4<double>(a[0],-a[1],-a[2],-a[3],
			 a[1], a[0], a[3],-a[2],
			 a[2],-a[3], a[0], a[1],
			 a[3], a[2],-a[1], a[0]);
}


class Quaternion{
public:
  double data[4];

  Quaternion(){data[0]=1.0; data[1] = data[2] = data[3] = 0.0;}
  Quaternion(const double& _s,const dvec3& _v){
    data[0]=_s; data[1] = _v[0]; data[2] = _v[1]; data[3] = _v[2];
  }
  Quaternion(const double& _x,const double& _y,const double& _z, const double& _w){
    data[0] = _x; data[1] = _y; data[2] = _z; data[3] = _w;
  }
  Quaternion(const Quaternion& src){
    data[0] = src[0]; data[1] = src[1]; data[2] = src[2]; data[3] = src[3];}
  Quaternion(const double& src){
    data[0] = data[1] = data[2] = data[3] = src;
  }

  Quaternion operator+(const Quaternion& rhs) const {
    return Quaternion(data[0]+rhs[0],
		      data[1]+rhs[1],
		      data[2]+rhs[2],
		      data[3]+rhs[3]);
  }
  const Quaternion& operator+=(const Quaternion& rhs){
    (*this) = (*this) + rhs;
    return (*this);
  }
  Quaternion operator-(const Quaternion& rhs) const {
    return Quaternion(data[0]-rhs[0],
		      data[1]-rhs[1],
		      data[2]-rhs[2],
		      data[3]-rhs[3]);
  }
  const Quaternion& operator-=(const Quaternion& rhs){
    (*this) = (*this) - rhs;
    return (*this);
  }
  Quaternion operator*(const Quaternion& rhs) const {
    const double s0 = data[0];
    const double s1 = rhs[0];
    const dvec3 v0  = dvec3(data[1],data[2],data[3]);
    const dvec3 v1  = dvec3(rhs[1],rhs[2],rhs[3]);
    return Quaternion(s0 * s1 - v0 % v1, v1*s0 + v0*s1 + (v0^v1));
  }
  Quaternion operator*(const dvec3& _rhs) const { 
    Quaternion rhs(0.0,_rhs);
    return (*this)*rhs;
  }
  friend Quaternion operator*(const dvec3& _lhs,const Quaternion& rhs){
    Quaternion lhs(0.0,_lhs);
    return lhs*rhs;
  }
  Quaternion operator*(const double& rhs) const { 
    return Quaternion(data[0]*rhs,
		      data[1]*rhs,
		      data[2]*rhs,
		      data[3]*rhs);
  }
  const Quaternion& operator*=(const double& rhs){
    (*this) = (*this) * rhs;
    return (*this);
  }
  friend Quaternion operator*(const Matrix4<double>& lhs,const Quaternion& rhs){ 
    return Quaternion(Quaternion(lhs.xx,lhs.xy,lhs.xz,lhs.xw) % rhs,
		      Quaternion(lhs.yx,lhs.yy,lhs.yz,lhs.yw) % rhs,
		      Quaternion(lhs.zx,lhs.zy,lhs.zz,lhs.zw) % rhs,
		      Quaternion(lhs.wx,lhs.wy,lhs.wz,lhs.ww) % rhs);
  }
  double operator%(const Quaternion& rhs) const {
    return data[0]*rhs[0] + data[1]*rhs[1] + data[2]*rhs[2] + data[3]*rhs[3];
  }

  const Quaternion& operator=(const Quaternion& rhs){ 
    data[0] = rhs[0];
    data[1] = rhs[1];
    data[2] = rhs[2];
    data[3] = rhs[3];

    return (*this);
  }
  const Quaternion& operator=(const double& rhs){
    data[0] = rhs;
    data[1] = rhs;
    data[2] = rhs;
    data[3] = rhs;

    return (*this);
  }

  const double& operator[](const int i) const {return data[i];}
        double& operator[](const int i)       {return data[i];}

  friend Quaternion normalize(const Quaternion& q){
    double ni = 1.0/norm(q);
    return Quaternion(q[0]*ni,
		      q[1]*ni,
		      q[2]*ni,
		      q[3]*ni);
  }

  friend Quaternion conj(const Quaternion& q){
    return Quaternion( q[0],
		      -q[1],
		      -q[2],
		      -q[3]);
  }
  friend double norm(const Quaternion& q){return sqrt(q%q);}

  friend Matrix3<double> Rotate(const Quaternion& q){
    return Matrix3<double>(q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
			   2.0*(q[1]*q[2] - q[0]*q[3]),
			   2.0*(q[1]*q[3] + q[0]*q[2]),

			   2.0*(q[1]*q[2] + q[0]*q[3]),
			   q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3],
			   2.0*(q[2]*q[3] - q[0]*q[1]),

			   2.0*(q[1]*q[3] - q[0]*q[2]),
			   2.0*(q[2]*q[3] + q[0]*q[1]),
			   q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3]);
  }
};

#endif
