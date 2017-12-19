#ifndef H_VECTOR2
#define H_VECTOR2

#ifdef __CUDACC__
#define __ATTRIBUTES__ __host__ __device__ __forceinline__
#else
#define __ATTRIBUTES__ inline
#endif

#include <cstdlib>

template <class T>
class Vector2{
public:
  T x,y;

  __ATTRIBUTES__
  Vector2():x((T)0),y((T)0){}
  __ATTRIBUTES__
  Vector2(const T s):x(s),y(s){}
  __ATTRIBUTES__
  Vector2(const T _x,const T _y):x(_x),y(_y){}
  __ATTRIBUTES__
  Vector2(const T xy[4]):x(xy[0]),y(xy[1]){}
  __ATTRIBUTES__
  Vector2(const Vector2 &r):x(r.x),y(r.y){}

  __ATTRIBUTES__
  const Vector2& operator = (const Vector2& rhs){
    x = rhs.x;
    y = rhs.y;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector2& operator = (const T& s){
    x = y = s;
    return (*this);
  }

  // vector x vextor
  __ATTRIBUTES__
  Vector2 operator + (const Vector2& rhs) const {
    return Vector2(x + rhs.x,
		   y + rhs.y);
  }
  __ATTRIBUTES__
  Vector2 operator - (const Vector2& rhs) const {
    return Vector2(x - rhs.x,
		   y - rhs.y);
  }
  __ATTRIBUTES__
  Vector2 operator * (const Vector2& rhs) const {
    return Vector2(x * rhs.x,
		   y * rhs.y);
  }
  __ATTRIBUTES__
  Vector2 operator / (const Vector2& rhs) const {
    return Vector2(x / rhs.x,
		   y / rhs.y);
  }
  __ATTRIBUTES__
  const Vector2& operator += (const Vector2& rhs){
    (*this) = (*this) + rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector2& operator -= (const Vector2& rhs){
    (*this) = (*this) - rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector2& operator *= (const Vector2& rhs){
    (*this) = (*this) * rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector2& operator /= (const Vector2& rhs){
    (*this) = (*this) / rhs;
    return (*this);
  }
  // vector x scholar
  __ATTRIBUTES__
  Vector2 operator + (const T& rhs) const {
    return Vector2(x + rhs,
		   y + rhs);
  }
  __ATTRIBUTES__
  Vector2 operator - (const T& rhs) const {
    return Vector2(x - rhs,
		   y - rhs);
  }
  __ATTRIBUTES__
  Vector2 operator * (const T& rhs) const {
    return Vector2(x * rhs,
		   y * rhs);
  }
  __ATTRIBUTES__
  friend Vector2 operator * (T& lhs,const Vector2& rhs){
    return rhs*lhs;
  }
  __ATTRIBUTES__
  Vector2 operator / (const T& rhs) const {
    return Vector2(x / rhs,
		   y / rhs);
  }
  __ATTRIBUTES__
  const Vector2& operator += (const T& rhs){
    (*this) = (*this) + rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector2& operator -= (const T& rhs){
    (*this) = (*this) - rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector2& operator *= (const T& rhs){
    (*this) = (*this) * rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector2& operator /= (const T& rhs){
    (*this) = (*this) / rhs;
    return (*this);
  }


  // inner product
  __ATTRIBUTES__
  T operator % (const Vector2& rhs) const {
    return (x*rhs.x) + (y*rhs.y);
  }

  //cast
  template <typename U>
  __ATTRIBUTES__
  operator Vector2<U> () const {
    return Vector2<U> (static_cast<U>(x),
		       static_cast<U>(y));
  }
  __ATTRIBUTES__
  T max() const {
    return x>y ? x : y;
  }
  __ATTRIBUTES__
  T min() const {
    return x<y ? x : y;
  }
  __ATTRIBUTES__
  const T& operator[](const int i) const {
    if(i>=2 || i<0){
      std::cerr << "error: invalid index(= " << i <<  ") for Vector2." << std::endl;
      exit(-1);
    }
    if(i==0) return x;
    if(i==1) return y;
  }
  __ATTRIBUTES__
  T& operator[](const int i){
    if(i>=2 || i<0){
      std::cerr << "error: invalid index(= " << i <<  ") for Vector2." << std::endl;
      exit(-1);
    }
    if(i==0) return x;
    if(i==1) return y;
  }
  __ATTRIBUTES__
  friend T norm(const Vector2& v){
    return sqrt(v%v);
  }
  __ATTRIBUTES__
  bool operator == (const Vector2& rhs) const {
    return ((x==rhs.x) || (y==rhs.y));
  }
  __ATTRIBUTES__
  bool operator != (const Vector2& rhs) const {
    return ((x!=rhs.x) || (y!=rhs.y));
  }
};

template <>
__ATTRIBUTES__
Vector2<float> Vector2<float>::operator / (const float& s) const {
  const float si = 1.f / s;
  return Vector2(x*si, y*si);
}

template <>
__ATTRIBUTES__
Vector2<double> Vector2<double>::operator / (const double& s) const {
  const double si = 1.f / s;
  return Vector2(x*si, y*si);
}

#undef __ATTRIBUTES__

#endif
