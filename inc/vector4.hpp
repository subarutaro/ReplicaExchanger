#ifndef H_VECTOR4
#define H_VECTOR4

#ifdef __CUDACC__
#define __ATTRIBUTES__ __host__ __device__ __forceinline__
#else
#define __ATTRIBUTES__ inline
#endif

template <class T>
class Vector4{
public:
  T x,y,z,w;

  __ATTRIBUTES__
  Vector4():x((T)0),y((T)0),z((T)0),w((T)0){}
  __ATTRIBUTES__
  Vector4(const T s):x(s),y(s),z(s),w(s){}
  __ATTRIBUTES__
  Vector4(const T _x,const T _y,const T _z,const T _w):x(_x),y(_y),z(_z),w(_w){}
  __ATTRIBUTES__
  Vector4(const T xyzw[4]):x(xyzw[0]),y(xyzw[1]),z(xyzw[2]),w(xyzw[3]){}
  __ATTRIBUTES__
  Vector4(const Vector4 &r):x(r.x),y(r.y),z(r.z),w(r.w){}

  __ATTRIBUTES__
  const Vector4& operator = (const Vector4& rhs){
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    w = rhs.w;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector4& operator = (const T& s){
    x = y = z = w = s;
    return (*this);
  }

  // vector x vextor
  __ATTRIBUTES__
  Vector4 operator + (const Vector4& rhs) const {
    return Vector4(x + rhs.x,
		   y + rhs.y,
		   z + rhs.z,
		   w + rhs.w);
  }
  __ATTRIBUTES__
  Vector4 operator - (const Vector4& rhs) const {
    return Vector4(x - rhs.x,
		   y - rhs.y,
		   z - rhs.z,
		   w - rhs.w);
  }
  __ATTRIBUTES__
  Vector4 operator * (const Vector4& rhs) const {
    return Vector4(x * rhs.x,
		   y * rhs.y,
		   z * rhs.z,
		   w * rhs.w);
  }
  __ATTRIBUTES__
  Vector4 operator / (const Vector4& rhs) const {
    return Vector4(x / rhs.x,
		   y / rhs.y,
		   z / rhs.z,
		   w / rhs.w);
  }
  __ATTRIBUTES__
  const Vector4& operator += (const Vector4& rhs){
    (*this) = (*this) + rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector4& operator -= (const Vector4& rhs){
    (*this) = (*this) - rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector4& operator *= (const Vector4& rhs){
    (*this) = (*this) * rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector4& operator /= (const Vector4& rhs){
    (*this) = (*this) / rhs;
    return (*this);
  }
  // vector x scholar
  __ATTRIBUTES__
  Vector4 operator + (const T& rhs) const {
    return Vector4(x + rhs,
		   y + rhs,
		   z + rhs,
		   w + rhs);
  }
  __ATTRIBUTES__
  Vector4 operator - (const T& rhs) const {
    return Vector4(x - rhs,
		   y - rhs,
		   z - rhs,
		   w - rhs);
  }
  __ATTRIBUTES__
  Vector4 operator * (const T& rhs) const {
    return Vector4(x * rhs,
		   y * rhs,
		   z * rhs,
		   w * rhs);
  }
  __ATTRIBUTES__
  friend Vector4 operator * (T& lhs,const Vector4& rhs){
    return rhs*lhs;
  }
  __ATTRIBUTES__
  Vector4 operator / (const T& rhs) const {
    return Vector4(x / rhs,
		   y / rhs,
		   z / rhs,
		   w / rhs);
  }
  __ATTRIBUTES__
  const Vector4& operator += (const T& rhs){
    (*this) = (*this) + rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector4& operator -= (const T& rhs){
    (*this) = (*this) - rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector4& operator *= (const T& rhs){
    (*this) = (*this) * rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector4& operator /= (const T& rhs){
    (*this) = (*this) / rhs;
    return (*this);
  }


  // inner product
  __ATTRIBUTES__
  T operator % (const Vector4& rhs) const {
    return (x*rhs.x) + (y*rhs.y) + (z*rhs.z) + (w*rhs.w);
  }

  //cast
  template <typename U>
  __ATTRIBUTES__
  operator Vector4<U> () const {
    return Vector4<U> (static_cast<U>(x),
		       static_cast<U>(y),
		       static_cast<U>(z),
		       static_cast<U>(w));
  }
  __ATTRIBUTES__
  T max() const {
    const T max1 = x>y ? x : y;
    const T max2 = z>w ? z : w;
    return max1>max2 ? max1 : max2;
  }
  __ATTRIBUTES__
  T min() const {
    const T min1 = x<y ? x : y;
    const T min2 = z<w ? z : w;
    return min1<min2 ? min1 : min2;
  }
  __ATTRIBUTES__
  const T& operator[](const int i) const {
    if(i>=4 || i<0){
      std::cerr << "error: invalid index(= " << i <<  ") for Vector4." << std::endl;
      exit(-1);
    }
    if(i==0) return x;
    if(i==1) return y;
    if(i==2) return z;
    if(i==3) return w;
  }
  __ATTRIBUTES__
  T& operator[](const int i){
    if(i>=4 || i<0){
      std::cerr << "error: invalid index(= " << i <<  ") for Vector4." << std::endl;
      exit(-1);
    }
    if(i==0) return x;
    if(i==1) return y;
    if(i==2) return z;
    if(i==3) return w;
  }
  __ATTRIBUTES__
  friend T norm(const Vector4& v){
    return sqrt(v%v);
  }
  __ATTRIBUTES__
  bool operator == (const Vector4& rhs) const {
    return ((x==rhs.x) || (y==rhs.y) || (z==rhs.z) || (w==rhs.w));
  }
  __ATTRIBUTES__
  bool operator != (const Vector4& rhs) const {
    return ((x!=rhs.x) || (y!=rhs.y) || (z!=rhs.z) || (w!=rhs.w));
  }
};

template <>
__ATTRIBUTES__
Vector4<float> Vector4<float>::operator / (const float& s) const {
  const float si = 1.f / s;
  return Vector4(x*si, y*si, z*si, w*si);
}

template <>
__ATTRIBUTES__
Vector4<double> Vector4<double>::operator / (const double& s) const {
  const double si = 1.f / s;
  return Vector4(x*si, y*si, z*si, w*si);
}

#undef __ATTRIBUTES__

#endif
