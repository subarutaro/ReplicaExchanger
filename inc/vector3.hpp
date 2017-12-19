#ifndef H_VECTOR3
#define H_VECTOR3

#ifdef __CUDACC__
#define __ATTRIBUTES__ __host__ __device__ __forceinline__
#else
#define __ATTRIBUTES__ inline
#endif

template <class T>
class Vector3{
public:
  T x,y,z,w;

  __ATTRIBUTES__
  Vector3():x((T)0),y((T)0),z((T)0){}
  __ATTRIBUTES__
  Vector3(const T s):x(s),y(s),z(s){}
  __ATTRIBUTES__
  Vector3(const T _x,const T _y,const T _z):x(_x),y(_y),z(_z){}
  __ATTRIBUTES__
  Vector3(const T xyz[4]):x(xyz[0]),y(xyz[1]),z(xyz[2]){}
  __ATTRIBUTES__
  Vector3(const Vector3 &r):x(r.x),y(r.y),z(r.z){}

  __ATTRIBUTES__
  const Vector3& operator = (const Vector3& rhs){
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector3& operator = (const T& s){
    x = y = z = s;
    return (*this);
  }

  // vector x vextor
  __ATTRIBUTES__
  Vector3 operator + (const Vector3& rhs) const {
    return Vector3(x + rhs.x,
		   y + rhs.y,
		   z + rhs.z);
  }
  __ATTRIBUTES__
  Vector3 operator - (const Vector3& rhs) const {
    return Vector3(x - rhs.x,
		   y - rhs.y,
		   z - rhs.z);
  }
  __ATTRIBUTES__
  Vector3 operator * (const Vector3& rhs) const {
    return Vector3(x * rhs.x,
		   y * rhs.y,
		   z * rhs.z);
  }
  __ATTRIBUTES__
  Vector3 operator / (const Vector3& rhs) const {
    return Vector3(x / rhs.x,
		   y / rhs.y,
		   z / rhs.z);
  }
  __ATTRIBUTES__
  const Vector3& operator += (const Vector3& rhs){
    (*this) = (*this) + rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector3& operator -= (const Vector3& rhs){
    (*this) = (*this) - rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector3& operator *= (const Vector3& rhs){
    (*this) = (*this) * rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector3& operator /= (const Vector3& rhs){
    (*this) = (*this) / rhs;
    return (*this);
  }
  // vector x scholar
  __ATTRIBUTES__
  Vector3 operator + (const T& rhs) const {
    return Vector3(x + rhs,
		   y + rhs,
		   z + rhs);
  }
  __ATTRIBUTES__
  Vector3 operator - (const T& rhs) const {
    return Vector3(x - rhs,
		   y - rhs,
		   z - rhs);
  }
  __ATTRIBUTES__
  Vector3 operator * (const T& rhs) const {
    return Vector3(x * rhs,
		   y * rhs,
		   z * rhs);
  }
  __ATTRIBUTES__
  friend Vector3 operator * (const T& lhs,const Vector3& rhs){
    return rhs*lhs;
  }
  __ATTRIBUTES__
  Vector3 operator / (const T& rhs) const {
    return Vector3(x / rhs,
		   y / rhs,
		   z / rhs);
  }
  __ATTRIBUTES__
  const Vector3& operator += (const T& rhs){
    (*this) = (*this) + rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector3& operator -= (const T& rhs){
    (*this) = (*this) - rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector3& operator *= (const T& rhs){
    (*this) = (*this) * rhs;
    return (*this);
  }
  __ATTRIBUTES__
  const Vector3& operator /= (const T& rhs){
    (*this) = (*this) / rhs;
    return (*this);
  }


  // inner product
  __ATTRIBUTES__
  T operator % (const Vector3& rhs) const {
    return (x*rhs.x) + (y*rhs.y) + (z*rhs.z);
  }
  // outer product
  __ATTRIBUTES__
  T operator ^ (const Vector3& rhs) const {
    return Vector3(y * rhs.z - z * rhs.y,
		   z * rhs.x - x * rhs.z,
		   x * rhs.y - y * rhs.x);
  }

  //cast
  template <typename U>
  __ATTRIBUTES__
  operator Vector3<U> () const {
    return Vector3<U> (static_cast<U>(x),
		       static_cast<U>(y),
		       static_cast<U>(z));
  }
  __ATTRIBUTES__
  T max() const {
    const T max = x>y ? x : y;
    return max>z ? max : z;
  }
  __ATTRIBUTES__
  T min() const {
    const T min = x<y ? x : y;
    return min<z ? min : z;
  }
  __ATTRIBUTES__
  const T& operator[](const int i) const {
    if(i>=3 || i<0){
      std::cerr << "error: invalid index(= " << i <<  ") for Vector3." << std::endl;
      exit(-1);
    }
    if(i==0) return x;
    if(i==1) return y;
    if(i==2) return z;
  }
  __ATTRIBUTES__
  T& operator[](const int i){
    if(i>=3 || i<0){
      std::cerr << "error: invalid index(= " << i <<  ") for Vector3." << std::endl;
      exit(-1);
    }
    if(i==0) return x;
    if(i==1) return y;
    if(i==2) return z;
  }
  __ATTRIBUTES__
  friend T norm(const Vector3& v){
    return sqrt(v%v);
  }
  __ATTRIBUTES__
  bool operator == (const Vector3& rhs) const {
    return ((x==rhs.x) || (y==rhs.y) || (z==rhs.z));
  }
  __ATTRIBUTES__
  bool operator != (const Vector3& rhs) const {
    return ((x!=rhs.x) || (y!=rhs.y) || (z!=rhs.z));
  }
};

template <>
__ATTRIBUTES__
Vector3<float> Vector3<float>::operator / (const float& s) const {
  const float si = 1.f / s;
  return Vector3(x*si, y*si, z*si);
}

template <>
__ATTRIBUTES__
Vector3<double> Vector3<double>::operator / (const double& s) const {
  const double si = 1.f / s;
  return Vector3(x*si, y*si, z*si);
}

#undef __ATTRIBUTES__

#endif
