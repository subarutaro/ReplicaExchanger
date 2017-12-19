#ifndef H_CONDITION
#define H_CONDITION

#include <iostream>
#include <sstream>
#include <string>
#include "vector.hpp"

class Condition{
public:
  double beta,T;
  const Condition& operator=(const Condition& c){
    beta = c.beta;
    T = 1.0/beta;
    return c;
  }
  void set(Vector<double,1> c){
    beta = c[0];
    T = 1.0/beta;
  }
  const std::string filename() const {
    std::stringstream strs;
    strs.setf(std::ios::fixed);
    strs  << "T" << std::setprecision(6) << T;
    strs << ".dat";
    return strs.str();
  }

  friend std::ostream& operator<<(std::ostream& os,const Condition& c){
    os << c.T;
    return os;
  }
};

class Result{
public:
  double pot,kin,tot;
  double T,P;

  double getLogBoltzmannFactor(const Condition& c){
    //return -pot/c.T;
    return 0.0;
  }
  const Result operator=(const Result& r){
    pot = r.pot;
    kin = r.kin;
    tot = r.tot;
    T = r.T;
    P = r.P;
    return r;
  }
  friend std::ostream& operator<<(std::ostream& os,const Result& r){
    os << r.pot << " " << r.kin << " " << r.tot << " " << r.T << " " << r.P;
    return os;
  }
};

#endif
