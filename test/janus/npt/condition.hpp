#ifndef H_CONDITION
#define H_CONDITION

#include <iostream>
#include <sstream>
#include <string>
#include "vector.hpp"

class Condition{
public:
  double T,Beta;
  double P,PBeta;
  const Condition& operator=(const Condition& c){
    Beta  = c.Beta;
    PBeta = c.PBeta;
    T = 1.0/c.Beta;
    P = c.PBeta / c.Beta;
    return *this;
  }
  void set(Vector<double,2> c){
    Beta  = c[0];
    PBeta = c[1];
    T = 1.0/Beta;
    P = PBeta / Beta;
  }
  const std::string getPrefix() const {
    std::stringstream strs;
    strs.setf(std::ios::fixed);
    strs  << "P" << std::setprecision(6) << P;
    strs  << "T" << std::setprecision(6) << T;
    return strs.str();
  }

  friend std::ostream& operator<<(std::ostream& os,const Condition& c){
    os << c.P << " " << c.T;
    return os;
  }
};

class Result{
public:
  double pot,kin,tot;
  double vol;
  double T,P;

  double getLogBoltzmannFactor(const Condition& c){
    return -(pot*c.Beta + vol*c.PBeta);
  }
  const Result operator=(const Result& r){
    pot = r.pot;
    kin = r.kin;
    tot = r.tot;
    vol = r.vol;
    T   = r.T;
    P   = r.P;
    return *this;
  }
  friend std::ostream& operator<<(std::ostream& os,const Result& r){
    os << r.pot << " " << r.kin << " " << r.vol << " " << r.T << " " << r.P;
    return os;
  }
};

#endif
