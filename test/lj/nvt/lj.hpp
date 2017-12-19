#ifndef HPP_LJ
#define HPP_LJ

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>

#include <cmath>
#include <cassert>

#include "condition.hpp"
#include <vector.hpp>

double random_number(const double max){
  return max*((double)rand() / (double)RAND_MAX - 0.5);
}

template <class TCondition,class TResult>
class LJ{
private:
  bool isInitialized = false;
public:
  dvec3 *r = nullptr;
  dvec3 *v = nullptr;
  dvec3 *f = nullptr;

  double volume;
  dvec3  side,sideh;
  dvec3  virial;

  int    N = 256;
  double T = 1.0;
  double density = 1.0;

  double M = 1.0;
  double dt = 0.0005;
  double rcut = 4.5;


  double pot,tot;
  dvec3  kin;
  double Ttmp;
  dvec3  Ptmp;

  void initialize(const int _N){
    if(isInitialized){
      fprintf(stderr,"error: this object is already initialized!\n");
      return;
    }
    N = _N;
    if(r == nullptr) r = new dvec3[N];
    if(v == nullptr) v = new dvec3[N];
    if(f == nullptr) f = new dvec3[N];
    isInitialized = true;
  }

  void initialize(const std::string filename = "input.txt"){
    if(isInitialized){
      fprintf(stderr,"error: this object is already initialized!\n");
      return;
    }
    std::ifstream ifs(filename);
    assert(!ifs.fail());
    for(std::string line;std::getline(ifs,line);){
      std::stringstream strs(line);
      std::string opt,val;
      strs >> opt >> val;
      if(opt[0] == '#') continue;
      if(opt == "dt"){
	dt = atof(val.c_str());
	fprintf(stderr,"dt\t%lf\n",dt);
      }else if(opt == "rcut"){
	rcut = atof(val.c_str());
	fprintf(stderr,"rcut\t%lf\n",rcut);
      }else if(opt == "N"){
	N = atoi(val.c_str());
	fprintf(stderr,"N\t%d\n",N);
      }else if(opt == "T"){
	T = atof(val.c_str());
	fprintf(stderr,"T\t%lf\n",T);
      }else if(opt == "density"){
	density = atof(val.c_str());
	fprintf(stderr,"density\t%lf\n",density);
      }else{
	fprintf(stderr,"error: undefined option %s\n",opt.c_str());
	exit(-1);
      }
    }
    if(r == nullptr) r = new dvec3[N];
    if(v == nullptr) v = new dvec3[N];
    if(f == nullptr) f = new dvec3[N];
    isInitialized = true;
  }

  ~LJ(){
    if(r != nullptr) delete[] r;
    if(v != nullptr) delete[] v;
    if(f != nullptr) delete[] f;
  }

  void setCubicFCC(){
    int nside = 1;
    while(4*nside*nside*nside < N) nside++;
    side   = powf((double)N/density,1./3.);
    volume = side[0]*side[1]*side[2];
    sideh  = side * 0.5;
    if(min(sideh) > rcut) rcut = min(sideh);

    const dvec3  unit  = side / (double)nside;
    const dvec3  unith = unit*0.5;
    const double units[4][3] = {
      {     0.0,      0.0,     0.0},
      {unith[0], unith[1],     0.0},
      {     0.0, unith[1], unith[2]},
      {unith[0],      0.0, unith[2]}
    };
    int count = 0;
    for(int i=0;i<nside;i++)
      for(int j=0;j<nside;j++)
	for(int k=0;k<nside;k++)
	  for(int l=0;l<4;l++){
	    if(count == N) break;
	    r[count][0] = i*unit[0] + units[l][0] - sideh[0];
	    r[count][1] = j*unit[1] + units[l][1] - sideh[1];
	    r[count][2] = k*unit[2] + units[l][2] - sideh[2];
	    count++;
	}
  }


  double calcKineticEnergy(){
    kin = 0.0;
    for(int i=0;i<N;i++){
      kin += v[i]*v[i];
    }
    kin *= 0.5;
    return sum(kin);
  }

  void scaleVelocity(){
    const double kin = calcKineticEnergy();
    const double scale = sqrt(T / (kin / (1.5*N)));
    for(int i=0;i<N;i++) v[i] *= scale;
  }

  void setRandomVelocity(){
    for(int i=0;i<N;i++){
      v[i][0] = random_number(1.0);
      v[i][1] = random_number(1.0);
      v[i][2] = random_number(1.0);
    }
    removeTotalMomentum();
  }

  void removeTotalMomentum(){
    dvec3 mom = 0.0;
    for(int i=0;i<N;i++) mom += v[i];
    mom /= (double)N;
    for(int i=0;i<N;i++) v[i] -= mom;
  }

  void setRandomVelocityAtTargetTemperature(){
    setRandomVelocity();
    scaleVelocity();
  }

  void calcForce(){
    const double rcut2 = rcut*rcut;
    for(int i=0;i<N;i++) f[i] = 0.0;
    virial = 0.0;

#pragma omp parallel for
    for(int i=0;i<N;i++){
      const dvec3 ri = r[i];
      const dvec3 vi = v[i];
      dvec3 fi = 0.0;
      for(int j=i+1;j<N;j++){
	dvec3 rij = ri - r[j];
	for(int d=0;d<3;d++){
	  if(rij[d] <= -sideh[d]) rij[d] += side[d];
	  if(rij[d] >   sideh[d]) rij[d] -= side[d];
	}
	const double r2 = sum(rij*rij);
	if(r2 > rcut2) continue;
	const double r02i = 1.0/r2;
	const double r06i = r02i*r02i*r02i;
	const double r12i = r06i*r06i;
	const double fc = 24.*(2. * r12i - r06i) * r02i;
	const dvec3 ftmp = rij * fc;
	fi += ftmp;
	#pragma omp atomic
	f[j] -= ftmp;

	virial[0] += ftmp[0] * rij[0];
	virial[1] += ftmp[1] * rij[1];
	virial[2] += ftmp[2] * rij[2];
      }
      f[i] += fi;
    }
  }

  double calcPotentialEnergy(){
    const double rcut2 = rcut*rcut;
    pot = 0.0;
#pragma omp parallel for
    for(int i=0;i<N;i++){
      const dvec3 ri = r[i];
      for(int j=i+1;j<N;j++){
	dvec3 rij = ri - r[j];
	for(int d=0;d<3;d++){
	  if(rij[d] <= -sideh[d]) rij[d] += side[d];
	  if(rij[d] >   sideh[d]) rij[d] -= side[d];
	}
	const double r2 = sum(rij*rij);
	if(r2 > rcut2) continue;
	const double r2i = 1.0 / r2;
	const double r6i = r2i*r2i*r2i;
	pot += 4.0*r6i*(r6i - 1.0);
      }
    }
    return pot;
  }
  void calcEnergy(){
    calcPotentialEnergy();
    calcKineticEnergy();
    tot = pot+sum(kin);
  }

  double calcPressure(){
    calcKineticEnergy();
    Ptmp = (virial + kin) / (3.0*volume);
    return sum(Ptmp);
  }

  void integrateOneStep(){
    const double dth = 0.5*dt;
    for(int i=0;i<N;i++){
      v[i] += f[i]*dth;
      r[i] += v[i]*dt;
      for(int d=0;d<3;d++){
	if(r[i][d] <= -sideh[d]) r[i][d] += side[d];
	if(r[i][d] >   sideh[d]) r[i][d] -= side[d];
      }
    }
    calcForce();
    for(int i=0;i<N;i++){
      v[i] += f[i]*dth;
    }
  }

  void setTemperature(const double _T){T=_T;};
  void setDensity(const double _density){density=_density;};
  void setDeltaTime(const double _dt){dt=_dt;};
  void setCutoffRadius(const double _rcut){rcut=_rcut;};

  double getPotentialEnergy() const {return pot;}
  double getKineticEnergy() const {return sum(kin);}
  double getTotalEnergy() const {return tot;}
  double getTargetTemperature() const {return T;}
  double getCurrentTemperature() const {return sum(kin) / (1.5*N);}
  double getCurrentPressure() const {return sum(Ptmp);}
  double getVolume()const {return volume;};
  TResult getResult() const {
    TResult ret;
    ret.pot = pot;
    ret.kin = sum(kin);
    ret.tot = tot;

    ret.T = getCurrentTemperature();
    ret.P = getCurrentPressure();
    return ret;
  }
  TCondition getCondition() const {
    TCondition ret;
    ret.beta = 1.0/T;
    ret.T = T;
    return ret;
  }
  void setCondition(TCondition& c){
    setTemperature(c.T);
  }

  void outputCDV(std::string filename){
    FILE *fp;
    if((fp = fopen(filename.c_str(),"w"))==NULL){
      fprintf(stderr,"error: opening %s failed\n",filename.c_str());
      exit(-1);
    }
    fprintf(fp,"'box_sx= %lf,box_sy= %lf,box_sz= %lf,box_ex= %lf,box_ey= %lf,box_ez= %lf\n",-sideh[0],-sideh[1],-sideh[2],sideh[0],sideh[1],sideh[2]);
    for(int i=0;i<N;i++){
      fprintf(fp,"%d 0 %lf %lf %lf\n",i,r[i][0],r[i][1],r[i][2]);
    }
  }
};

#ifdef LJ_TEST
int main(int argc,char **argv){
  LJ<Condition,Result> lj;
  lj.initialize("input.txt");
  lj.setCubicFCC();
  lj.setRandomVelocityAtTargetTemperature();
  lj.calcForce();

  const int nstep    = 2000;
  const int nstep_eq = 2000;
  for(int s=-nstep_eq;s<nstep;s++){
    if(s<0){
      lj.scaleVelocity();
    }
    lj.integrateOneStep();
    if(s%10 == 0){
      lj.calcEnergy();
      Result r = lj.getResult();
      fprintf(stdout,"%d %lf %lf %lf %lf %lf\n",
	      s,r.pot,r.kin,r.tot,r.T,r.P);
    }
  }
}
#endif

#endif
 
