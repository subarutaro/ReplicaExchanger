#ifndef HPP_Janus
#define HPP_Janus

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>

#include <cmath>
#include <cassert>

#include "condition.hpp"
#include "quaternion.hpp"
#include <vector.hpp>

double random_number(const double max){
  return max*((double)rand() / (double)RAND_MAX);
}

#define NPATCH_MAX 4
class PatchyInfo{
public:
  int npatch = 1;
  dvec3 patch[NPATCH_MAX];
  double tm[NPATCH_MAX];

  double coef_r;
  double coef_a[NPATCH_MAX][NPATCH_MAX];
  double coef_v;
  double dj = 1.0; // size of solvent particle
  double sol_ratio = 0.0;

  PatchyInfo(){}
  PatchyInfo(const std::string filename){
    std::ifstream ifs(filename);
    {
      std::string line; std::getline(ifs,line);
      while(line[0]=='#');
      std::stringstream strs(line);
      std::string tag,val;
      strs >> tag >> val;
      if(tag != "npatch:"){
	fprintf(stderr,"error: npatch must be defined at first!\n");
	exit(-1);
      }
      npatch = atoi(val.c_str());
    }
    for(std::string line;std::getline(ifs,line);){
      if(line[0] == '#') continue;
      std::string tag,val;
      std::stringstream strs(line);
      strs >> tag;
      if(tag == "patch:"){
	fprintf(stderr,"patch:");
	for(int i=0;i<npatch;i++){
	  strs >> val;
	  patch[i][0] = atof(val.c_str());
	  strs >> val;
	  patch[i][1] = atof(val.c_str());
	  strs >> val;
	  patch[i][2] = atof(val.c_str());
	  fprintf(stderr," %lf %lf %lf ",patch[i][0],patch[i][1],patch[i][2]);
	}
	fprintf(stderr,"\n");
      }else if(tag == "tm:"){
	fprintf(stderr,"tm:");
	for(int i=0;i<npatch;i++){
	  strs >> val;
	  tm[i] = atof(val.c_str()) / 180. * M_PI;
	  fprintf(stderr," %lf",tm[i]);
	}
	fprintf(stderr,"\n");
      }else if(tag == "coef_r:"){
	strs >> val;
	coef_r = atof(val.c_str());
	fprintf(stderr,"coef_r: %lf\n",coef_r);
      }else if(tag == "coef_a:"){
	fprintf(stderr,"coef_a:");
	for(int i=0;i<npatch;i++){
	  for(int j=0;j<npatch;j++){
	    strs >> val;
	    coef_a[i][j] = atof(val.c_str());
	    fprintf(stderr," %lf",coef_a[i][j]);
	  }
	}
	fprintf(stderr,"\n");
      }else if(tag == "coef_v:"){
	strs >> val;
	coef_v = atof(val.c_str());
	fprintf(stderr,"coef_v: %lf\n",coef_v);
      }else if(tag == "sol_ratio:"){
	strs >> val;
	sol_ratio = atof(val.c_str());
	fprintf(stderr,"sol_ratio: %lf\n",sol_ratio);
      }else if(tag == "dj:"){
	strs >> val;
	dj = atof(val.c_str());
	fprintf(stderr,"dj: %lf\n",dj);
      }else{
	fprintf(stderr,"error: undefined variable for PachyInfo!\n");
	exit(-1);
      }
    }
  }
  PatchyInfo(const PatchyInfo& pinfo){
    npatch = pinfo.npatch;
    for(int i=0;i<npatch;i++){
      patch[i] = pinfo.patch[i];
      tm[i] = pinfo.tm[i];
      for(int j=0;j<npatch;j++){
	coef_a[i][j] = pinfo.coef_a[i][j];
      }
    }
    coef_r = pinfo.coef_r;
    coef_v = pinfo.coef_v;
  }
};

template <class TCondition,class TResult>
class Janus{
private:
  PatchyInfo pinfo;
  const double Rwall = 1.0; // width of slit or diameter of tube
  const double density_wall = 10.0; // density of wall
  const double coef_a_wall[2] = {100.,100.}; // attractive coefficient of {patchy, solvent}

public:

  dvec3      *r = nullptr;
  dvec3      *v = nullptr;
  dvec3      *f = nullptr;

  Quaternion *a = nullptr;
  dvec3      *w = nullptr;
  dvec3      *t = nullptr;

  int *type = nullptr;

  double area,vol;
  double hight;
  dvec3  side,sideh;
  dvec3  virial;

  double g;

     int N = 256;
  double T = 1.0;
  double Tprev = 0.0;
  double P = 1.0;

  double Q = 1.0, zeta = 0.0, logs = 0.0;
  double M = 10.0, Pv = 0.0;

  double dt = 0.0005;

  double pot,tot;
  dvec3  kin;
  double rot;
  double Ttmp;
  dvec3  Ptmp;

  void initialize(const int _N,const PatchyInfo& _pinfo){
    N = _N;
    if(r==nullptr) r = new dvec3[N];
    if(v==nullptr) v = new dvec3[N];
    if(f==nullptr) f = new dvec3[N];
    if(a==nullptr) a = new Quaternion[N];
    if(w==nullptr) w = new dvec3[N];
    if(t==nullptr) t = new dvec3[N];

    if(type==nullptr) type = new int[N];

    pinfo = _pinfo;
  }

  void initialize(const std::string filename = "input.txt"){
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
      }else if(opt == "N"){
	N = atoi(val.c_str());
	fprintf(stderr,"N\t%d\n",N);
      }else if(opt == "T"){
	T = atof(val.c_str());
	fprintf(stderr,"T\t%lf\n",T);
      }else if(opt == "P"){
	P = atof(val.c_str());
	fprintf(stderr,"P\t%lf\n",P);
      }else{
	fprintf(stderr,"error: undefined option %s\n",opt.c_str());
	exit(-1);
      }
    }
    if(r == nullptr)r = new dvec3[N];
    if(v == nullptr)v = new dvec3[N];
    if(f == nullptr)f = new dvec3[N];
    if(a == nullptr)a = new Quaternion[N];
    if(w == nullptr)w = new dvec3[N];
    if(t == nullptr)t = new dvec3[N];
    if(type == nullptr)type = new int[N];
  }

  ~Janus(){
    if(r != nullptr) delete[] r;
    if(v != nullptr) delete[] v;
    if(f != nullptr) delete[] f;
    if(a != nullptr) delete[] a;
    if(w != nullptr) delete[] w;
    if(t != nullptr) delete[] t;
    if(type != nullptr) delete[] type;
  }

  void setParticlePlane(const double density = 1.0){
    side = sqrt((double)N/density);
    area  = side[0] * side[1];
    hight = Rwall;
    vol = area*hight;
    side[2] = hight;
    sideh   = side*0.5;
    int nunit = 1;
    while(nunit*nunit < N) nunit++;

    const dvec3 unit_size = side/(double)nunit;
    int ip=0;
    for(int i=0; i<nunit; i++){
      for(int j=0; j<nunit; j++){
	if(ip == N) continue;
	r[ip][0] = (i*unit_size[0])/side[0];
	r[ip][1] = (j*unit_size[1])/side[1];
	r[ip][2] = 0.0;
	v[ip][0] = random_number(2.0) - 1.0;
	v[ip][1] = random_number(2.0) - 1.0;
	v[ip][2] = random_number(2.0) - 1.0;

	a[ip][0] = 1.0;
	a[ip][1] = 0.0;
	a[ip][2] = 0.0;
	a[ip][3] = 0.0;
	w[ip][0] = random_number(2.0) - 1.0;
	w[ip][1] = random_number(2.0) - 1.0;
	w[ip][2] = random_number(2.0) - 1.0;

	type[ip] = 0;

	g += 6.0;
	ip++;
      }
    }
    assert(ip == N);

    const int Nsol = N * pinfo.sol_ratio;
    for(int count=0;count<Nsol;){
      const int pivot = rand()%N;
      if(type[pivot] == 0){
	type[pivot] = 1;
	g -= 3.0;
	count++;
      }
    }

    for(int i=0; i<N; i++){
      r[i][0] -= 0.5;
      r[i][1] -= 0.5;

      if(r[i][0] >= 0.5) r[i][0] -= 1.0;
      if(r[i][1] >= 0.5) r[i][1] -= 1.0;
      if(r[i][0] < -0.5) r[i][0] += 1.0;
      if(r[i][1] < -0.5) r[i][1] += 1.0;

      assert(-0.5 <= r[i][0] && r[i][0] < 0.5);
      assert(-0.5 <= r[i][1] && r[i][1] < 0.5);
    }

    dvec3  cm_vel = 0.0;
    for(int i=0; i<N; i++){
      cm_vel  += v[i];
    }
    cm_vel /= (double)N;
    for(int i=0; i<N; i++){
      v[i] -= cm_vel;
    }
  }

  void adjustVolume(){
    calcForce();
    calcKineticEnergy();
    calcPressure();
    const double acceptable_error = 0.01;
    double diff = 0.5*(Ptmp[0]+Ptmp[1]) - P;
    while(fabs(diff) > acceptable_error){
      if(diff>0) updateVolume(area*1.01);
      else       updateVolume(area*0.99);
      calcForce();
      calcPressure();
      diff = 0.5*(Ptmp[0]+Ptmp[1]) - P;
      //printf("diff: %lf %lf %lf %lf\n",area,0.5*(Ptmp[0]+Ptmp[1]),diff,vol);
    }
    #ifdef FMRES_DEBUG_PRINT
    printf("area is adjusted: A= %lf, Ptmp= %lf, P=  %lf\n",area,0.5*(Ptmp[0]+Ptmp[1]),P);
    #endif
  }

  double calcKineticEnergy(){
    calcTranslationalKineticEnergy();
    calcRotationalKineticEnergy();
    return sum(kin) + rot;
  }

  double calcTranslationalKineticEnergy(){
    kin = 0.0;
    for(int i=0;i<N;i++){
      const dvec3 tmp = v[i] / side;
      kin += tmp*tmp;
    }
    kin *= 0.5;
    return sum(kin);
  }
  double calcRotationalKineticEnergy(){
    rot = 0.0;
    for(int i=0;i<N;i++){
      if(type[i]!=1) rot += w[i]%w[i];
    }
    rot *= 0.5;
    return rot;
  }

  void scaleVelocity(){
    const double vscale = sqrt(T / (sum(kin) / (1.5*N)));
    const double wscale = sqrt(T / (rot / (0.5*(g-3*N))));
    //printf("scale: %lf %lf\n",vscale,wscale);
    for(int i=0;i<N;i++){
      v[i] *= vscale;
      if(type[i] == 0) w[i] *= wscale;
    }

    Pv   = 0.0;
    zeta = 0.0;
    logs = 0.0;
  }

  void scaleVelocityAfterExchange(){
    if(T == Tprev || Tprev == 0.0) return;
    const double scale = sqrt(T / Tprev);
    //printf("scale: %lf %lf\n",vscale,wscale);
    for(int i=0;i<N;i++){
      v[i] *= scale;
      if(type[i] == 0) w[i] *= scale;
    }
    Pv   = 0.0;
    zeta = 0.0;
    logs = 0.0;
  }

  void setRandomVelocity(){
    for(int i=0;i<N;i++){
      v[i][0] = random_number(1.0);
      v[i][1] = random_number(1.0);
      v[i][2] = random_number(1.0);
      w[i][0] = random_number(1.0);
      w[i][1] = random_number(1.0);
      w[i][2] = random_number(1.0);
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
    calcKineticEnergy();
    scaleVelocity();
  }

  void calcForce(){
    const double ph = M_PI*0.5;
    for(int i=0;i<N;i++){
      f[i] = 0.0;
      t[i] = 0.0;
    }
    virial = 0.0;

#pragma omp parallel for
    for(int i=0;i<N;i++){
      const dvec3 ri = r[i];
      const Quaternion ai = a[i];
      const int type_i = type[i];
      dvec3 fi = 0.0;
      double u = 0.0;
      const double dj = pinfo.dj;
      const double dij[2][2] = {{1.0,(1.0+dj)*0.5},{(1.0+dj)*0.5,dj}};
      for(int j=0;j<N;j++){
	if(i==j) continue;
	dvec3 dr = ri - r[j];
	if(dr[0] >= 0.5) dr[0] -= 1.0;
	if(dr[0] < -0.5) dr[0] += 1.0;
	if(dr[1] >= 0.5) dr[1] -= 1.0;
	if(dr[1] < -0.5) dr[1] += 1.0;
	dr *= side;
	const double r2 = dr%dr;
	const double d = dij[type_i][type[j]];
	if(r2 > d) continue;
	const double rinv = 1.0/sqrt(r2);
	const double r2i = rinv*rinv;
	const double rd  = r2*rinv/d;
	const double wij = 1.0 - rd;
	// repulsive force
	dvec3 ftmp = dr*(pinfo.coef_r * wij * rinv);
	dvec3 ttmp = 0.0;
	#if 1
	if(type_i == 0 && type[j] == 0){
	  // attractive force and torque
	  const Quaternion  aj = a[j];
	  for(int k=0;k<pinfo.npatch;k++){
	    const dvec3 ni = Rotate(ai)*pinfo.patch[k];
	    const double cos_tm_i = cos(pinfo.tm[k]);
	    for(int l=0;l<pinfo.npatch;l++){
	      dvec3 nj = Rotate(aj)*pinfo.patch[l];
	      const double cos_ti = -(ni % dr) * rinv;
	      const double cos_tj =  (nj % dr) * rinv;
	      const double cos_tm_j = cos(pinfo.tm[l]);

	      if(cos_ti < cos_tm_i || cos_tj < cos_tm_j) continue;
	      const double ti = acos(cos_ti);
	      const double tj = acos(cos_tj);

	      const double cos_phtitm = cos(ph*ti/pinfo.tm[k]);
	      const double cos_phtjtm = cos(ph*tj/pinfo.tm[l]);

	      const double ff = cos_phtitm * cos_phtjtm;
	      const double fv  = (ff==0.0) ? 0.0 : powf(ff, pinfo.coef_v);
	      const double fvi = (ff==0.0) ? 0.0 : powf(ff, pinfo.coef_v - 1.0);

	      const dvec3 dcostidr = -(ni - dr*((ni%dr)*r2i))*rinv;
	      const dvec3 dcostjdr =  (nj - dr*((nj%dr)*r2i))*rinv;

	      double dtidcossin;
	      if(cos_ti*cos_ti != 1.0)
		dtidcossin = - sin(ph*ti/pinfo.tm[k]) / sqrt(1.0 - cos_ti*cos_ti);
	      else
		dtidcossin = - ph/pinfo.tm[k];

	      double dtjdcossin;
	      if(cos_tj*cos_tj != 1.0)
		dtjdcossin = - sin(ph*tj/pinfo.tm[l]) / sqrt(1.0 - cos_tj*cos_tj);
	      else
		dtjdcossin = - ph/pinfo.tm[l];

	      const double tmpi = dtidcossin * cos_phtjtm;
	      const double tmpj = dtjdcossin * cos_phtitm;

	      ftmp += dr*(pinfo.coef_a[k][l] * fv * (0.5 - rd)*rinv);
	      ftmp -= (dcostidr*(tmpi/pinfo.tm[k]) + dcostjdr*(tmpj/pinfo.tm[l]))*(0.5*pinfo.coef_a[k][l] * rd*(1.0 - rd) * pinfo.coef_v * fvi * ph);
	      const double ttmp2 = 0.5*ph*pinfo.coef_a[k][l]/pinfo.tm[k]*rd*(1.0 - rd)*pinfo.coef_v*fvi*tmpi;
	      const dvec3 g = dr*ttmp2;
	      t[i] += ni ^ g;
	      //ttmp += nj ^ g;
	    }
	  }
	}
	#endif
	fi += ftmp;
	/*
#pragma omp atomic
	{
	  f[j] -= ftmp;
	  t[j] += ttmp;
	}
	//*/
	virial += ftmp * dr;
      } // j loop
      // wall force
      const double Rup = 0.5*Rwall - ri[2];
      const double Rdw = 0.5*Rwall + ri[2];
      if(Rup <= 1.0)
	fi[2] -= 0.16666666666 * M_PI * density_wall * coef_a_wall[type[i]]
	  * (1.0 - 2.0 * Rup + 2.0 * Rup*Rup*Rup - Rup*Rup*Rup*Rup);
      if(Rdw <= 1.0)
	fi[2] += 0.16666666666 * M_PI * density_wall * coef_a_wall[type[i]]
	  * (1.0 - 2.0 * Rdw + 2.0 * Rdw*Rdw*Rdw - Rdw*Rdw*Rdw*Rdw);

      f[i] += fi;
    }
    virial *= 0.5;
  }

  double calcPotentialEnergy(){
    const double ph = 0.5*M_PI;
    pot = 0.0;
#pragma omp parallel for reduction(+:pot)
    for(int i=0;i<N;i++){
      const dvec3 ri = r[i];
      const Quaternion ai = a[i];
      const int type_i = type[i];
      const double dj = pinfo.dj;
      const double dij[2][2] = {{1.0,(1.0+dj)*0.5},{(1.0+dj)*0.5,dj}};
      for(int j=i+1;j<N;j++){
	dvec3 dr = ri - r[j];
	if(dr[0] >= 0.5) dr[0] -= 1.0;
	if(dr[0] < -0.5) dr[0] += 1.0;
	if(dr[1] >= 0.5) dr[1] -= 1.0;
	if(dr[1] < -0.5) dr[1] += 1.0;
	dr *= side;
	const double r2 = dr%dr;
	const double d = dij[type_i][type[j]];
	if(r2 > d) continue;
	const double rinv = 1.0/sqrt(r2);
	const double r2i = rinv*rinv;
	const double rd  = r2*rinv/d;
	const double wij = 1.0 - rd;
	pot += 0.5 * pinfo.coef_r * wij*wij;

	if(type_i == 0 && type[j] == 0){
	  // attractive force and torque
	  const Quaternion  aj = a[j];
	  for(int k=0;k<pinfo.npatch;k++){
	    const dvec3 ni = Rotate(ai)*pinfo.patch[k];
	    const double cos_tm_i = cos(pinfo.tm[k]);
	    for(int l=0;l<pinfo.npatch;l++){
	      dvec3 nj = Rotate(aj)*pinfo.patch[l];
	      const double cos_ti = -(ni % dr) * rinv;
	      const double cos_tj =  (nj % dr) * rinv;
	      const double cos_tm_j = cos(pinfo.tm[l]);

	      if(cos_ti < cos_tm_i || cos_tj < cos_tm_j) continue;
	      const double ti = acos(cos_ti);
	      const double tj = acos(cos_tj);

	      const double cos_phtitm = cos(ph*ti/pinfo.tm[k]);
	      const double cos_phtjtm = cos(ph*tj/pinfo.tm[l]);

	      const double ff = cos_phtitm * cos_phtjtm;
	      const double fv  = (ff==0.0) ? 0.0 : powf(ff, pinfo.coef_v);
	      const double fvi = (ff==0.0) ? 0.0 : powf(ff, pinfo.coef_v - 1.0);

	      const dvec3 dcostidr = -(ni - dr*((ni%dr)*r2i))*rinv;
	      const dvec3 dcostjdr =  (nj - dr*((nj%dr)*r2i))*rinv;

	      double dtidcossin;
	      if(cos_ti*cos_ti != 1.0)
		dtidcossin = - sin(ph*ti/pinfo.tm[k]) / sqrt(1.0 - cos_ti*cos_ti);
	      else
		dtidcossin = - ph/pinfo.tm[k];

	      double dtjdcossin;
	      if(cos_tj*cos_tj != 1.0)
		dtjdcossin = - sin(ph*tj/pinfo.tm[l]) / sqrt(1.0 - cos_tj*cos_tj);
	      else
		dtjdcossin = - ph/pinfo.tm[l];

	      const double tmpi = dtidcossin * cos_phtjtm;
	      const double tmpj = dtjdcossin * cos_phtitm;

	      pot -= 0.5*fv*pinfo.coef_a[k][l]*rd*(1.0 - rd);
	    }
	  }
	}
      } // j loop
    }
    return pot;
  }
  void calcEnergy(){
    calcPotentialEnergy();
    calcKineticEnergy();
    tot  = pot+sum(kin);
    tot += 0.5*zeta*zeta*Q + 3.0*N*T*logs;
    tot += 0.5*Pv*Pv/M + P*vol;
  }

  void calcPressure(){
    calcKineticEnergy();
    //printf("press: %lf %lf %lf\n",sum(kin),rot,sum(virial));
    Ptmp = (virial + kin*2.0) / vol;
  }

  void updateVolume(const double _area){
    area  = _area;
    vol   = area * hight;
    side  = sqrt(area);
    side[2] = hight;
    sideh = side*0.5;
  }

  void integrateOneStep(){
    const double dth = 0.5*dt;
    for(int i=0;i<N;i++){
      v[i] *= exp(-zeta*dth);
      v[i] += f[i]*side*dth;
      if(type[i] == 0){
	w[i] *= exp(-zeta*dth);
	w[i] += t[i] * dth;
      }
    }
    Pv += (0.5*(virial[0]+virial[1])/area - P*hight)*dth;

    updateVolume(area + Pv/M*dt);

    for(int i=0;i<N;i++){
      r[i] += (v[i]/(side*side))*dt;
      if(r[i][0] >= 0.5) r[i][0] -= 1.0;
      if(r[i][0] < -0.5) r[i][0] += 1.0;
      if(r[i][1] >= 0.5) r[i][1] -= 1.0;
      if(r[i][1] < -0.5) r[i][1] += 1.0;
#if 1
      if(type[i] == 0){
	//printf("a %d: %lf %lf %lf %lf\n",i,a[i][0],a[i][1],a[i][2],a[i][3]);
	//printf("w %d: %lf %lf %lf\n",i,w[i][0],w[i][1],w[i][2]);
	//printf("t %d: %lf %lf %lf\n",i,t[i][0],t[i][1],t[i][2]);
	const Quaternion w1  = S(a[i]) * Quaternion(0.0, w[i]) * 0.5;
	const Quaternion a1  = normalize(a[i] + w1 * dth * 2.0);
	const Quaternion a1h = normalize(a[i] + w1 * dth);
	const Quaternion w2  = S(a1h) * Quaternion(0.0, w[i]) * 0.5;
	const Quaternion a2  = normalize(a1h + w2*dth);
	//printf("a1: %lf %lf %lf %lf\n",a1[0],a1[1],a1[2],a1[3]);
	//printf("a1h: %lf %lf %lf %lf\n",a1h[0],a1h[1],a1h[2],a1h[3]);
	//printf("a2: %lf %lf %lf %lf\n",a2[0],a2[1],a2[2],a2[3]);
	a[i] = normalize(a2 * 2.0 - a1);
	assert(fabs(a[i]%a[i] - 1.0)<1e-6);
      }
#endif
    }
    calcKineticEnergy();
    zeta += (2.0*(sum(kin)+rot) - g*T)/Q * dt;
    logs += zeta*dt;
    Pv   += (kin[0]+kin[1])/vol*dt;

    updateVolume(area + Pv/M*dt);
    calcForce();

    Pv += (0.5*(virial[0]+virial[1])/area - P*hight)*dth;
    for(int i=0;i<N;i++){
      v[i] += f[i]*side*dth;
      v[i] *= exp(-zeta*dth);
      if(type[i] == 0){
	w[i] += t[i] * dth;
	w[i] *= exp(-zeta*dth);
      }
    }
  }

  void setPressure(const double _P){P=_P;};
  void setTemperature(const double _T){Tprev=T;T=_T;};
  void setDeltaTime(const double _dt){dt=_dt;};

  double getPotentialEnergy()    const {return pot;}
  double getKineticEnergy()      const {return sum(kin);}
  double getTotalEnergy()        const {return tot;}
  double getTargetTemperature()  const {return T;}
  double getTargetPressure()     const {return P;}
  double getCurrentTemperature() const {return (sum(kin)+rot) / (0.5*g);}
  double getCurrentPressure()    const {return 0.5*(Ptmp[0]+Ptmp[1]);}
  double getVolume()             const {return vol;};
  TResult getResult() const {
    TResult ret;
    ret.pot = pot;
    ret.kin = sum(kin)+rot;
    ret.tot = tot;
    ret.vol = vol;

    ret.T = getCurrentTemperature();
    ret.P = getCurrentPressure();
    return ret;
  }
  TCondition getCondition() const {
    TCondition ret;
    ret.Beta  = 1.0/T;
    ret.PBeta = P/T;
    ret.T = T;
    ret.P = P;
    return ret;
  }
  void setCondition(TCondition& c){
    setTemperature(c.T);
    setPressure(c.P);
  }

  void outputCDV(std::string filename) const {
    FILE *fp;
    if((fp = fopen(filename.c_str(),"w"))==NULL){
      fprintf(stderr,"error: opening %s failed\n",filename.c_str());
      exit(-1);
    }
    fprintf(fp,"'box_sx= %lf,box_sy= %lf,box_sz= %lf,box_ex= %lf,box_ey= %lf,box_ez= %lf\n",-sideh[0],-sideh[1],-sideh[2],sideh[0],sideh[1],sideh[2]);
    for(int i=0;i<N;i++){
      dvec3 rtmp = r[i] * side;
      fprintf(fp,"%d %d %lf %lf %lf\n",i,type[i],rtmp[0],rtmp[1],rtmp[2]);
    }
    fclose(fp);
  }
  void outputCheckPoint(std::string filename) const {
    FILE *fp;
    if((fp = fopen(filename.c_str(),"wb"))==NULL){
      fprintf(stderr,"error: opening %s failed\n",filename.c_str());
      exit(-1);
    }
    fwrite(&N,sizeof(int),1,fp);
    fwrite(&P,sizeof(double),1,fp);
    fwrite(&T,sizeof(double),1,fp);

    fwrite(r,sizeof(dvec3),N,fp);
    fwrite(v,sizeof(dvec3),N,fp);
    fwrite(f,sizeof(dvec3),N,fp);
    fwrite(a,sizeof(Quaternion),N,fp);
    fwrite(w,sizeof(dvec3),N,fp);
    fwrite(t,sizeof(dvec3),N,fp);
    fwrite(type,sizeof(int),N,fp);

    fwrite(&area,sizeof(double),1,fp);
    fwrite(&vol,sizeof(double),1,fp);
    fwrite(&hight,sizeof(double),1,fp);
    fwrite(&side,sizeof(dvec3),1,fp);
    fwrite(&sideh,sizeof(dvec3),1,fp);
    fwrite(&virial,sizeof(dvec3),1,fp);
    fwrite(&g,sizeof(double),1,fp);
    fwrite(&Q,sizeof(double),1,fp);
    fwrite(&zeta,sizeof(double),1,fp);
    fwrite(&logs,sizeof(double),1,fp);
    fwrite(&M,sizeof(double),1,fp);
    fwrite(&Pv,sizeof(double),1,fp);
    fwrite(&dt,sizeof(double),1,fp);
    fwrite(&pot,sizeof(double),1,fp);
    fwrite(&tot,sizeof(double),1,fp);
    fwrite(&kin,sizeof(dvec3),1,fp);
    fwrite(&rot,sizeof(double),1,fp);
    fwrite(&Ttmp,sizeof(double),1,fp);
    fwrite(&Ptmp,sizeof(dvec3),1,fp);
    fwrite(&pinfo,sizeof(PatchyInfo),1,fp);

    fclose(fp);
  }
  void readCheckPoint(std::string filename){
    FILE *fp;
    if((fp = fopen(filename.c_str(),"rb"))==NULL){
      fprintf(stderr,"error: opening %s failed\n",filename.c_str());
      exit(-1);
    }

    int Ntmp;
    fread(&Ntmp,sizeof(int),1,fp);
    if(N != Ntmp){
      fprintf(stderr,"error: N(=%d) != Ntmp(%d)\n",N,Ntmp);
      exit(EXIT_FAILURE);
    }
    fread(&P,sizeof(double),1,fp);
    fread(&T,sizeof(double),1,fp);

    assert(r != nullptr);
    assert(v != nullptr);
    assert(f != nullptr);
    assert(a != nullptr);
    assert(w != nullptr);
    assert(t != nullptr);
    assert(type != nullptr);
    fread(r,sizeof(dvec3),N,fp);
    fread(v,sizeof(dvec3),N,fp);
    fread(f,sizeof(dvec3),N,fp);
    fread(a,sizeof(Quaternion),N,fp);
    fread(w,sizeof(dvec3),N,fp);
    fread(t,sizeof(dvec3),N,fp);
    fread(type,sizeof(int),N,fp);

    fread(&area,sizeof(double),1,fp);
    fread(&vol,sizeof(double),1,fp);
    fread(&hight,sizeof(double),1,fp);
    fread(&side,sizeof(dvec3),1,fp);
    fread(&sideh,sizeof(dvec3),1,fp);
    fread(&virial,sizeof(dvec3),1,fp);
    fread(&g,sizeof(double),1,fp);
    fread(&Q,sizeof(double),1,fp);
    fread(&zeta,sizeof(double),1,fp);
    fread(&logs,sizeof(double),1,fp);
    fread(&M,sizeof(double),1,fp);
    fread(&Pv,sizeof(double),1,fp);
    fread(&dt,sizeof(double),1,fp);
    fread(&pot,sizeof(double),1,fp);
    fread(&tot,sizeof(double),1,fp);
    fread(&kin,sizeof(dvec3),1,fp);
    fread(&rot,sizeof(double),1,fp);
    fread(&Ttmp,sizeof(double),1,fp);
    fread(&Ptmp,sizeof(dvec3),1,fp);
    fread(&pinfo,sizeof(PatchyInfo),1,fp);

    fclose(fp);
  }
};

#ifdef Janus_TEST
int main(int argc,char **argv){
  Janus<Condition,Result> janus;
  janus.initialize("input.txt");
  janus.setCubicFCC();
  janus.setRandomVelocityAtTargetTemperature();
  janus.calcForce();

  const int nstep    = 2000;
  const int nstep_eq = 2000;
  for(int s=-nstep_eq;s<nstep;s++){
    if(s<0){
      janus.scaleVelocity();
    }
    janus.integrateOneStep();
    if(s%10 == 0){
      janus.calcEnergy();
      Result r = janus.getResult();
      fprintf(stdout,"%d %lf %lf %lf %lf %lf\n",
	      s,r.pot,r.kin,r.tot,r.T,r.P);
    }
  }
}
#endif

#endif
