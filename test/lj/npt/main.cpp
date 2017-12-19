#include <iostream>

#include "lj.hpp"

#include <vector.hpp>
#include <replica_exchanger.hpp>

int main(int argc,char **argv){
  RE::Initialize(argc,argv);

  RE::ReplicaExchanger<2,LJ<Condition,Result>,Condition,Result> re;
  Vector<int,2> ndim;
  ndim[0] = 16;
  ndim[1] = 10;
  re.Initialize(ndim);
  for(int i=0;i<re.getLocalNumberOfReplicas();i++){
    re[i].initialize(256);
    re[i].setCubicFCC();
  }
  Vector<double,2> max,min;
  max[0] = 1.0; max[1] = 1.0;
  min[0] = 0.5; min[1] = 0.5;
  re.generateConditionRegularInterval(max,min);
  re.copyREMInfoToReplicas();
  for(int i=0;i<re.getLocalNumberOfReplicas();i++){
    re[i].setRandomVelocityAtTargetTemperature();
    re[i].calcForce();
    re[i].calcEnergy();
    re[i].calcPressure();
    re[i].adjustVolume();
  }
  const int ninterval = 10000;
  for(int s=-1000;s<ninterval;s++){
    printf("step %d\n",s);
    re.copyREMInfoToReplicas();
    if(s>=0) re.exchangeReplicasWithNeighbor();
    for(int i=0;i<re.getLocalNumberOfReplicas();i++){
      for(int ss=0;ss<100;ss++){
	if(s<0){
	  re[i].scaleVelocity();
	  re[i].adjustVolume();
	}
	re[i].integrateOneStep();
      }
      re[i].calcEnergy();
      re[i].calcPressure();
    }
    re.getREMInfoFromReplicas();
    if(s>=0) re.Output();
  }

  RE::Finalize();
}
