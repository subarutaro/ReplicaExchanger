#include <iostream>

#include "lj.hpp"

#include <vector.hpp>
#include <replica_exchanger.hpp>

int main(int argc,char **argv){
  RE::Initialize(argc,argv);

  RE::ReplicaExchanger<1,LJ<Condition,Result>,Condition,Result> re;
  Vector<int,1> ndim;
  ndim[0] = 36;
  re.Initialize(ndim);
  for(int i=0;i<re.getLocalNumberOfReplicas();i++){
    re[i].initialize(108);
    re[i].setCubicFCC();
  }
  Vector<double,1> max,min;
  max[0] = 1.0/0.5;
  min[0] = 1.0/1.5;
  re.generateConditionRegularInterval(max,min);
  for(int i=0;i<re.getLocalNumberOfReplicas();i++){
    re[i].setRandomVelocityAtTargetTemperature();
    re[i].calcForce();
    re[i].calcEnergy();
    re[i].calcPressure();
  }

  const int ninterval = 10000;
  for(int s=-1000;s<ninterval;s++){
    re.copyREMInfoToReplicas();
    if(s>=0) re.exchangeReplicasWithNeighbor();
    for(int i=0;i<re.getLocalNumberOfReplicas();i++){
      for(int ss=0;ss<500;ss++){
	re[i].scaleVelocity();
	re[i].integrateOneStep();
      }
      re[i].calcEnergy();
      re[i].calcPressure();
    }
    re.getREMInfoFromReplicas();
    if(s>=0)re.Output();
  }

  RE::Finalize();
}
