#include <iostream>

#include "janus.hpp"

#include <vector.hpp>
#include <replica_exchanger.hpp>

int main(int argc,char **argv){
  RE::Initialize(argc,argv);

  int ninterval = 10000;
  int nstep = 100;
  int ninterval_eq = 1000;
  int chk_interval = ninterval;
  int cdv_interval = ninterval;
  std::string patchy_info_file = "3patch.inp";
  std::string input_prefix = "./",output_prefix = "./";
  bool flag_read = false;
  int N = 100;
  Vector<int,2> ndim(1,1);
  Vector<double,2> max(1.0,1.0),min(1.0,1.0);

  for(int count=1;count < argc;count++){
    std::string opt = argv[count];
    if(opt == "--ninterval" || opt == "-n"){
      ninterval = atoi(argv[++count]);
      fprintf(stderr,"ninterval:\t%d\n",ninterval);
      continue;
    }
    if(opt == "--ninterval_eq" || opt == "-e"){
      ninterval_eq = atoi(argv[++count]);
      fprintf(stderr,"ninterval_eq:\t%d\n",ninterval_eq);
      continue;
    }
    if(opt == "--nstep" || opt == "-s"){
      nstep = atoi(argv[++count]);
      fprintf(stderr,"nstep:\t\t%d\n",nstep);
      continue;
    }
    if(opt == "--patchy_info_file" || opt == "-f"){
      patchy_info_file = argv[++count];
      fprintf(stderr,"patchy_info_file:\t%s\n",patchy_info_file.c_str());
      continue;
    }
    if(opt == "-N"){
      N = atoi(argv[++count]);
      fprintf(stderr,"N:\t\t%d\n",N);
      continue;
    }
    if(opt == "--dimension" || opt == "-d"){
      ndim[0] = atoi(argv[++count]);
      assert(count < argc);
      ndim[1] = atoi(argv[++count]);
      fprintf(stderr,"ndim:\t\t%d %d\n",ndim[0],ndim[1]);
      continue;
    }
    if(opt == "--max"){
      max[0] = atof(argv[++count]);
      assert(count < argc);
      max[1] = atof(argv[++count]);
      fprintf(stderr,"max:\t\t%lf %lf\n",max[0],max[1]);
      continue;
    }
    if(opt == "--min"){
      min[0] = atof(argv[++count]);
      assert(count < argc);
      min[1] = atof(argv[++count]);
      fprintf(stderr,"min:\t\t%lf %lf\n",min[0],min[1]);
      continue;
    }
    if(opt == "--input" || opt == "-i"){
      input_prefix = argv[++count];
      flag_read = true;
      fprintf(stderr,"input_prefix:\t%s\n",input_prefix.c_str());
      fprintf(stderr,"read mode is on\n");
      continue;
    }
    if(opt == "--output" || opt == "-o"){
      output_prefix = argv[++count];
      fprintf(stderr,"output_prefix:\t%s\n",output_prefix.c_str());
      fprintf(stderr,"read mode is on\n");
      continue;
    }
    if(opt == "--snapshot-interval" || opt == "-S"){
      cdv_interval = atoi(argv[++count]);
      fprintf(stderr,"snapshot_interval:\t%d\n",cdv_interval);
      continue;
    }
    if(opt == "--checkpoint-interval" || opt == "-c"){
      chk_interval = atoi(argv[++count]);
      fprintf(stderr,"checkpoint_interval:\t%d\n",chk_interval);
      continue;
    }
    if(opt == "--help" || opt == "-h"){
      fprintf(stderr,"available options:\n");
      fprintf(stderr,"--input               | -i: prefix of input chk files\n");
      fprintf(stderr,"--output              | -o: prefix of output files\n");
      fprintf(stderr,"--ninterval           | -n: number of intervals for replica exchange\n");
      fprintf(stderr,"--ninterval_eq        | -e: number of intervals for equilibration\n");
      fprintf(stderr,"--nstep               | -s: number of steps in a interval\n");
      fprintf(stderr,"--patchy_info_file    | -p: file name of patchy info file\n");
      fprintf(stderr,"-N:                         number of particles\n");
      fprintf(stderr,"--max:                      max Beta and PBeta(e.g. --max 1.0 2.0)\n");
      fprintf(stderr,"--min:                      min Beta and PBeta(e.g. --min 0.5 1.0)\n");
      fprintf(stderr,"--snapshot_interval   | -S: snapshot is generated at every this interval\n");
      fprintf(stderr,"--checkpoint_interval | -c: checkpoint is generated at every this interval\n");
      exit(0);
    }
    fprintf(stderr,"error: undefined option %s\n",opt.c_str());
    exit(EXIT_FAILURE);
  }
  if(max[0] < min[0]) max[0] = min[0];
  if(max[1] < min[1]) max[1] = min[1];

  PatchyInfo pinfo(patchy_info_file);
  RE::ReplicaExchanger<2,Janus<Condition,Result>,Condition,Result> re;
  re.Initialize(ndim);
  for(int i=0;i<re.getLocalNumberOfReplicas();i++){
    re[i].initialize(N,pinfo);
  }
  re.generateConditionRegularInterval(max,min);
  re.copyREMInfoToReplicas();
  for(int i=0;i<re.getLocalNumberOfReplicas();i++){
    if(flag_read){
      re[i].readCheckPoint(input_prefix+re[i].getCondition().getPrefix()+".chk");
    }else{
      re[i].setParticlePlane();
      re[i].setRandomVelocityAtTargetTemperature();
      re[i].calcForce();
      re[i].calcEnergy();
      re[i].calcPressure();
      re[i].adjustVolume();
    }
  }

  int next_chk = 0;
  int next_cdv = 0, cdv_step = 0;
  for(int s=-ninterval_eq;s<ninterval;s++){
    re.copyREMInfoToReplicas();
    if(s>=0){
      re.exchangeReplicasWithNeighbor();
    }
    for(int i=0;i<re.getLocalNumberOfReplicas();i++){
      if(s<0){
	re[i].adjustVolume();
	re[i].calcKineticEnergy();
	re[i].scaleVelocity();
      }else{
	re[i].scaleVelocityAfterExchange();
      }
      for(int ss=0;ss<nstep;ss++){
	re[i].integrateOneStep();
      }
      re[i].calcEnergy();
      re[i].calcPressure();
    }
    re.getREMInfoFromReplicas();
    if(s>=0) re.Output(output_prefix);

    if(next_chk == s){
      for(int i=0;i<re.getLocalNumberOfReplicas();i++){
	re[i].outputCheckPoint(output_prefix+re[i].getCondition().getPrefix()+".chk");
      }
      next_chk += chk_interval;
    }
    if(next_cdv == s){
      for(int i=0;i<re.getLocalNumberOfReplicas();i++){
	std::stringstream strs;
	strs << output_prefix << re[i].getCondition().getPrefix();
	strs << "_" << std::setfill('0') << std::setw(8) << cdv_step++;
	strs << ".cdv";
	re[i].outputCDV(strs.str());
      }
      next_cdv += cdv_interval;
    }
  }

  for(int i=0;i<re.getLocalNumberOfReplicas();i++){
    re[i].outputCheckPoint(output_prefix+re[i].getCondition().getPrefix()+".chk");
  }
  RE::Finalize();
}
