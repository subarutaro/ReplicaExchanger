#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#ifdef FMRES_MPI_PARALLEL
#include <mpi.h>
#endif

#include "vector.hpp"

#include <cmath>
#include <cassert>
#include <algorithm>

double random_number(){
  return (double)rand()/(double)RAND_MAX;
}

namespace RE{
  void Initialize(int argc,char **argv){
    #ifdef FMRES_DEBUG_PRINT
    fprintf(stderr,"note: initializing ReplicaExchanger\n");
    #endif
    #ifdef FMRES_MPI_PARALLEL
    MPI::Init(argc,argv);
    #endif
  }
  void Finalize(){
    #ifdef FMRES_DEBUG_PRINT
    fprintf(stderr,"note: finializing ReplicaExchanger\n");
    #endif
#ifdef FMRES_MPI_PARALLEL
    MPI::Finalize();
#endif
  }

  class AcceptanceRatio{
  private:
    int naccept = 0;
    int ndeny   = 0;
  public:
    void Accept(){naccept++;}
    void Deny(){ndeny++;}

    double get(){(double)naccept / (double)(naccept + ndeny);}
  };

  template<int NDIM,class TReplica,class TCondition,class TResult>
  class ReplicaExchanger{
  private:
    typedef Vector<double,NDIM> dvec;
    typedef Vector<   int,NDIM> ivec;
    class REMInfo{
    public:
      typedef TCondition Condition;
      typedef TResult    Result;
      int cindex,rindex;
      Condition condition;
      Result    result;

      REMInfo(){}
      REMInfo(const REMInfo& r){
	cindex    = r.cindex;
	rindex    = r.rindex;
	condition = r.condition;
	result    = r.result;
      }
    };
    REMInfo *reminfo  = nullptr;
    TReplica *replica = nullptr;
    AcceptanceRatio *accept_ratio[NDIM];

    std::fstream *output = nullptr;

    int *nreplica, *offset;
    int nreplica_total;
    int nreplica_local, offset_local;
    int nprocs, myrank;
    ivec ndim;

  public:
    TReplica& operator[](const int i){ return replica[i]; }
    const TReplica& operator[](const int i) const { return replica[i]; }

    ReplicaExchanger(){
      for(int i=0;i<NDIM;i++) accept_ratio[i] = nullptr;
    };
    ~ReplicaExchanger(){
      if(replica != nullptr) delete[] replica;
      if(reminfo != nullptr) delete[] reminfo;
      for(int i=0;i<NDIM;i++){
	if(accept_ratio[i] != nullptr)
	  delete[] accept_ratio[i];
      }

      if(output != nullptr){
	for(int i=0;i<nreplica_local;i++){
	  if(output[i].is_open()) output[i].close();
	}
	delete[] output;
      }
    }

    void Initialize(const ivec _ndim){
      ndim = _ndim;
      nreplica_total = prod(ndim);

      reminfo = new REMInfo[nreplica_total];
      for(int i=0;i<nreplica_total;i++){
	reminfo[i].rindex = reminfo[i].cindex = i;
      }

      for(int i=0;i<NDIM;i++){
	accept_ratio[i] = new AcceptanceRatio[nreplica_total];
      }

#ifdef FMRES_MPI_PARALLEL
      nprocs = MPI::COMM_WORLD.Get_size();
      myrank = MPI::COMM_WORLD.Get_rank();
      nreplica = new int[nprocs];
      offset   = new int[nprocs+1];
      offset[0] = 0;
      for(int i=0;i<nprocs;i++){
	nreplica[i] = nreplica_total / nprocs;
	if(i < nreplica_total%nprocs){
	  nreplica[i]++;
	}
	offset[i+1] = offset[i] + nreplica[i];
      }
      nreplica_local = nreplica[myrank];
      offset_local = offset[myrank];
      assert(nreplica_total == offset[nprocs]);
#else
      nprocs = 1;
      myrank = 0;
      nreplica = new int[1];
      offset = new int[1];
      nreplica[0] = nreplica_local = nreplica_total;
      offset[0] = offset_local = 0;
#endif
      replica = new TReplica[nreplica_local];

#ifdef FMRES_DEBUG_PRINT
      fprintf(stderr,"nprocs:\t%d\n",nprocs);
      fprintf(stderr,"myrank:\t%d\n",myrank);
      fprintf(stderr,"offset:\t%d\n",offset_local);
      fprintf(stderr,"nreplica_total:\t%d\n",nreplica_total);
      fprintf(stderr,"nreplica_local:\t%d\n",nreplica_local);
#endif
    }

    void getREMInfoFromReplicas(){
#ifdef FMRES_MPI_PARALLEL
      REMInfo send_buffer[nreplica_local];
      for(int i=0;i<nreplica_local;i++){
	const int index = i + offset_local;
	//reminfo[index].condition = replica[i].getCondition();
	reminfo[index].result = replica[i].getResult();

        send_buffer[i] = reminfo[index];
      }
      int send_count[nprocs];
      int send_offset[nprocs];
      for(int i=0;i<nprocs;i++){
        send_count[i]  = nreplica[i] * sizeof(REMInfo) / sizeof(float);
	send_offset[i] = offset[i]   * sizeof(REMInfo) / sizeof(float);
      }
      MPI::COMM_WORLD.Allgatherv(send_buffer, sizeof(REMInfo)/sizeof(float)*nreplica_local, MPI_FLOAT,
	reminfo,send_count,send_offset, MPI_FLOAT);
#else
      for(int i=0;i<nreplica_total;i++){
        reminfo[i].result = replica[i].getResult();
      }
#endif
    }

    void copyREMInfoToReplicas(){
      //sort replica information using replica id
      for(int i=0;i<nreplica_local;i++){
	const int index = offset_local + i;
        replica[i].setCondition(reminfo[index].condition);
      }
    }

    void exchangeReplicasWithNeighbor(){
      //sort replica information using condition id
      std::sort(reminfo,reminfo+nreplica_total,
		[](const REMInfo& a,const  REMInfo& b){return a.cindex < b.cindex;});

      static int dim = 0;
      static int oddoreven = 0;
      //exchange conditions between neighbor
      if(myrank == 0){
	for(int i=0;i<nreplica_total;i++){
	  ivec indices = index2Indices(i,ndim);
	  if(indices[dim]%2 != oddoreven || indices[dim] == ndim[dim]-1) continue;
	  const int a = i;
	  indices[dim]++;
	  if(indices[dim]>=ndim[dim]) indices[dim] -= ndim[dim];
	  const int b = indices2Index(indices,ndim);
	  #ifdef FMRES_DEBUG_PRINT
	  fprintf(stderr,"trying to exchange condition %d and %d\n",a,b);
	  #endif
	  REMInfo& ri_a = reminfo[a];
	  REMInfo& ri_b = reminfo[b];
	  const double bf_aa = ri_a.result.getLogBoltzmannFactor(ri_a.condition);
	  const double bf_bb = ri_b.result.getLogBoltzmannFactor(ri_b.condition);
	  const double bf_ab = ri_a.result.getLogBoltzmannFactor(ri_b.condition);
	  const double bf_ba = ri_b.result.getLogBoltzmannFactor(ri_a.condition);
	  //const double prob = bf_ab*bf_ba / (bf_aa*bf_bb);
	  const double delta = bf_aa - bf_ba + bf_bb - bf_ab;
	  const double prob = exp(-delta);
          #ifdef FMRES_DEBUG_PRINT
	  //fprintf(stderr,"exchange prob = %lf : %lf %lf %lf %lf\n",prob,bf_aa,bf_bb,bf_ab,bf_ba);
	  #endif
	  if(prob >= random_number()){
	    #ifdef FMRES_DEBUG_PRINT
	    fprintf(stderr,"replica %d and %d is exchanged\n",ri_a.rindex,ri_b.rindex);
	    #endif
	    std::swap(ri_a.rindex,ri_b.rindex);
	    accept_ratio[dim][a].Accept();
	  }else{
	    accept_ratio[dim][a].Deny();
	  }
	}
	dim++;
	if(dim==NDIM){
	  oddoreven = (oddoreven==0) ? 1 : 0;
	  dim = 0;
	}
      }
      #ifdef FMRES_DEBUG_PRINT
      static int s=0;
      if(myrank == 0){
	printf("rindex %d",s);
	for(int i=0;i<nreplica_total;i++){
	  printf(" %d",reminfo[i].rindex);
	}
	printf("\n");
      }
      #endif
      #ifdef FMRES_MPI_PARALLEL
      MPI::COMM_WORLD.Bcast(reminfo,nreplica_total*sizeof(REMInfo)/sizeof(float),MPI_FLOAT,0);
      /*
      #ifdef FMRES_DEBUG_PRINT
      MPI::COMM_WORLD.Barrier();
      if(myrank == 1){
        for(int j=0;j<nreplica_total;j++){
	  printf("rank %d: replica %d T = %lf, ri = %d, ci = %d after exchange\n",myrank,j,reminfo[j].condition.T,reminfo[j].rindex,reminfo[j].cindex);
        }
      }
      MPI::COMM_WORLD.Barrier();
      #endif // DEBUG
      //*/
      #endif

      std::sort(reminfo,reminfo+nreplica_total,
		[](const REMInfo& a,const REMInfo& b){return a.rindex < b.rindex;});

      #ifdef FMRES_DEBUG_PRINT
      if(myrank == 0){
	printf("cindex %d",s);
	for(int i=0;i<nreplica_total;i++){
	  printf(" %d",reminfo[i].cindex);
	}
	printf("\n");
      }
      s++;
      #endif
    }

    void generateConditionRegularInterval(dvec max,dvec min){
      #ifdef FMRES_DEBUG_PRINT
      fprintf(stderr,"note: generating condition at regular interval\n");
      #endif
      std::sort(reminfo,reminfo+nreplica_total,
		[](const REMInfo& a,const REMInfo& b){return a.cindex < b.cindex;});
      dvec diff;
      for(int i=0;i<NDIM;i++){
	if(ndim[i]>1) diff[i] = (max[i] - min[i]) / (double)(ndim[i] - 1);
	else diff[i] = 0.0;
      }
      for(int i=0;i<nreplica_total;i++){
	ivec indices = index2Indices(i,ndim);
	dvec c = min + diff * indices;
	reminfo[i].condition.set(c);
	#ifdef FMRES_DEBUG_PRINT
	std::cerr << reminfo[i].condition << std::endl;
	#endif
      }
      std::sort(reminfo,reminfo+nreplica_total,
		[](const REMInfo& a,const REMInfo& b){return a.rindex < b.rindex;});
    }

    void Output(std::string prefix = ""){
      #ifdef FMRES_DEBUG_PRINT
      for(int i=0;i<nreplica_total;i++){
      printf("rank %d: reminfo %d ri = %d, ci = %d, T = %lf\n",myrank,i,reminfo[i].rindex,reminfo[i].cindex,reminfo[i].condition.T);
      }
      #ifdef FMRES_MPI_PARALLEL
      MPI::COMM_WORLD.Barrier();
      #endif
      #endif

      std::sort(reminfo,reminfo+nreplica_total,
		[](const REMInfo& a,const REMInfo& b){return a.cindex < b.cindex;});

      #ifdef FMRES_DEBUG_PRINT
      #ifdef FMRES_MPI_PARALLEL
      MPI::COMM_WORLD.Barrier();
      #endif
      for(int i=0;i<nreplica_total;i++){
      printf("rank %d: reminfo %d ri = %d, ci = %d, T = %lf\n",myrank,i,reminfo[i].rindex,reminfo[i].cindex,reminfo[i].condition.T);
      }
      #endif
#if 0 // output binary
      if(output == nullptr){
	output = new std::fstream[nreplica_local];
	for(int i=0;i<nreplica_local;i++){
	  const int index = i + offset_local;
	  assert(index == reminfo[index].cindex);
	  output[i].open(prefix+reminfo[index].condition.getPrefix()+".dat",std::ios::out | std::ios::binary);
	  assert(!output[i].fail());
	  output[i].write((char*)&reminfo[index].condition,sizeof(TCondition));
	}
      }
      for(int i=0;i<nreplica_local;i++){
	const int index = i + offset_local;
	assert(index == reminfo[index].cindex);
	output[i].write((char*)&reminfo[index].result,sizeof(TResult));
      }
#else // output ascii
      if(output == nullptr){
	output = new std::fstream[nreplica_local];
	for(int i=0;i<nreplica_local;i++){
	  const int index = i + offset_local;
	  assert(index == reminfo[index].cindex);
	  output[i].open(prefix+reminfo[index].condition.getPrefix()+".dat",std::ios::out);
	  assert(!output[i].fail());
	  output[i] << std::fixed;
	  output[i] << reminfo[index].condition << std::endl;
	}
      }
      for(int i=0;i<nreplica_local;i++){
	const int index = i + offset_local;
	assert(index == reminfo[index].cindex);
	output[i] << reminfo[index].condition << " " << reminfo[index].result << std::endl;
      }
#endif
      std::sort(reminfo,reminfo+nreplica_total,
		[](const REMInfo& a,const REMInfo& b){return a.rindex < b.rindex;});
    }

    void OutputAcceptanceRatio(std::string prefix = "./"){
      if(myrank == 0){
	std::sort(reminfo,reminfo+nreplica_total,
		  [](const REMInfo& a,const REMInfo& b){return a.cindex < b.cindex;});

	std::string filename = prefix + "acceptance_ratio.dat";
	std::ofstream ofs(filename);
	for(int d=0;d<NDIM;d++){
	  for(int i=0;i<nreplica_total;i++){
	    ivec index = index2Indices(i,ndim);
	    index[d]++;
	    if(index[d]>=ndim[d]) index[d] -= ndim[d];
	    const int j = indices2Index(index,ndim);
	    ofs << reminfo[i].condition << " <--> " << reminfo[j].condition << " " << accept_ratio[d][i].get() << std::endl;
	  }
	}

	std::sort(reminfo,reminfo+nreplica_total,
		  [](const REMInfo& a,const REMInfo& b){return a.rindex < b.rindex;});
      }
    }

    const int getTotalNumberOfReplicas() const { return nreplica_total; }
    const int getLocalNumberOfReplicas() const { return nreplica_local; }
    int getRank() const {
#ifdef FMRES_MPI_PARALLEL
      return MPI::COMM_WORLD.Get_rank();
#else
      return 0;
#endif
    }
    int getNumberOfProcs() const {
#ifdef FMRES_MPI_PARALLEL
      return MPI::COMM_WORLD.Get_size();
#else
      return 1;
#endif
    }
  };
};
