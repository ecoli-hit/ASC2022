#include "map_aparam.h"
#ifndef _OPENMP
#include <omp.h> 
#endif
template <typename FPTYPE>
void deepmd::map_aparam_cpu (
    FPTYPE * output,
    const FPTYPE * aparam,
    const int * nlist,
    const int & nloc,
    const int & nnei,
    const int & numb_aparam
    )
//
//	output:	nloc x nnei x numb_aparam
//	aparam:	nall x numb_aparam
//	nlist:	nloc x nnei
//
{
  int tmp = nloc * nnei * numb_aparam;
  #pragma omp paraller for simd 
  for (int ii = 0; ii < tmp; ++ii){
    output[ii] = 0.;
  }

  // loop over loc atoms
  for (int ii = 0; ii < nloc; ++ii){
    int i_idx = ii;	
    // loop over neighbor atoms
    for (int jj = 0; jj < nnei; ++jj){
      int j_idx = nlist[i_idx * nnei + jj];
      if (j_idx < 0) continue;
      // loop over elements of aparam
      int tmp1 = ii * nnei * numb_aparam + jj * numb_aparam;
      int tmp2 = j_idx * numb_aparam;
      for (int dd = 0; dd < numb_aparam; ++dd){
	      output[tmp1 + dd] = aparam[tmp2 + dd];
      }
    }
  }  
}

template
void deepmd::map_aparam_cpu<double> (
    double * output,
    const double * aparam,
    const int * nlist,
    const int & nloc,
    const int & nnei,
    const int & numb_aparam
    );

template
void deepmd::map_aparam_cpu<float> (
    float * output,
    const float * aparam,
    const int * nlist,
    const int & nloc,
    const int & nnei,
    const int & numb_aparam
    );

