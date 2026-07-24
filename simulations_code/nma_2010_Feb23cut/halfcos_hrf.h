/*  halfcos_hrf.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(halfcos_hrf_h)
#define halfcos_hrf_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "miscmaths/miscmaths.h"
#include "utils/log.h"
//#include "fftw/fftw.h"

using namespace NEWMAT;
using namespace MISCMATHS;
using namespace Utilities;

namespace Nma {
    
  class Halfcos_Hrf
    {
    public:

      // constructor
      Halfcos_Hrf(){}

      Halfcos_Hrf(double pres, float ptr)
	:res(pres),
	 tr(ptr),
	 nsecs(0)
	{
	}

      Halfcos_Hrf(const Halfcos_Hrf& hrf)
	{
	  *this=hrf;
	}

      const Halfcos_Hrf& operator=(const Halfcos_Hrf& hrf)
      {
	res=hrf.res;
	tr=hrf.tr;	
	nsecs=hrf.nsecs;
	nscans=hrf.nscans;
	ntpts=hrf.ntpts;
	zeropad=hrf.zeropad;
	subsampledtpts=hrf.subsampledtpts;

	return *this;
      }

      virtual ~Halfcos_Hrf()
      {
// 	delete [] in;
// 	delete [] out;
      }

      void set_nsecs(double pnsecs); 

      void fft_hrf(float m1, float m2, float m3, float m4, float c1, float c2, ColumnVector& ffthrf_real, ColumnVector& ffthrf_imag) const;

      void fft_stimulus(const ColumnVector& stimulus, ColumnVector& fftstimulus_real,  ColumnVector& fftstimulus_imag, double& stim_mean) const;

      void convolve_hrf_fft(const ColumnVector& fftstimulus_real,  const ColumnVector& fftstimulus_imag, double stim_mean, const ColumnVector& ffthrf_real, const ColumnVector& ffthrf_imag, ColumnVector& response) const;
      //      void convolve_hrf_fftw(const ColumnVector& fftstimulus_real,  const ColumnVector& fftstimulus_imag, double stim_mean, const ColumnVector& ffthrf_real, const ColumnVector& ffthrf_imag, ColumnVector& response) const;

      void hrf_ifft(const ColumnVector& ffthrf_real, const ColumnVector& ffthrf_imag, ColumnVector& response) const;
          
      int get_nscans() const {return nscans;}
      float get_tr() const {return tr;}
 
    private:    
      
      void halfcos_fft(ColumnVector& preal, ColumnVector& pimag, float ybot, float ytop, float xleft, float xright, int flipud) const;

      int establishZeroPadding(int in) const
      {
	return (int)pow(2,ceil(log(in)/log(2)));
      }

      double res; // in secs
      float tr; // in secs
      double nsecs; // in secs

      int nscans; // num of time points at fmri low res
      int ntpts; // num of time points at high res
      int subsampledtpts; // num of timepoints in zeropadded fmri low res
      
      int zeropad;
   
      mutable ColumnVector real;
      mutable ColumnVector imag;

//       fftw_complex *in, *out;
//       fftw_plan plan;
        
    };

}   

#endif







