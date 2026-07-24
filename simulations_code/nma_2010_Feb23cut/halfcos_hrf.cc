/*  halfcos_hrf.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT  */

#include "miscmaths/miscmaths.h"
#include "utils/tracer_plus.h"
#include "utils/log.h"
#include "halfcos_hrf.h"
//#include "fftw/fftw.h"

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace Utilities;

namespace Nma {

  void Halfcos_Hrf::fft_stimulus(const ColumnVector& stimulus, ColumnVector& fftstimulus_real,  ColumnVector& fftstimulus_imag, double& stim_mean) const
  {
    Tracer_Plus trace("Halfcos_Hrf::fft_stimulus");
    
    if(ntpts!=stimulus.Nrows())
      {
	OUT(ntpts);
	OUT(stimulus.Nrows());
	cout << "Stimulus had wrong number of timepoints" << endl;
	throw Exception(string("Stimulus had wrong number of timepoints").data());
      }
    
    // Remove and store mean
    stim_mean=MISCMATHS::mean(stimulus).AsScalar();
    
    real.ReSize(zeropad);
    real = 0;
    real.Rows(1,ntpts) = stimulus - stim_mean;
  
    ColumnVector zeros(zeropad); zeros=0;
    
    FFT(real, zeros, fftstimulus_real, fftstimulus_imag);
    //RealFFT(real, fftstimulus_real, fftstimulus_imag);
    
//     write_ascii_matrix(stimulus,LogSingleton::getInstance().appendDir("stimulus"));
//     write_ascii_matrix(real,LogSingleton::getInstance().appendDir("real"));
   
  }   

  void Halfcos_Hrf::set_nsecs(double pnsecs) {
    Tracer_Plus trace("Halfcos_Hrf::set_nsecs");

	nsecs=pnsecs;
	if(std::abs((tr/res)-::round(tr/res))/::round(tr/res)>1e-3)
	  {
	    LOGOUT(abs((tr/res)-::round(tr/res)));
	    LOGOUT(tr/res);
	    LOGOUT(::round(tr/res));
	    LOGOUT(nsecs/res);
	    LOGOUT(tr);
	    LOGOUT(res);
	    LOGOUT("In set_nsecs(), tr/Resolution is not a whole number and would require interpolation - which is currently not supported");
	    throw Exception(string("tr/Resolution is not a whole number and would require interpolation - which is currently not supported").data());
	  }

	if(std::abs((nsecs/res)-::round(nsecs/res))/::round(nsecs/res)>1e-3)
	  {
	    LOGOUT(abs((nsecs/res)-::round(nsecs/res)));
	    LOGOUT(nsecs/res);
	    LOGOUT(::round(nsecs/res));
	    LOGOUT(nsecs);
	    LOGOUT(res);
	    LOGOUT("In set_nsecs(), nsecs/Resolution is not a whole number and would require interpolation - which is currently not supported");
	    throw Exception(string("nsecs/Resolution is not a whole number and would require interpolation - which is currently not supported").data());
	  }
	
	nscans = int(::round(nsecs/tr));
	ntpts = int(::round(nsecs/res));
	
	LOGOUT(nscans);
	LOGOUT(nsecs);

	// setup zeropadding
	zeropad = establishZeroPadding(ntpts);
	subsampledtpts = int(::round(zeropad/(tr/res)));

// 	in = new fftw_complex[zeropad];
// 	out = new fftw_complex[zeropad];
//         for(int c=0;c<zeropad;c++){
// 	  c_re(in[c]) = 0;
// 	  c_re(out[c]) = 0;
// 	  c_im(in[c]) = 0;
// 	  c_im(out[c]) = 0;
// 	}
// 	plan = fftw_create_plan(zeropad, FFTW_BACKWARD, FFTW_MEASURE);
	
  }

  void Halfcos_Hrf::halfcos_fft(ColumnVector& preal, ColumnVector& pimag, float ybot, float ytop, float xleft, float xright, int flipud) const
  {
    Tracer_Plus trace("Halfcos_Hrf::halfcos_fft");

    if(nsecs==0)
      throw Exception(string("nsecs not set").data());

    if(flipud==-1)
      {
	float tmp=ybot;
	ybot=ytop;
	ytop=tmp;
      }

    preal=0;
    pimag=0;

    if(xright!=xleft && ytop!=ybot)
      {
	float ntow = 2*M_PI/(zeropad*res);

	for(int n = 2; n <= zeropad/2; n++)
	  {
	    //    cerr << "w=" << w << ",ybot=" << ybot << ",ytop=" << ytop << "xleft=" << xleft << "xleft=" << xleft << "flipud=" << flipud << endl;
	    float w = ntow*(n-1);

	    float tmp1 = (2*w*(Sqr(M_PI)-Sqr(w*(xright-xleft))));
	    float tmp2 = Sqr(M_PI)*(ytop+ybot);
	    float tmp3 = 2*Sqr(w*(xright-xleft));
	    float a1 = (ybot*tmp3-tmp2)/tmp1;
	    float a2 = (ytop*tmp3-tmp2)/tmp1;

	    //    cerr << "a2=" << a2 << endl;
	    float b1 = -M_PI/2.0-xright*w;
	    //    cerr << "b1=" << b1 << endl;
	    float b2 = -M_PI/2.0-xleft*w;
	    //    cerr << "b2=" << b2 << endl;
	    
	    preal(n) = a1*cos(b1)-a2*cos(b2);
	    pimag(n) = a1*sin(b1)-a2*sin(b2);
	  }
      }
  }

  
  void Halfcos_Hrf::fft_hrf(float m1, float m2, float m3, float m4, float c1, float c2, ColumnVector& ffthrf_real, ColumnVector& ffthrf_imag) const
  {
    Tracer_Plus trace("Halfcos_Hrf::fft_hrf");
    
   if(nsecs==0)
      throw Exception(string("nsecs not set").data());

    ColumnVector hrfw1_real(zeropad/2);
    ColumnVector hrfw1_imag(zeropad/2);
    halfcos_fft(hrfw1_real,hrfw1_imag,-1*c1,0,0,m1,1);
    ColumnVector hrfw2_real(zeropad/2);
    ColumnVector hrfw2_imag(zeropad/2);
    halfcos_fft(hrfw2_real,hrfw2_imag,-1*c1,1,m1,m1+m2,-1);
    ColumnVector hrfw3_real(zeropad/2);
    ColumnVector hrfw3_imag(zeropad/2);
    halfcos_fft(hrfw3_real,hrfw3_imag,-1*c2,1,m1+m2,m1+m2+m3,1);
    ColumnVector hrfw4_real(zeropad/2);
    ColumnVector hrfw4_imag(zeropad/2);
    halfcos_fft(hrfw4_real,hrfw4_imag,-1*c2,0,m1+m2+m3,m1+m2+m3+m4,-1);

    ffthrf_real.ReSize(zeropad);
    ffthrf_real=0;
    ffthrf_imag.ReSize(zeropad);
    ffthrf_imag=0;

    for(int n = 2; n <= zeropad/2; n++)
      {	       
	ffthrf_real(n) = hrfw1_real(n)+hrfw2_real(n)+hrfw3_real(n)+hrfw4_real(n);
	ffthrf_imag(n) = hrfw1_imag(n)+hrfw2_imag(n)+hrfw3_imag(n)+hrfw4_imag(n);       
      }    
  }

  void Halfcos_Hrf::convolve_hrf_fft(const ColumnVector& fftstimulus_real,  const ColumnVector& fftstimulus_imag, double stim_mean, const ColumnVector& ffthrf_real, const ColumnVector& ffthrf_imag, ColumnVector& response) const
  {    
    Tracer_Plus trace("Halfcos_Hrf::convolve_hrf_fft");        
    
   if(nsecs==0)
      throw Exception(string("nsecs not set").data());
   
   // Convolve
   FFTI(SP(fftstimulus_real, ffthrf_real)-SP(fftstimulus_imag, ffthrf_imag), 	 SP(fftstimulus_imag, ffthrf_real)+SP(fftstimulus_real, ffthrf_imag),	 real, imag);
   //RealFFTI(SP(fftstimulus_real, ffthrf_real)-SP(fftstimulus_imag, ffthrf_imag), 	 SP(fftstimulus_imag, ffthrf_real)+SP(fftstimulus_real, ffthrf_imag),	 real);
   
//    OUT(real.Nrows());
//    OUT(nscans);

//    OUT(zeropad);
//    OUT(subsampledtpts);

   // take first ntpts and subsample
   float ss = float(zeropad)/float(subsampledtpts);

   response.ReSize(nscans);

//    write_ascii_matrix(real,LogSingleton::getInstance().appendDir("response_hr"));

   // sub sample    
   for(int scan = 1; scan <= nscans; scan++)
     {
	response(scan) = real(int(::round(scan*ss)));
     }
   

//    write_ascii_matrix(response,LogSingleton::getInstance().appendDir("response"));
     
   response += stim_mean;
   
   //  write_ascii_matrix(response,LogSingleton::getInstance().appendDir("response"));
   //     write_ascii_matrix(getSubsampledStimulus(0),LogSingleton::getInstance().appendDir("stimulus"));    
   
  }

// void Halfcos_Hrf::convolve_hrf_fftw(const ColumnVector& fftstimulus_real,  const ColumnVector& fftstimulus_imag, double stim_mean, const ColumnVector& ffthrf_real, const ColumnVector& ffthrf_imag, ColumnVector& response) const
//   {    
//     Tracer_Plus trace("Halfcos_Hrf::convolve_hrf_fftw");        
    
//    if(nsecs==0)
//       throw Exception(string("nsecs not set").data());
   
//     fftw_one(plan, in, out);

//    // Convolve
// //    FFTI(SP(fftstimulus_real, ffthrf_real)-SP(fftstimulus_imag, ffthrf_imag), 	 SP(fftstimulus_imag, ffthrf_real)+SP(fftstimulus_real, ffthrf_imag),	 real, imag);
//    //RealFFTI(SP(fftstimulus_real, ffthrf_real)-SP(fftstimulus_imag, ffthrf_imag), 	 SP(fftstimulus_imag, ffthrf_real)+SP(fftstimulus_real, ffthrf_imag),	 real);
   
//    // take first ntpts and subsample
//    float ss = zeropad/subsampledtpts;
   
//    response.ReSize(nscans);
   
//    for(int scan = 1; scan <= nscans; scan++)
//      {
//        response(scan) = real(int(scan*ss));
//      }
   
//    response += stim_mean;
   
//    //  write_ascii_matrix(response,LogSingleton::getInstance().appendDir("response"));
//    //     write_ascii_matrix(getSubsampledStimulus(0),LogSingleton::getInstance().appendDir("stimulus"));    
   
//   }

  void Halfcos_Hrf::hrf_ifft(const ColumnVector& ffthrf_real, const ColumnVector& ffthrf_imag, ColumnVector& response) const
  {    
    Tracer_Plus trace("Halfcos_Hrf::hrf_ifft");        
    
   if(nsecs==0)
      throw Exception(string("nsecs not set").data());

   ColumnVector real_highres(ffthrf_real.Nrows());
   real_highres=0;
   ColumnVector imag_highres(ffthrf_imag.Nrows());
   imag_highres=0;
   
   FFTI(ffthrf_real, ffthrf_imag, real_highres, imag_highres);
   //RealFFTI(ffthrf_real, ffthrf_imag, real_highres);
   
   // take first ntpts
   response.ReSize(ntpts);
   
   for(int t = 1; t <= ntpts; t++)
     {
       response(t) = real_highres(t);
     }
   
   //  write_ascii_matrix(response,LogSingleton::getInstance().appendDir("response"));
   
  }
  
}
