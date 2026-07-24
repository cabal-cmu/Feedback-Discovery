/*  mcmc_hmc.h

    Mark Woolrich FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(mcmc_hmc_h)
#define mcmc_hmc_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "nmaoptions.h"
#include "newimage/newimageall.h"
#include "newmat.h"
#include "mcmc_component.h"
#include "mcmc_mh.h"

using namespace NEWIMAGE;
using namespace NEWMAT;

namespace Nma {

  class Mcmc_hmc : public Mcmc
    {
    public:

      // constructor
      Mcmc_hmc(vector<Mcmc_Component*> pcomponents, Mcmc_Log_Likelihood& pmcmc_log_likelihood, int pnsamps, int pnburnin, int pnjumps_per_sample, int pdebuglevel, bool poutput_samples)	:
      Mcmc(pcomponents,pmcmc_log_likelihood,pnsamps,pnburnin,pnjumps_per_sample,-1,-1,pdebuglevel,poutput_samples)
	{
	  setup();
	}

      // load data from file in from file and set up starting values
      void setup();

      // runs the chain
      void run();

      void jump(double& energy);

      void initialise_epsilons();

      double calc_energy();
      
      bool in_bounds();

      void calc_gradient(ColumnVector& gradient, double energy);

      void set_values(const RowVector& values);

      // Destructor
      virtual ~Mcmc_hmc() {
	//	LOGOUT("~Mcmc_hmc start");
	//	LOGOUT("~Mcmc_hmc end");
      }

    private:    
 
      Mcmc_hmc();
      const Mcmc_hmc& operator=(Mcmc_hmc& mcmc_mh);     
      Mcmc_hmc(Mcmc_hmc& mcmc_mh);

      int nparams;

      RowVector values;

      ColumnVector epsilons;

      vector<string> parameter_names;

      int naccepted;
      int nrejected;
    };  

}   
#endif

