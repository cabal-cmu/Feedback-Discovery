/*  mcmc_mh.h

    Mark Woolrich FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(mcmc_mh_h)
#define mcmc_mh_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "nmaoptions.h"
#include "newimage/newimageall.h"
#include "newmat.h"
#include "mcmc_component.h"
#include "vb_component.h"

using namespace NEWIMAGE;
using namespace NEWMAT;

namespace Nma {

  template <class T> int write_num(const T& num, string p_fname)
  {
    ofstream out;
    out.open(p_fname.c_str(), ios::out);
    
    if(!out)
      {
	cerr << "Unable to open " << p_fname << endl;
	return -1;
      }

    out << num << endl;

    out.close();

    return 0;
  }

  class Mcmc
    {
    public:
      
      // constructor
      Mcmc(vector<Mcmc_Component*> pcomponents, Vb_Component& pvb_comp, Mcmc_Log_Likelihood& pmcmc_log_likelihood, int pnsamps, int pnburnin, int pnjumps_per_sample, int psample_vb_jump_every, int pburnin_vb_jump_every, int pdebuglevel, bool poutput_samples)	:
	nsamps(pnsamps),
	nburnin(pnburnin),
	njumps_per_sample(pnjumps_per_sample),
	sample_vb_jump_every(psample_vb_jump_every),
	burnin_vb_jump_every(pburnin_vb_jump_every),
	components(pcomponents),
	vb_component(pvb_comp),
	ncomponents(pcomponents.size()),
	mcmc_log_likelihood(pmcmc_log_likelihood),
	debuglevel(pdebuglevel),
	output_samples(poutput_samples)
	{
	}

      // load data from file in from file and set up starting values
      virtual void setup()=0;

      // runs the chain
      virtual void run()=0;

      virtual void save();

      // normalised evidence = unnormalised_model_evidence*exp(model_evidence_normalisation)
      void evaluate_model_evidence();
      double get_unnormalised_model_evidence() const {return unnormalised_model_evidence;}
      double get_model_evidence_normalisation() const {return model_evidence_normalisation;}
      double get_aic() const {return aic;}
      double get_bic() const {return bic;}

      void set_values_to_sample_mean();
      void set_values_to_sample_map();

      const vector<Mcmc_Component*>& get_components() const {return components;}

      void find_parameter(const string& name, int& component_number, int& parameter_number) const;

      // Destructor
      virtual ~Mcmc() {
	//	LOGOUT("~Mcmc start");
	//	LOGOUT("~Mcmc end");
      }

      const vector<double>& get_energy_hist() const {return energy_hist;}
      const Mcmc_Log_Likelihood& get_mcmc_log_likelihood() const {return mcmc_log_likelihood;}

      int get_nsamps() const {return nsamps;}
      int get_nburnin() const {return nburnin;}

    protected:    
 
      int nsamps;  // number of samps not including burnin
      int nburnin; // number of samps in burnin
      int njumps_per_sample;
      int sample_vb_jump_every;
      int burnin_vb_jump_every;

      vector<Mcmc_Component*> components;
      Vb_Component& vb_component;

      int ncomponents;      
  
      Mcmc_Log_Likelihood& mcmc_log_likelihood;

      int debuglevel;

      bool output_samples;

      vector<double> energy_hist; // nsamples

      // normalised evidence = unnormalised_model_evidence*exp(model_evidence_normalisation)
      double unnormalised_model_evidence;
      double model_evidence_normalisation;

      double aic;
      double bic;

    private:
      Mcmc();
      const Mcmc& operator=(Mcmc& mcmc_mh);     
      Mcmc(Mcmc& mcmc_mh);

    };

  class Mcmc_Mh : public Mcmc
    {
    public:

      // constructor
      Mcmc_Mh(vector<Mcmc_Component*> pcomponents, Vb_Component& pvb_comp, Mcmc_Log_Likelihood& pmcmc_log_likelihood, int pnsamps, int pnburnin, int pnjumps_per_sample, int psample_vb_jump_every, int pburnin_vb_jump_every, int pdebuglevel, bool poutput_samples)	:
	Mcmc(pcomponents,pvb_comp,pmcmc_log_likelihood,pnsamps,pnburnin,pnjumps_per_sample, psample_vb_jump_every,pburnin_vb_jump_every,pdebuglevel,poutput_samples)
	{
	  setup();
	}

      // load data from file in from file and set up starting values
      void setup();

      // runs the chain
      void run();

      // Destructor
      virtual ~Mcmc_Mh() {
	//	LOGOUT("~Mcmc_Mh start");
	//	LOGOUT("~Mcmc_Mh end");
      }

    private:    
 
      Mcmc_Mh();
      const Mcmc_Mh& operator=(Mcmc_Mh& mcmc_mh);     
      Mcmc_Mh(Mcmc_Mh& mcmc_mh);

    };
}   
#endif

