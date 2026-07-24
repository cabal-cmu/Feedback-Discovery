/*  mcmc_component.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */
  
#if !defined(mcmc_component_h)
#define mcmc_component_h

#include "utils/tracer_plus.h"
#include "miscmaths/miscprob.h"

using namespace Utilities;
using namespace MISCMATHS;

namespace Nma {

  class Mcmc_Component
    {
    public:
      // Constructor
      Mcmc_Component(const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pcov_mode, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision);

      Mcmc_Component(const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pcov_mode, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision, const ColumnVector& pprior_isard);

      virtual ~Mcmc_Component(){}

      // setup
      virtual void setup();

      // move to next point in MCMC
      virtual double jump(double old_energy);
      virtual double vb_jump(double old_energy);

      // sample at the current point in the chain
      virtual void sample();

      // save data to logger dir
      virtual void save();
      
      const vector<float>& get_gradient_step_sizes() const {return gradient_step_sizes;}

      void update_proposal_stds();
      void update_proposal_cov(int njumps_per_update);

      virtual double calc_prior_energy(const ColumnVector& pvalues, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision){
	ColumnVector tmp=pvalues-pprior_mean;
	return 0.5*(tmp.t()*pprior_precision*tmp).AsScalar();
 
      }

      virtual double calc_energy() = 0;

      void set_values_to_sample_mean();

      void set_values_to_sample(int index);
      
      virtual bool in_bounds() const {
	bool ret=true;
	for(int n=0; n<nparams && ret; n++)
	  ret=(values[n]<param_max[n] && values[n]>param_min[n]);	

	return ret;
      }

      const vector<string>& get_parameter_names() const {return parameter_names;}
      const string& get_name() const {return name;}
      const vector<vector<float> >& get_samples() const {return samples;} // num_params*num_samps
      float get_proposal_std() const {return proposal_std;}
      const SymmetricMatrix& get_proposal_cov() const {return proposal_cov;}
      const vector<float>& get_values() const {return values;}

      const vector<float>& get_param_min() const {return param_min;}
      const vector<float>& get_param_max() const {return param_max;}

      void set_debuglevel(int pdebuglevel) {debuglevel=pdebuglevel;}

      virtual void set_values_vec(const vector<float>& values_vec) {values=values_vec;}
     
      // store calculations for old values when new values are proposed
      virtual void store_old() = 0; 

      // restore calculations for old values when proposal is rejected
      virtual void restore_old() = 0;

      virtual void append_gradient(vector<float>& grad_vec, double energy);

      virtual void calc_jacobian(Matrix& jacob);
      virtual void calc_forward_model(ColumnVector& fm) = 0;
      virtual void calc_noise_precision(ColumnVector& prec)  = 0;
      virtual const ColumnVector& get_vb_data() = 0;

      const SymmetricMatrix& get_prior_precision() const {return prior_precision;}
      const ColumnVector& get_prior_mean() const {return prior_mean;}
      const ColumnVector& get_prior_isard() const {return prior_isard;}

      virtual void set_values() = 0;

      bool get_do_vb() const {return do_vb;}
      bool get_save_out() const {return save_out;}

    protected:
  
      string name;
      vector<string> parameter_names;

      int nparams;
      
      vector<vector<float> > samples; // num_params*num_samps

      vector<float> values;
      vector<float> param_min;
      vector<float> param_max;
      vector<float> gradient_step_sizes;

      Mvnormrandm normrand;

      int nacceptedtotal;
      int nrejectedtotal;   
      int naccepted;
      int nrejected;
      int naccepted_cov;
      int nrejected_cov;

      vector<vector<double> > ss_values;   // num_params*num_jumps_per_update_proposal

      float proposal_std;
      SymmetricMatrix proposal_cov;

      int cov_mode;

      int debuglevel;

      bool output_samples;

      Matrix jacobian;
      ColumnVector noise_precision;
      ColumnVector forward_model;

      ColumnVector prior_mean;
      SymmetricMatrix prior_precision;
      ColumnVector prior_isard;

      bool do_vb;

      bool save_out;

      Matrix jacob;
      ColumnVector fm;
      ColumnVector noise_prec;

    private:      
      Mcmc_Component();
      const Mcmc_Component& operator=(Mcmc_Component&);     
      Mcmc_Component(Mcmc_Component&);

      friend void copy_proposal_distributions(const vector<Mcmc_Component*>& from, vector<Mcmc_Component*> to);
    }; 

  class Mcmc_Log_Likelihood
  {
  public:
    // Constructor
    Mcmc_Log_Likelihood(){}
    
    virtual ~Mcmc_Log_Likelihood(){}
    
    virtual const double evaluate() = 0; // sets and returns log_likelihood

    virtual const double calc_energy() = 0; // sets and returns energy (log posterior)

    virtual void sample(){      
      log_likelihood_hist.push_back(log_likelihood);
    }

    const vector<double>& get_log_likelihood_hist() const {return log_likelihood_hist;}
    double get_log_likelihood() const {return log_likelihood;}

    virtual const int get_num_data_points() const = 0;

  protected:    
    vector<double> log_likelihood_hist;
    double log_likelihood;
    double energy;

  private:      
    const Mcmc_Log_Likelihood& operator=(Mcmc_Log_Likelihood&);     
    Mcmc_Log_Likelihood(Mcmc_Log_Likelihood&);
    
  };

  void copy_proposal_distributions(const vector<Mcmc_Component*>& from, vector<Mcmc_Component*> to);
}

#endif

