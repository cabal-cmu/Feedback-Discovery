/*  vb_component.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */
  
#if !defined(vb_component_h)
#define vb_component_h

#include "mcmc_component.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscprob.h"

using namespace Utilities;
using namespace MISCMATHS;

namespace Nma {

  class Vb_Component
    {
    public:
      // Constructor
      Vb_Component(const vector<Mcmc_Component*>& pmcmc_comps, int pdebuglevel);

      virtual ~Vb_Component(){}

      // setup
      void setup();
      
      // move to next point in VB
      double vb_jump(bool& accept);
      virtual double calc_prior_energy(const ColumnVector& pvalues, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision){
	ColumnVector tmp=pvalues-pprior_mean;
	return 0.5*(tmp.t()*pprior_precision*tmp).AsScalar();
      }

      virtual double calc_energy() = 0;
            
      void set_debuglevel(int pdebuglevel) {debuglevel=pdebuglevel;}
            
//       void append_gradient(vector<float>& grad_vec, double energy);
      
      void calc_jacobian(Matrix& jacob);

      virtual void calc_forward_model(ColumnVector& fm) = 0;
      virtual void calc_noise_precision(ColumnVector& prec) = 0;
      virtual const ColumnVector& get_vb_data() = 0;
    
      void initialise_vb();

//       // store calculations for old values when new values are proposed
//       virtual void store_old() = 0;
//       // restore calculations for old values when proposal is rejected
//       virtual void restore_old() = 0;

      //     void add_mcmc_component(Mcmc_Component* mcmc_comp) { mcmc_components.push_back(mcmc_comp);}
      
    protected:
  

      vector<Mcmc_Component*> mcmc_components;
      int nparams;
  
      Mvnormrandm normrand;

      int debuglevel;

      Matrix jacobian;
      ColumnVector noise_precision;
      ColumnVector forward_model;

      ColumnVector prior_mean;
      SymmetricMatrix prior_precision;
      ColumnVector prior_isard;

      Matrix jacob;
      ColumnVector fm;
      ColumnVector noise_prec;

    private:      
      Vb_Component();
      const Vb_Component& operator=(Vb_Component&);     
      Vb_Component(Vb_Component&);

    }; 


}

#endif

