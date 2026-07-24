/*  vb_component.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include "vb_component.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"
#include "newmat.h"
#include "utils/tracer_plus.h"

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWMAT;
using namespace MISCPLOT;
using namespace MISCPIC;

namespace Nma {

  void append_block(SymmetricMatrix& prior_precision, const SymmetricMatrix& add)
  {
    SymmetricMatrix tmp(prior_precision.Nrows()+add.Nrows());
    tmp=0;
    tmp.SymSubMatrix(1,prior_precision.Nrows())=prior_precision;
    tmp.SymSubMatrix(prior_precision.Nrows()+1,tmp.Nrows())=add;
    prior_precision=tmp;
  }

  Vb_Component::Vb_Component(const vector<Mcmc_Component*>& pmcmc_comps, int pdebuglevel)
    :         mcmc_components(pmcmc_comps),
	      debuglevel(pdebuglevel)
  {
    setup();
  }
  
  void Vb_Component::setup()
  {
    Tracer_Plus trace("Vb_Component::setup");

    //  setup ColumnVector prior_mean;

    // setup SymmetricMatrix prior_precision;

    nparams=0;
    for(unsigned int n=0; n<mcmc_components.size(); n++)
      {
	ColumnVector pm=mcmc_components[n]->get_prior_mean();
	prior_mean &= pm;
	append_block(prior_precision, mcmc_components[n]->get_prior_precision());
	ColumnVector pisard=mcmc_components[n]->get_prior_isard();
	prior_isard &= pisard;
	nparams+=pm.Nrows();
      }	
  }


  void Vb_Component::calc_jacobian(Matrix& jacobian)
  {
    Tracer_Plus trace("Vb_Component::calc_jacobian");

    calc_forward_model(fm);
    
    jacobian.ReSize(fm.Nrows(), nparams);
    int q=1;

    for(unsigned int n=0; n<mcmc_components.size(); n++)
      {
	mcmc_components[n]->store_old();

	vector<float> vals=mcmc_components[n]->get_values();
	vector<float> steps=mcmc_components[n]->get_gradient_step_sizes();
	for(unsigned int p=0; p<vals.size(); p++)
	  {
	    float old_val=vals[p];
	    vals[p]+=steps[p];
	    mcmc_components[n]->set_values_vec(vals);	  
	    mcmc_components[n]->set_values();

	    ColumnVector fm_grad;
	    
	    calc_forward_model(fm_grad);
	    jacobian.Column(q)=(fm_grad-fm)/steps[p];
	    q++;

	    vals[p]=old_val;	
	    mcmc_components[n]->set_values_vec(vals);	  
	  }


	mcmc_components[n]->set_values();	  
	mcmc_components[n]->restore_old();
      }

  }

  void Vb_Component::initialise_vb()
  {
    Tracer_Plus trace("Vb_Component::initialise_vb");

    // loop over comps extracting values
    ColumnVector values; 
    for(unsigned int n=0; n<mcmc_components.size(); n++)
      {
	ColumnVector vs=vector2ColumnVector(mcmc_components[n]->get_values());
	values &= vs;
      }

    for( int n=0; n<nparams; n++)
      {
	if(prior_isard(n+1))
	  {
	    prior_precision(n+1,n+1)=1.0/Sqr(values(n+1)-prior_mean(n+1));
	  }
      }
    calc_noise_precision(noise_prec);    
    //    OUT(size(noise_prec));
  }

  double Vb_Component::vb_jump(bool& accept)
  {
    Tracer_Plus trace("Vb_Component::vb_jump");
    
    accept=false;

    // backup old values
    vector<vector<float> > old_values_comps;
    for(unsigned int n=0; n<mcmc_components.size(); n++)
      {
	mcmc_components[n]->store_old();
	
	// save old values
	old_values_comps.push_back(mcmc_components[n]->get_values());
      }

    if(debuglevel>=6)
      LOGOUT("********************* start");
        
    // loop over comps extracting values
    ColumnVector old_values; 
    for(unsigned int n=0; n<mcmc_components.size(); n++)
      {
	ColumnVector vs=vector2ColumnVector(mcmc_components[n]->get_values());
	old_values &= vs;
	
// 	for(unsigned int p=0; p<mcmc_components[n]->get_values().size(); p++)
// 	  {
// 	    LogSingleton::getInstance().str() << mcmc_components[n]->get_parameter_names()[p] << "=" << mcmc_components[n]->get_values()[p] << endl;
// 	  }
      }
    

     if(debuglevel>=6)
       for(int n=0; n<nparams; n++)
	 LogSingleton::getInstance().str() << n << "_old=" << old_values(n+1) << endl;

    double old_lik_energy = calc_energy();         	
    double old_energy = old_lik_energy + calc_prior_energy(old_values,prior_mean,prior_precision);
 
    //     store_old();

    float maxfactor=1e10;
    float factor=1e-10;
    bool loop=true;
    // precompute stuff that can be precomputed

    // propose new parameters based on a nonlinear vb update
    //    OUT("Compute Jacobian");
    calc_jacobian(jacob);
    //   OUT(size(jacob));
    calc_forward_model(fm);
    //    OUT(size(fm));
    
    //Matrix tmp=jacob.t()*noise_prec_diag;
    Matrix tmp(jacob.Ncols(),noise_prec.Nrows());
    tmp=0;
    for( int r=1; r<=tmp.Nrows(); r++)
      for( int c=1; c<=tmp.Ncols(); c++)
	{
	  tmp(r,c)=jacob(c,r)*noise_prec(c);
	}
    
    SymmetricMatrix prec;
    //    prec<<jacob.t()*noise_prec*jacob+prior_precision;
    prec<<tmp*jacob+prior_precision;
    //    ColumnVector beta=jacob.t()*noise_prec*((get_vb_data()-fm)+jacob*vector2ColumnVector(old_values))+prior_precision*prior_mean;
    //    ColumnVector beta=tmp*((get_vb_data()-fm)+jacob*old_values)+prior_precision*prior_mean;
    ColumnVector beta2=tmp*(get_vb_data()-fm) + prior_precision*prior_mean;

    double new_lik_energy=0;         	
    double new_energy=0;
    ColumnVector values;

    while(loop)
      {	
	// update params

	// do L-M diagonal boosting	
	SymmetricMatrix prec_new=prec;
	for(int n=1; n<=prec.Nrows(); n++)      
	  {
	    prec_new(n,n)=prec(n,n)+factor*prec(n,n);
	  }

	//	values=prec_new.i()*beta;
	ColumnVector del=prec_new.i()*beta2;
	values=del+old_values;

	// test energy     
	loop=false;
	
	// loop over components inserting values into components
	int q=1;
	for(unsigned int n=0; n<mcmc_components.size(); n++)      
	  {	
	    // extract values for this component from updated vb values 
	    vector<float> vals;
	    for(unsigned int p=0; p<mcmc_components[n]->get_values().size(); p++)
	      {
		vals.push_back(values(q));
		q++;
	      }
	    
	    mcmc_components[n]->set_values_vec(vals);	  
	    mcmc_components[n]->set_values();
	  }

// 	// loop over components to see if values are in bounds
// 	for(unsigned int n=0; n<mcmc_components.size(); n++)      
// 	  {
// 	    if(!mcmc_components[n]->in_bounds())
// 	      {
// 		loop=true;
// 		OUT("OUT OF BOUNDS");
// 		for(unsigned int p=0; p<mcmc_components[n]->get_values().size(); p++)
// 		  {
// 		    LogSingleton::getInstance().str() << mcmc_components[n]->get_parameter_names()[p] << "=" << mcmc_components[n]->get_values()[p] << endl;
// 		  }
// 	      }
// 	  }
	
	if(debuglevel==6)
	  {    
	    //ColumnVector beta=tmp*((get_vb_data()-fm)+jacob*old_values)+prior_precision*prior_mean;

// 	    OUT(size(fm));
// 	    OUT(size(noise_prec));
// 	    OUT(size(jacob));

	    ColumnVector del=prec_new.i()*(jacob.t()*(SP(noise_prec,(get_vb_data()-fm))) + prior_precision*prior_mean);
	    OUT(del);
	    OUT(del+old_values);
	    OUT(values);

	    ColumnVector beta=tmp*((get_vb_data()-fm)+jacob*old_values)+prior_precision*prior_mean;
	    OUT(prec_new.i()*beta);
 
	    OUT((jacob.t()*(SP(noise_prec,(get_vb_data()-fm)))));
	    OUT(prec);
	    OUT(prec_new);
	    // 	    OUT(jacob.t()*noise_prec*((get_vb_data()-fm)+jacob*(old_values)));
	    OUT(prior_precision);
	    OUT(prior_mean);
	    OUT(beta2);

	    for(int n=0; n<nparams; n++)
	      LogSingleton::getInstance().str() << n << "_new=" << values(n+1) << endl;
	    for(int n=0; n<nparams; n++)
	      LogSingleton::getInstance().str() << n << "_old=" << old_values(n+1) << endl;
	    
	    for(unsigned int n=0; n<mcmc_components.size(); n++)
	      {
		OUT(mcmc_components[n]->get_values().size());
		OUT(vector2ColumnVector(mcmc_components[n]->get_values()));
		for(unsigned int p=0; p<mcmc_components[n]->get_values().size(); p++)
		  {
		    OUT(p);
		    LogSingleton::getInstance().str() << mcmc_components[n]->get_parameter_names()[p] << "=" << mcmc_components[n]->get_values()[p] << ", inbounds=" << mcmc_components[n]->in_bounds()<< endl;
		  }
	      }
	  }	

	if(!loop)
	  {
	    // calc new energy       	
	    new_lik_energy = calc_energy();         	
	    new_energy = new_lik_energy + calc_prior_energy(values,prior_mean,prior_precision);
	    	    
	    if(debuglevel==6)
	      {	
		OUT(old_lik_energy);
		OUT(new_lik_energy);
		LogSingleton::getInstance().str() << "old_energy()=" << old_energy << endl;
		LogSingleton::getInstance().str() << "new_energy()=" << new_energy << endl;
		LogSingleton::getInstance().str() << "old_energy()-new_energy()=" << old_energy-new_energy << endl;
		// 		for(int n=0; n<nparams; n++)
		// 		  LogSingleton::getInstance().str() << n << "_new=" << values(n+1) << endl;
		// 		LOGOUT(prec);
		// 		LOGOUT(prec_new);	  		
		
		// 	 write_ascii_matrix(get_vb_data(), LogSingleton::getInstance().appendDir("data"));	
		// 	 write_ascii_matrix(fm, LogSingleton::getInstance().appendDir("fm"));
		// 	 write_ascii_matrix(vector2ColumnVector(old_values), LogSingleton::getInstance().appendDir("old_values"));
		// 	 write_ascii_matrix(vector2ColumnVector(values), LogSingleton::getInstance().appendDir("values"));
		// 	 write_ascii_matrix(jacob, LogSingleton::getInstance().appendDir("jacob"));	    
	      }
	
	    if(new_energy<old_energy)
	      {
		accept=true;

		loop=false;		 
	      }
	    else 
	      {
		loop=true;

		for(unsigned int n=0; n<mcmc_components.size(); n++)
		  {
		    mcmc_components[n]->set_values_vec(old_values_comps[n]);
		    mcmc_components[n]->set_values();
		    mcmc_components[n]->restore_old();
		  }
	      }
	  }

	if(loop && factor<=maxfactor)
	  {
	    loop=true;
	    factor*=100;
	  }
	else
	  {
	    loop=false;
	  }

      }  // end of LM while loop

    // update ARD precisions
    SymmetricMatrix covar=prec.i();
    for( int n=0; n<nparams; n++)
	  {
	    if(prior_isard(n+1))
	      {
		prior_precision(n+1,n+1)=1.0/(Sqr(values(n+1)-prior_mean(n+1))+covar(n+1,n+1));
	      }
	  }
    
    // update noise precisions
    calc_noise_precision(noise_prec);    

    if(accept)
      {
	//	if(debuglevel==6)
	  {
	    OUT("ACCEPTED");
	    OUT(factor);
	  }
	old_lik_energy = new_lik_energy;
      }
    else // reject      
      {      
	//	if(debuglevel==6)
	  {
	    OUT("REJECTED");
	    OUT(factor);       	
	  }
      }
      
    if(debuglevel==6)
      LOGOUT("********************* end");
		
  return old_lik_energy;
  
}

 

}
