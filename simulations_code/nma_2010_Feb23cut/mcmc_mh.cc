/*  mcmc_mh.cc

    Mark Woolrich FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  CCOPYRIGHT  */

#include "mcmc_mh.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include <set>

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace Nma {

  ///////////////////////
  /// Mcmc fns

  void Mcmc::save()
  {
    Tracer_Plus trace("Mcmc::save");

    evaluate_model_evidence();

    write_num(get_unnormalised_model_evidence(),LogSingleton::getInstance().appendDir("unnormalised_model_evidence.txt"));
    write_num(get_model_evidence_normalisation(),LogSingleton::getInstance().appendDir("model_evidence_normalisation.txt"));
    write_num(get_aic(),LogSingleton::getInstance().appendDir("aic.txt"));
    write_num(get_bic(),LogSingleton::getInstance().appendDir("bic.txt"));

    write_ascii_matrix(vector2ColumnVector(energy_hist),LogSingleton::getInstance().appendDir("energy_hist.txt"));
    write_ascii_matrix(vector2ColumnVector(mcmc_log_likelihood.get_log_likelihood_hist()),LogSingleton::getInstance().appendDir("log_likelihood_hist.txt"));

    for(int i=0; i<ncomponents; i++)
      components[i]->save();
  }
  
  void Mcmc::set_values_to_sample_mean()
  {
    Tracer_Plus trace("Mcmc::set_values_to_sample_mean");
    
    for(int i=0; i<ncomponents; i++)
      components[i]->set_values_to_sample_mean();
  }

  void Mcmc::set_values_to_sample_map()
  {
    Tracer_Plus trace("Mcmc::set_values_to_sample_map");
    
    // first: find the MAP 
    int map_index=0;
    double min_en=1e32;
    for(unsigned int i=0; i<energy_hist.size(); i++)
      {
	if(energy_hist[i]<min_en)
	  {
	    map_index=i;
	    min_en=energy_hist[i];
	  }
      }

    for(int i=0; i<ncomponents; i++)      
      components[i]->set_values_to_sample(map_index+1);
  }

  void Mcmc::find_parameter(const string& name, int& component_number, int& parameter_number) const
  {    
    for(unsigned int i=0; i<components.size(); i++)
      {
	const vector<string>& parameter_names=components[i]->get_parameter_names();
	
	for(unsigned int n=0; n< parameter_names.size(); n++)
	  {
	    if(name==parameter_names[n])
	      {
		component_number=i;
		parameter_number=n;
	      }
	  }
      }
  }


  void Mcmc::evaluate_model_evidence()
  {
    Tracer_Plus trace("Mcmc::evaluate_model_evidence");
    
    // uses mcmc samples approach and also calcs AIC, BIC

    ColumnVector logliktmp=vector2ColumnVector(mcmc_log_likelihood.get_log_likelihood_hist());
    float minloglik=logliktmp.Minimum();

    minloglik+=100;

    // eqn 14 of: A Bayesian framework for global tractography, Jbabdi et al., NI (2007)
    // P(Y)=1/N ( \sum_s 1/P(Y|theta_s))^(-1) where s indexes posterior samples and P(Y|theta) is the likelihood
    //
    // To avoid overflow we actually use P(Y|theta_s)=exp(log(P(Y|theta_s))-model_evidence_normalisation)*exp(model_evidence_normalisation) 
    // giving P(Y)=exp(model_evidence_normalisation)*unnormalised_model_evidence=exp(MEN)*UME
    // where unnormalised_model_evidence = 1/N (\sum_s exp[-(log(P(Y|theta_s))-model_evidence_normalisation)])^(-1)
    //
    // so bayes factor between M1 and M2 is: UME_1/UME_2*exp(MEN_1-MEN_2)

    model_evidence_normalisation=minloglik;

    double evidence_sum=0;

    int count=0;
    for(int i=1; i<=logliktmp.Nrows(); i++)
      {       
	double tmp=(logliktmp(i)-model_evidence_normalisation);
	double tmp2=std::exp(-tmp);
	if(tmp2>0) 
	  {
// 	    LOGOUT(i);
// 	    LOGOUT(tmp);
// 	    LOGOUT(tmp2);
// 	    LOGOUT(minloglik);
// 	    LOGOUT(logliktmp(i));
// 	    LOGOUT(evidence_sum);
// 	    if(isnan(tmp2))
// 	      exit(1);

 	    count++;	    
	    evidence_sum += tmp2;
	  }

      }

    unnormalised_model_evidence=1.0/(evidence_sum/nsamps);
    
    LOGOUT(evidence_sum);
    LOGOUT(nsamps);
    LOGOUT(unnormalised_model_evidence);
    LOGOUT(model_evidence_normalisation);

    float percent_samps_with_nonzero_evidence=100*float(count)/float(nsamps);

    LOGOUT(percent_samps_with_nonzero_evidence);

    // now get MAP log likelihood to do AIC and BIC calculations:
    
    // first: find the MAP 
    int map_index=0;
    double min_en=1e32;
    for(unsigned int i=0; i<energy_hist.size(); i++)
      {
	if(energy_hist[i]<min_en)
	  {
	    map_index=i;
	    min_en=energy_hist[i];
	  }
      }
    
    float loglik=mcmc_log_likelihood.get_log_likelihood_hist()[map_index];

    // calc num params
    int nparams=0;
    for(int i=0; i<ncomponents; i++)
      nparams += components[i]->get_values().size();

    // num of data points
    int ndata_points=mcmc_log_likelihood.get_num_data_points();
    
    LOGOUT(nparams);
    LOGOUT(ndata_points);

    // aic
    float complexity_term=nparams*2;
    aic=-2*loglik+complexity_term;

    LOGOUT(aic);

    // bic
    complexity_term=nparams*std::log(ndata_points);
    bic=-2*loglik+complexity_term;

    LOGOUT(bic);
  }
  
  ////////////////////////
  // Mcmc_Mh fns

  void Mcmc_Mh::setup()
  {
    Tracer_Plus trace("Mcmc_Mh::setup");

  }

  void Mcmc_Mh::run()
  {
    Tracer_Plus trace("Mcmc_Mh::run");

    int jumps_per_proposal_update=200; //300
    // int jumps_per_vb_update=30;
    int jumps_per_cov_update=500;	//900
    int update_count=0;
    int update_cov_count=0;

    double energy_lik=components[0]->calc_energy();

    if(debuglevel==2)
      {
	for(int i=0; i<ncomponents; i++)
	  {
	    OUT(components[i]->get_name());
	    OUT(components[i]->get_proposal_std());
	    OUT(components[i]->get_proposal_cov());
	  }
      }

//     if(nburnin==0)
//       {
// 	// take samples		
// 	for(int i=0; i<ncomponents; i++)
// 	  components[i]->sample();
//       }

    int count_out=0;
    int count_vb=sample_vb_jump_every-1;

    int max_vb_iterations=50;
    
    bool accept=true;
    
    vb_component.initialise_vb();
    for(int n=1; n<=max_vb_iterations && accept; n++)    
      {
	energy_lik=vb_component.vb_jump(accept);
      }

    for(int n=1; n<=(nsamps+nburnin); n++)
      {
	if(n > count_out*(nsamps+nburnin)/100.0)
	  {
	    LogSingleton::getInstance().str() << " " << (count_out+1);
	    LogSingleton::getInstance().flush();
	    count_out++;
	  }

	for(int k=1; k<=njumps_per_sample; k++)
	  {
	    for(int i=0; i<ncomponents; i++)
	      energy_lik=components[i]->jump(energy_lik);

	    // THIS WILL NOT WORK WITH DECODING
	    //	    if(count_vb>jumps_per_vb_update && n<nburnin)
	    //	    if((count_vb>burnin_vb_jump_every && n<nburnin && burnin_vb_jump_every!=-1) || (count_vb>sample_vb_jump_every && n>=nburnin && sample_vb_jump_every!=-1))
	    if((count_vb>burnin_vb_jump_every && n<nburnin && burnin_vb_jump_every!=-1) || (count_vb>sample_vb_jump_every && n>=nburnin && sample_vb_jump_every!=-1))
	      {	    
		bool accept=true;

		vb_component.initialise_vb();
		for(int n=1; n<=max_vb_iterations && accept; n++)    
		  {
		    energy_lik=vb_component.vb_jump(accept);
		  }

// 		for(int i=0; i<ncomponents; i++)
// 		  {
// 		    if(components[i]->get_do_vb())
// 		      energy_lik=components[i]->vb_jump(energy);
// 		  }
		count_vb=0;
	      }
	    else
	      count_vb++;

	    //if(n<nburnin || sample_vb_jump_every!=-1)
	      {
		if(update_count>jumps_per_proposal_update)
		  {
		    if(debuglevel==2)
		    {
		      LOGOUT("UPDATE PROPOSAL STDS");
		      LOGOUT(n); LOGOUT(k);
		      LOGOUT(energy_lik);
		    }
		    
		    // update proposal distribution stds
		    for(int i=0; i<ncomponents; i++)
		      {
			components[i]->update_proposal_stds();				
		      }
		    update_count=0;
		  }
		else
		  update_count++;		
	      	    
		
		if(update_cov_count>jumps_per_cov_update)
		  {
		    if(debuglevel==2)
		      {
			LOGOUT("UPDATE PROPOSAL COV");
			LOGOUT(n); LOGOUT(k);
			LOGOUT(energy_lik);
		      }
		    
		    // update cov shape
		for(int i=0; i<ncomponents; i++)
		  {
		    components[i]->update_proposal_cov(jumps_per_cov_update);		
		    if(debuglevel==2) LOGOUT("COMPONENT"); 
		  }
		
		update_cov_count=0;
		  }
		else
		  update_cov_count++;
	      }	  
	  }
	
	if(n>nburnin)
	  {
	    // take samples		
	    for(int i=0; i<ncomponents; i++)
	      components[i]->sample();

	    energy_hist.push_back(energy_lik);

	    // sample log likelihood - will be used to calc evidence
	    mcmc_log_likelihood.evaluate();
	    mcmc_log_likelihood.sample();
	  }	
	
      }
    LogSingleton::getInstance().str() << endl;
    
  }



}
