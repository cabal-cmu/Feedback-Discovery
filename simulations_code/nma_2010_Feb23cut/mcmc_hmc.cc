/*  mcmc_hmc.cc

    Mark Woolrich FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  CCOPYRIGHT  */

#include "mcmc_hmc.h"
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

//   void Mcmc_hmc::update_proposal_stds()
//   {
//     Tracer_Plus trace("Mcmc_hmc::update_proposal_stds");
    
//     for(int p=1; p<=nparams; p++)
//       {
// 	//    double rejectionrate=0.5;   
// 	//    proposal_std *= rejectionrate/((1+nrejected)/double(1+naccepted+nrejected));
	
// 	// F=10; N=800; x=log(F)/log((N-1)/1); clear rs update;  rs=0:N; for i=1:length(rs), r=rs(i); a=N-r; if(a==0) a=1, r=N-1; end; if(r==0) r=1, a=N-1; end; update(i)=(a/r)^x; end; figure; semilogy(rs,update); for i=1:length(rs), r=rs(i); a=N-r; update(i)=0.5*(1+a+r)/(1+r); end; ho; semilogy(rs,update);
	
// 	int a=naccepted(p);
// 	int r=nrejected(p);

// 	float F=10; // this will be the update if nrejected=0; update will be 1/F if naccepted=0
	
// 	if(a+r<=2) {a=1; r=100;} // this will be because everything is out of bounds
// 	if(a==0) {a=1; r=r-1;}
// 	if(r==0) {r=1; a=a-1;}
// 	int N=a+r;        
	
// 	float x=std::log(F)/std::log((N-1.0)/1.0);
// 	float update=std::pow(float(a)/float(r),x);
// 	proposal_std(p) *= update;    
	
// 	if(proposal_std<1e-10) proposal_std(p)=1e-10;
	
// 	if(debuglevel==4)
// 	  {
// 	    LOGOUT("Update proposal_std");
// 	    LOGOUT(parameter_names[p-1]);  
// 	    LOGOUT(update);	   
// 	    LOGOUT(naccepted(p));
// 	    LOGOUT(nrejected(p));
// 	    LOGOUT(float(a)/float(r));
// 	    LOGOUT(in_bounds());
// 	    LOGOUT(proposal_std);
// 	    LOGOUT("End Update proposal_std");
// 	  }
	
// 	naccepted(p) = 0;
// 	nrejected(p) = 0;	
	
//       }

//   }

  void Mcmc_hmc::initialise_epsilons()
  {
    Tracer_Plus trace("Mcmc_hmc::initialise_epsilons");
  
    LOGOUT("Initialise epsilons");
       
    epsilons.ReSize(nparams);
    epsilons=0.001;

    // save old values
    RowVector old_values=values;
    RowVector values_tmo;
    RowVector values_tmt;

    // initial momentum is Normal(0,1) 
    ColumnVector momentum;    
    ColumnVector momentum_tmo;
    ColumnVector momentum_tmt;

    for(int i=0; i<ncomponents; i++)
      components[i]->store_old();

    // calc energy
    double energy = calc_energy();   

    // calc gradient
    ColumnVector gradient;
    calc_gradient(gradient, energy);
    ColumnVector old_gradient=gradient;

    float tol=0.01;
    float lam=0.01;

    bool finished=false;

    while(!finished)
      {
	// restore values
	ColumnVector momentum=normrnd(nparams,1,0,1);  
	gradient=old_gradient;
	values=old_values;
	values_tmo=values;
	momentum_tmo=momentum;
	
	// do two values steps to give estimates of values_tmo and values_tmt
	for(int i=1; i<=2; i++)
	  {
	    values_tmt=values_tmo;
	    values_tmo=values;
	    momentum_tmt=momentum_tmo;
	    momentum_tmo=momentum;
	    
	    momentum -=  SP(epsilons,gradient) / 2.0 ; // make half-step in momentum
	    values += SP(epsilons,momentum).t() ; // make step in x 
	    set_values(values); // make sure model object updates its parameters	   
	    calc_gradient(gradient, calc_energy()); // find new gradient 
	    momentum -= SP(epsilons,gradient) / 2.0 ; // make half-step in momentum 
	  }

	vector<bool> completed(nparams);
	for(int n=1; n<=nparams; n++)
	  {
	    // compare one step with two to assess epsilon
	    float values_twostep=values_tmt(n)+2*epsilons(n)*momentum_tmt(n);	   
	    float err=std::abs(values_twostep-values(n))/std::abs(values(n));
	    
	    // update eps so that err is tending towards tol   	 
	    float update=std::pow(tol/err,lam);
	    epsilons(n)=update*epsilons(n);
	    
	    if(debuglevel==3 )
	      {
		LOGOUT(err);
		LOGOUT(n);
		LogSingleton::getInstance().str() << "update=" << update << endl;       
		LOGOUT(epsilons(n));
		OUT(vector2ColumnVector(completed));
	      }
	   
	    completed[n-1] = (err<10*tol && err>0.1*tol);
	  }

	// test to see if finished - only finished if epsilons for all params are within acceptable tolerance range
	
	if(!finished)
	  {
	    finished=true;
	    for(int n=1; n<=nparams; n++)
	      {		
		if(!completed[n-1])
		  finished=false;
	      }	  
	  }
      }
    
    // restore stuff
    values=old_values;
    set_values(values);
    
    for(int i=0; i<ncomponents; i++)
      components[i]->restore_old();    
    
    LOGOUT(epsilons);
    LOGOUT("Finished");
  }

  void Mcmc_hmc::jump(double& energy)
  {
    Tracer_Plus trace("Mcmc_hmc::jump");
    
    if(debuglevel==3 )
      {	
	LOGOUT("********************* start");
      } 
    // save old values
    RowVector old_values=values;
    RowVector values_tmo;
    RowVector values_tmt;

    for(int i=0; i<ncomponents; i++)
      components[i]->store_old();

    // initial momentum is Normal(0,1) 
    ColumnVector momentum=SP(epsilons/(mean(epsilons).AsScalar()),normrnd(nparams,1,0,1));    
    ColumnVector momentum_tmo;
    ColumnVector momentum_tmt;

    // calc energy
    energy = calc_energy();
    double old_energy = energy;
 
    // calc gradient
    ColumnVector gradient;
    calc_gradient(gradient, energy);
    ColumnVector old_gradient=gradient;

    // evaluate H(x,p)  
    double kinetic_energy=(momentum.t()*momentum).AsScalar()/2.0;  
    double H = kinetic_energy + energy;

    int nsubjumps=int(unifrnd(1,1,5,15).AsScalar());
    
    bool inbounds=true;

    values_tmo=values;
    momentum_tmo=momentum;

    for(int i=1; i<=nsubjumps && inbounds; i++)
      {
	values_tmt=values_tmo;
	values_tmo=values;
	momentum_tmt=momentum_tmo;
	momentum_tmo=momentum;

	momentum -=  SP(epsilons,gradient) / 2.0 ; // make half-step in momentum
	values += SP(epsilons,momentum).t() ; // make step in values 
	set_values(values); // make sure model object updates its parameters

	if(!in_bounds()) 
	  {
	    inbounds=false; // check if in bounds
	    values=values_tmo;
	    set_values(values); // make sure model object updates its parameters
	    momentum=momentum_tmo;
	    
	    if(debuglevel==3 )
	      {
		LogSingleton::getInstance().str() << "Out of bounds." << endl;
		LogSingleton::getInstance().str() << "i=" << i << endl;
	      }
	    
	    break;
	  }

	calc_gradient(gradient, calc_energy()); // find new gradient 
	momentum -= SP(epsilons,gradient) / 2.0 ; // make half-step in momentum 

	if(i>3)
	  {
	    for(int n=1; n<=nparams; n++)
	      {
		/////// compare one step with two to assess epsilon
		float values_twostep=values_tmt(n)+2*epsilons(n)*momentum_tmt(n);
		
		float err=std::abs(values_twostep-values(n))/std::abs(values(n));
		
		////// update eps so that err is tending towards tol   
		
		float tol=0.01;
		float lam=0.01;
		float update=std::pow(tol/err,lam);
		epsilons(n)=update*epsilons(n);
		
// 		if(debuglevel==3 )
// 		  {
// 		    LOGOUT(err);
// 		    LOGOUT(n);
// 		    LogSingleton::getInstance().str() << "update=" << update << endl;       
// 		    LOGOUT(epsilons(n));
// 		  } 
	      }
	  }
      }

    bool accept=false;

    energy = calc_energy();   
    double kinetic_energy_new=(momentum.t()*momentum).AsScalar()/2.0;  
    double Hnew = kinetic_energy_new + energy;
    double dH = Hnew - H ; 
    
    if ( dH < 0 ) accept = true ; 
    else if ( unifrnd().AsScalar() < std::exp(-dH) ) accept = true ; 
    else accept = false ;
    
    if(debuglevel==3 )
      {
	LogSingleton::getInstance().str() << "new_H()=" << Hnew << endl;
	LogSingleton::getInstance().str() << "new_kinetic_energy()=" << kinetic_energy_new << endl;
	LogSingleton::getInstance().str() << "new_energy()=" << energy << endl;
      } 
                 
    if(debuglevel==3 )
      {
	for(int n=0; n<nparams; n++)
	  LogSingleton::getInstance().str() << parameter_names[n] << "=" << values(n+1) << endl;
	LogSingleton::getInstance().str() << "old_H()=" << H << endl;
	LogSingleton::getInstance().str() << "old_kinetic_energy()=" << kinetic_energy << endl;
	LogSingleton::getInstance().str() << "old_energy()=" << old_energy << endl;
	for(int n=0; n<nparams; n++)
	  LogSingleton::getInstance().str() << parameter_names[n] << " (old) =" << old_values(n+1) << endl;
	LogSingleton::getInstance().str() << "accept=" << accept << endl;
	LOGOUT(in_bounds());
	
	LOGOUT("********************* end");
      } 
  
    if(accept)
      {		  	 
	naccepted++;

	bool nochange=true;
	for(int n=0; n<nparams; n++)
	  {
	    double del=values(n+1)-old_values(n+1);
	    if(del!=0) nochange=false;
	  }

	if(nochange && inbounds)
	  {
	    LOGOUT("no change, but the jump is being accepted");
	    
// 	    LOGOUT(epsilons);

	    LogSingleton::getInstance().str() << "old_energy()=" << old_energy << endl;
	    LogSingleton::getInstance().str() << "new_energy()=" << energy << endl;

	  }      

      }
    else // reject      
      {
	if(!isnan(energy))
	  {
	    nrejected++;
	  }

	values=old_values;
	set_values(values);

	for(int i=0; i<ncomponents; i++)
	  components[i]->restore_old();

	energy = old_energy;
      }

  }

  double Mcmc_hmc::calc_energy()
  {
    Tracer_Plus trace("Mcmc_hmc::calc_energy");
    
    return mcmc_log_likelihood.calc_energy();
  }

  void Mcmc_hmc::calc_gradient(ColumnVector& gradient, double energy)
  {
    Tracer_Plus trace("Mcmc_hmc::calc_gradient");
    
    vector<float> grad_vec;
    for(int i=0; i<ncomponents; i++)
      {
	//	OUT(components[i]->get_name());
	components[i]->append_gradient(grad_vec, energy);
      }
    gradient=vector2ColumnVector(grad_vec);

  }

  bool Mcmc_hmc::in_bounds()
  {
    Tracer_Plus trace("Mcmc_hmc::in_bounds");

    bool inbounds=true;
    for(int i=0; i<ncomponents && inbounds; i++)
      {
	inbounds=components[i]->in_bounds();

	if(!inbounds)
	  if(debuglevel==3)
	    {
	      LogSingleton::getInstance().str() << "Out of bounds." << endl;
	      LogSingleton::getInstance().str() << "component=" << components[i]->get_name() << endl;
	    }
      }
    return inbounds;
  }


  void Mcmc_hmc::set_values(const RowVector& values)
  {
    Tracer_Plus trace("Mcmc_hmc::set_values");

    int index=1;
    for(int i=0; i<ncomponents; i++)
      {
	vector<float> comp_values;

	for(unsigned int j=0; j<components[i]->get_parameter_names().size(); j++)
	  {	    
	    comp_values.push_back(values(index));
	    index++;
	  }

	components[i]->set_values_vec(comp_values);
	components[i]->set_values(); 
      }

  }

  void Mcmc_hmc::setup()
  {
    Tracer_Plus trace("Mcmc_hmc::setup");

    naccepted =0;
    nrejected =0;

    vector<float> values_vec;

    for(int i=0; i<ncomponents; i++)
      {
	vector<string> names=components[i]->get_parameter_names();
	const vector<float>& comp_values=components[i]->get_values();
	for(unsigned int j=0; j<names.size(); j++)
	  {
	    values_vec.push_back(comp_values[j]);
	    parameter_names.push_back(names[j]);
	  }
      }
    values=vector2ColumnVector(values_vec).t();

    nparams=parameter_names.size();

    // set epsilons 
    initialise_epsilons();

//     int index=1;
//     for(int i=0; i<ncomponents; i++)
//       {
// 	vector<string> names=components[i]->get_parameter_names();
// 	for(int j=0; j<names.size(); j++)
// 	  {
// 	    epsilons(index)=(components[i]->get_param_max()[j]-components[i]->get_param_min()[j])/100.0;
// 	    index++;
// 	  }
//       }

  }

  void Mcmc_hmc::run()
  {
    Tracer_Plus trace("Mcmc_hmc::run");

    OUT(debuglevel);

//     int jumps_per_proposal_update=200; //300
//     int update_count=0;

//     if(nburnin==0)
//       {
// 	// take samples		
// 	for(int i=0; i<ncomponents; i++)
// 	  components[i]->sample();
//       }

    double energy;

    int count_out=0;
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
	    jump(energy);
	    
	   //  if(update_count>jumps_per_proposal_update)
// 	      {
// 		if(debuglevel==2)
// 		  {
// 		    LOGOUT("UPDATE PROPOSAL STDS");
// 		    LOGOUT(n); LOGOUT(k);
// 		    LOGOUT(energy);
// 		  }
		
// 		// update proposal distribution stds
// 		for(int i=0; i<ncomponents; i++)
// 		  {
// 		    update_proposal_stds();				
// 		  }
// 		update_count=0;
// 	      }
// 	    else
// 	      update_count++;	
	    
	    
	  }	  
	
	if(n>nburnin)
	  {
	    // take samples		
	    for(int i=0; i<ncomponents; i++)
	      components[i]->sample();

	    energy_hist.push_back(energy);

	    // sample log likelihood - will be used to calc evidence
	    mcmc_log_likelihood.evaluate();
	    mcmc_log_likelihood.sample();
	  }	
	
      }

    LogSingleton::getInstance().str() << endl;

    OUT(naccepted);
    OUT(nrejected);
    
  }

}
