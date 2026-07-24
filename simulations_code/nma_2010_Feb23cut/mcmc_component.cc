/*  mcmc_component.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include "mcmc_component.h"
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

  Mcmc_Component::Mcmc_Component(const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pcov_mode, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
    :name(pname),
     parameter_names(pparameter_names),     
     nparams(pparameter_names.size()),
     samples(pparameter_names.size()),
     values(pinitial_values), 
     param_min(pparam_min),
     param_max(pparam_max),
     gradient_step_sizes(pgradient_step_sizes),
     ss_values(pparameter_names.size()),
     proposal_std(1),
     proposal_cov(pparameter_names.size()),
     cov_mode(pcov_mode),
     debuglevel(pdebuglevel),
     output_samples(poutput_samples),
     prior_mean(pprior_mean),
     prior_precision(pprior_precision),
     do_vb(true),
     save_out(true)
  {
    prior_isard=prior_mean; 
    prior_isard=0;
  }
  
  Mcmc_Component::Mcmc_Component(const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pcov_mode, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision, const ColumnVector& pprior_isard)
    :name(pname),
     parameter_names(pparameter_names),     
     nparams(pparameter_names.size()),
     samples(pparameter_names.size()),
     values(pinitial_values), 
     param_min(pparam_min),
     param_max(pparam_max),
     gradient_step_sizes(pgradient_step_sizes),
     ss_values(pparameter_names.size()),
     proposal_std(1),
     proposal_cov(pparameter_names.size()),
     cov_mode(pcov_mode),
     debuglevel(pdebuglevel),
     output_samples(poutput_samples),
     prior_mean(pprior_mean),
     prior_precision(pprior_precision),
     prior_isard(pprior_isard),
     do_vb(true),
     save_out(true)
  {
  }
  
  void Mcmc_Component::setup()
  {
    Tracer_Plus trace("Mcmc_Component::setup");

    proposal_cov=0;

    for(int n=0; n<nparams; n++)
      {
	if(debuglevel==2)
	  {
	    LOGOUT(param_max[n]);
	    LOGOUT(param_min[n]);
	  }

	if(param_max[n]==param_min[n])
	  {
	    proposal_cov(n+1,n+1)=1;
	    LOGOUT(parameter_names[n]);
	    LOGOUT(param_max[n]);
	    LOGOUT(param_min[n]);
	    LOGOUT("WARNING (param_max[n]==param_min[n]) in Mcmc_Component::setup");
	  }
	else
	  proposal_cov(n+1,n+1)=(param_max[n]-param_min[n])/(100.0*Sqr(nparams));

// 	OUT(parameter_names[n]);
// 	OUT(nparams);
// 	OUT((param_max[n]-param_min[n])/(100.0*Sqr(nparams)));
     }
    
    if(cov_mode==3)
      {
	proposal_cov<<IdentityMatrix(nparams);
      }

    double ssnormalise=1.0;
    for(int n=0; n<nparams; n++)
      {
	ssnormalise*=pow(proposal_cov(n+1,n+1),1.0/nparams);
      }

    proposal_std=ssnormalise;
   
    ssnormalise=1.0;//sqrt(prod(diag(cov)));
    for(int k=1; k<=proposal_cov.Nrows(); k++)
      {
	ssnormalise*=pow(proposal_cov(k,k),1.0/nparams);
      }
    proposal_cov=(proposal_cov/ssnormalise);	     	      
    
//     if(debuglevel==5)// && parameter_names[0]=="ipl_art_hrf_2")
//       {
// 	LOGOUT("In Mcmc_Component::setup");
// 	LOGOUT(parameter_names[0]);    
// 	LOGOUT(ssnormalise);
// 	LOGOUT(proposal_std);
// 	LOGOUT(proposal_cov);
//       }

    normrand.setcovar(proposal_cov*Sqr(proposal_std));

    // setup counters
    nacceptedtotal = 0; nrejectedtotal = 0;
    naccepted = 0; nrejected = 0;
    naccepted_cov = 0; nrejected_cov = 0;

    ss_values.resize(nparams);

    // initialise energy    
    set_values();
    calc_energy();

//     do_vb=false;
//     if(do_vb)
//       {
// 	// setup vb stuff
// 	//data   
// 	setup_data();
//       }
  }

  void Mcmc_Component::sample()
  {
    Tracer_Plus trace("Mcmc_Component::sample");
      
    for(int n=0; n<nparams; n++)
      samples[n].push_back(values[n]);
  }

  void Mcmc_Component::append_gradient(vector<float>& grad_vec, double energy_lik)
  {
    Tracer_Plus trace("Mcmc_Component::append_gradient");
    
    float step=1e-6;
    store_old();
    double energy=energy_lik+calc_prior_energy(vector2ColumnVector(values),prior_mean,prior_precision);

    for(int n=0; n<nparams; n++)
      {
	float old_val=values[n];
	values[n]+=step;
	set_values();	  

	float energy_grad=calc_energy()+calc_prior_energy(vector2ColumnVector(values),prior_mean,prior_precision);

	grad_vec.push_back((energy_grad-energy)/step);

	values[n]=old_val;
      }

    set_values();	  
    restore_old();
  }

  void Mcmc_Component::calc_jacobian(Matrix& jacobian)
  {
    Tracer_Plus trace("Mcmc_Component::calc_jacobian");

    float step=1e-4;
    store_old();
    calc_forward_model(fm);
    
    jacobian.ReSize(fm.Nrows(), nparams);

    for(int n=0; n<nparams; n++)
      {
	float old_val=values[n];
	values[n]+=step;
	set_values();	  

	ColumnVector fm_grad;
	
	calc_forward_model(fm_grad);
	jacobian.Column(n+1)=(fm_grad-fm)/step;

	values[n]=old_val;	
      }

    set_values();	  
    restore_old();

  }


  void Mcmc_Component::update_proposal_stds()
  {
    Tracer_Plus trace("Mcmc_Component::update_proposal_stds");
   
    //    double rejectionrate=0.5;   
    //    proposal_std *= rejectionrate/((1+nrejected)/double(1+naccepted+nrejected));

    // F=10; N=800; x=log(F)/log((N-1)/1); clear rs update;  rs=0:N; for i=1:length(rs), r=rs(i); a=N-r; if(a==0) a=1, r=N-1; end; if(r==0) r=1, a=N-a; end; update(i)=(a/r)^x; end; figure; semilogy(rs,update); for i=1:length(rs), r=rs(i); a=N-r; update(i)=0.5*(1+a+r)/(1+r); end; ho; semilogy(rs,update);

    int a=naccepted;
    int r=nrejected;
    float F=10; // this will be the update if nrejected=0; update will be 1/F if naccepted=0

    if(a+r<=2) {a=1; r=100;} // this will be because everything is out of bounds
    if(a==0) {a=1; r=r-1;}
    if(r==0) {r=1; a=a-1;}
    int N=a+r;        

    float x=std::log(F)/std::log((N-1.0)/1.0);
    float update=std::pow(float(a)/float(r),x);
    proposal_std *= update;    

    if(proposal_std<1e-10) proposal_std=1e-10;

    if(debuglevel==4)// && parameter_names[0]=="ipl_art_hrf_2")
      {
	LOGOUT("Update proposal_std");
	LOGOUT(parameter_names[0]);  
	LOGOUT(update);
	LOGOUT(x);
	
	LOGOUT(naccepted);
	LOGOUT(nrejected);
	LOGOUT(float(a)/float(r));
	LOGOUT(in_bounds());
	for(int n=0; n<nparams; n++)
	  LogSingleton::getInstance().str() << parameter_names[n] << "=" << values[n] << endl;

	LOGOUT(proposal_std);
	LOGOUT(proposal_cov);
	LOGOUT("End Update proposal_std");
      }

    try{
      normrand.setcovar(proposal_cov*Sqr(proposal_std));	
    }
    catch(Exception& cock)
      {
	LOGOUT("Update proposal");
	LOGOUT(parameter_names[0]);    
	LOGOUT(naccepted);
	LOGOUT(nrejected);
	LOGOUT(proposal_std);
	LOGOUT(proposal_cov*Sqr(proposal_std));
	exit(1);
      }
 
    naccepted = 0;
    nrejected = 0;	
  }

  void Mcmc_Component::update_proposal_cov(int njumps_per_update) 
  {
    Tracer_Plus("Mcmc_Component::update_proposal_cov");       

    if(debuglevel==4)// && parameter_names[0]=="ipl_art_hrf_2")
      {   
	LOGOUT("Update proposal_cov");
	LOGOUT(parameter_names[0]);   
	// 	LOGOUT(proposal_std);
	// 	LOGOUT(proposal_cov);
	LOGOUT(naccepted_cov);
	LOGOUT(nrejected_cov);
      }
    
    SymmetricMatrix proposal_cov_old(proposal_cov);
    
    // 0 for full cov update
    // 1 for same corr for all
    // 2 for no corr update but individual variance update
    // 3 for no corr and no individual variance update
    //bool cov_mode_local=cov_mode; 

    //int max_params=8;

//     if(nparams>=max_params && !global_cov) 
//       {
// 	LOGOUT("Warning. Too many params in MCMC component to update proposal correlation. Will use same correlation for all params in the MCMC component.");
// 	global_cov_local=true;
//       }

    int subsample_num=5;

    if(naccepted_cov/float(subsample_num)-nparams>10)// && nrejected_cov>njumps_per_update/5.0 && nparams>1)
      {	  
	ColumnVector mns(nparams);
	mns=0;
	vector<ColumnVector> ss_values_subsampled(nparams);

	for(int n=0; n<nparams; n++)
	  {
	    // subsample ss_values

	    vector<double> ss_values_subsampled_vec;
 
	    int tmpcount=0;
	    for(int k=0; k<int(ss_values[n].size()); k++)
	      {
		if(tmpcount>subsample_num)
		  {
		    ss_values_subsampled_vec.push_back(ss_values[n][k]);
		    
		    tmpcount=0;
		  }
		else
		  {
		    tmpcount++;
		  }
	      }

	    // store in vector 
	    ss_values_subsampled[n]=vector2ColumnVector(ss_values_subsampled_vec);

	    // calc means
	    mns(n+1)=mean(ss_values_subsampled[n]).AsScalar();
	  }

	proposal_cov.ReSize(nparams);
	proposal_cov=0;

	int N=ss_values_subsampled[0].Nrows();

	for(int k=1; k<=N; k++)
	  {
	    ColumnVector tmp(nparams);
	    tmp=0;		  
	    for(int n=0; n<nparams; n++)
	      tmp(n+1)=ss_values_subsampled[n](k);
	    
	    SymmetricMatrix tmp2;
	    tmp2 << (tmp-mns)*(tmp-mns).t();	   
	    
	    proposal_cov+=tmp2;		
	  }	

// 	if(N-nparams<20 && debuglevel==4)
// 	  {
// 	    LOGOUT("number of jumps per proposal cov update too low");
// 	    LOGOUT(N);
// 	    LOGOUT(nparams);
// 	    LOGOUT(naccepted_cov/float(subsample_num));
// 	    LOGOUT(naccepted_cov);
// 	  }
	
	proposal_cov=proposal_cov/(N-nparams);   

	if(cov_mode==1)
	  {
	    // calc average cross corr
	    float av=0;
	    float av_cov=0;

	    for(int n=1; n<=nparams; n++)
	      {
		av+=proposal_cov(n,n);
		for(int n2=n+1; n2<=nparams; n2++)		
		  av_cov+=proposal_cov(n,n2);
	      }

	    av/=nparams;
	    av_cov/=(nparams*(nparams-1)/2.0);

	    for(int n=1; n<=nparams; n++)
	      {
		//proposal_cov(n,n)=av;
		for(int n2=n+1; n2<=nparams; n2++)
		  proposal_cov(n,n2)=av_cov;
	      }
	  }
	else if(cov_mode==2)
	  {
	    for(int n=1; n<=nparams; n++)
	      {
		for(int n2=n+1; n2<=nparams; n2++)		
		  proposal_cov(n,n2)=0;
	      }
	  }
	else if(cov_mode==3)
	  {
	    proposal_cov<<IdentityMatrix(nparams);
	  }

	double ssnormalise=1.0;//sqrt(prod(diag(cov)));
	//double ssnormalise=0.0;
	for(int k=1; k<=proposal_cov.Nrows(); k++)
	  {
	    ssnormalise*=pow(proposal_cov(k,k),1.0/nparams);
	    //ssnormalise+=sqrt(proposal_cov(k,k))/nparams;
	  }
	//ssnormalise=Sqr(ssnormalise);
	proposal_cov=(proposal_cov/ssnormalise);	     	      
	
	// take weighted average of this and last proposal cov
	proposal_cov=0.5*proposal_cov+0.5*proposal_cov_old;

	try{
	  normrand.setcovar(proposal_cov*Sqr(proposal_std));	
	}
	catch(Exception& cock)
	  {
	    LOGOUT("proposal cov no good");
	    LOGOUT(parameter_names[0]); 
	    LOGOUT(ssnormalise);
	    LOGOUT(N);
	    LOGOUT(proposal_std);
	    LOGOUT(proposal_cov);
	    LOGOUT(proposal_cov_old);
	    LOGOUT(naccepted_cov);
	    LOGOUT(naccepted);
	    LOGOUT(nrejected_cov);
	    LOGOUT(nrejected);
	    LOGOUT(calc_energy());
	    
	    for(int n=0; n<nparams; n++)
	      {
		OUT(parameter_names[n]);
		LOGOUT(param_max[n]);
		LOGOUT(param_min[n]);
		ColumnVector cols=vector2ColumnVector(samples[n]);
		
		write_ascii_matrix(cols, LogSingleton::getInstance().appendDir(parameter_names[n]+"_samples"));	

		cols=vector2ColumnVector(ss_values[n]);
		write_ascii_matrix(cols, LogSingleton::getInstance().appendDir(parameter_names[n]+"_ss"));	
		
	      }
	    exit(1);
	  }  

	if(debuglevel==4)// && parameter_names[0]=="ipl_art_hrf_2")
	  {   
	    //	    LOGOUT("Update proposal_cov");
	    LOGOUT(N);
	    LOGOUT(proposal_std);
	    LOGOUT(proposal_cov);
	  } 
      }
    else
      {
      	if(debuglevel==4)// && parameter_names[0]=="ipl_art_hrf_2") 
	  cout<<"There are NOT enough acceptances to update cov"<< endl;
      }
  
    naccepted_cov = 0;
    nrejected_cov = 0;
    for(int n=0; n<nparams; n++)
      ss_values[n].clear();
  }
  
  void Mcmc_Component::save()
  {
    Tracer_Plus trace("Mcmc_Component::save");

    if(save_out)
      for(int n=0; n<nparams; n++)
	{
	  ColumnVector cols=vector2ColumnVector(samples[n]);
	  write_ascii_matrix(mean(cols), LogSingleton::getInstance().appendDir(parameter_names[n]+"_mean"));
	  write_ascii_matrix(var(cols), LogSingleton::getInstance().appendDir(parameter_names[n]+"_var"));
	  
	  if(output_samples)
	    {	      
	      write_ascii_matrix(cols, LogSingleton::getInstance().appendDir(parameter_names[n]+"_samples"));	
	    }
	}
  }  

  double Mcmc_Component::vb_jump(double old_lik_energy)
  {
    Tracer_Plus trace("Mcmc_Component::vb_jump");

//     if(debuglevel==5)
//       LOGOUT("********************* start")

    // save old values
    vector<float> old_values=values;
    store_old();

    double old_prior_energy = calc_prior_energy(vector2ColumnVector(values),prior_mean,prior_precision);
    double old_energy = old_lik_energy+old_prior_energy;

    // propose new parameters based on a nonlinear vb update
    //    OUT("Compute Jacobian");
    calc_jacobian(jacob);
    //   OUT(size(jacob));
    calc_forward_model(fm);
    //    OUT(size(fm));
    calc_noise_precision(noise_prec);    
    //    OUT(size(noise_prec));

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
    //    ColumnVector beta=jacob.t()*noise_prec*((get_vb_data()-fm)+jacob*vector2ColumnVector(old_values))+prior_precision*prior_mean;
    prec<<tmp*jacob+prior_precision;
    ColumnVector beta=tmp*((get_vb_data()-fm)+jacob*vector2ColumnVector(old_values))+prior_precision*prior_mean;

    values=ColumnVector2vector(prec.i()*beta);
    set_values();

    // test energy     
    bool accept=true;
    // calc new energy
    double new_lik_energy = calc_energy();         	
    double new_prior_energy = calc_prior_energy(vector2ColumnVector(values),prior_mean,prior_precision);
    double new_energy = new_lik_energy+new_prior_energy;

    bool inbounds=true;

    if(!in_bounds())
      {
	inbounds=false;
	accept=false;
      }

    if(debuglevel==6)
      {	
// 	OUT(prec);
// 	OUT(jacob.t()*noise_prec*((get_vb_data()-fm)+jacob*vector2ColumnVector(old_values)));
// 	OUT(prior_precision*prior_mean);
	OUT(beta);
	
	for(int n=0; n<nparams; n++)
	  LogSingleton::getInstance().str() << parameter_names[n] << " (old)=" << old_values[n] << endl;
	for(int n=0; n<nparams; n++)
	  LogSingleton::getInstance().str() << parameter_names[n] << "=" << values[n] << endl;
	
	LogSingleton::getInstance().str() << "old_energy()=" << old_energy << endl;
	LogSingleton::getInstance().str() << "new_energy()=" << new_energy << endl;
	LogSingleton::getInstance().str() << "accept=" << accept << endl;
	LOGOUT(in_bounds());
	LOGOUT(samples[0].size());
	
	LOGOUT("********************* end");
	
	// 	 write_ascii_matrix(get_vb_data(), LogSingleton::getInstance().appendDir("data"));	
	// 	 write_ascii_matrix(fm, LogSingleton::getInstance().appendDir("fm"));
	// 	 write_ascii_matrix(vector2ColumnVector(old_values), LogSingleton::getInstance().appendDir("old_values"));
	// 	 write_ascii_matrix(vector2ColumnVector(values), LogSingleton::getInstance().appendDir("values"));
	// 	 write_ascii_matrix(jacob, LogSingleton::getInstance().appendDir("jacob"));
	
      }
    
    if(accept && new_energy<old_energy)
      {
	accept=true;
      }
    else 
      {
	accept=false;
      }

    if(accept)
      {
	//	LOGOUT("ACCEPTED");
	old_lik_energy = new_lik_energy;
      }
    else // reject      
      {      
	//	LOGOUT("REJECTED");
	values=old_values;
	set_values();
	restore_old();
      }

    return old_lik_energy;
      
  }

  double Mcmc_Component::jump(double old_lik_energy)
  {
    Tracer_Plus trace("Mcmc_Component::jump");

//     if(debuglevel==5)
//       LOGOUT("********************* start")

    // save old values
    vector<float> old_values=values;
    store_old();
  
    double old_prior_energy = calc_prior_energy(vector2ColumnVector(values),prior_mean,prior_precision);
    double old_energy = old_lik_energy+old_prior_energy;

//      double old_energy_debug=calc_energy(prior_mean,prior_precision);

//     if(old_energy!=old_energy_debug)
//       {
// 	LOGOUT("shiTTTTTT");
// 	LOGOUT(old_energy);
// 	LOGOUT(old_energy_debug);
// 	LOGOUT(old_energy-old_energy_debug);
//       }

//      old_energy=old_energy_debug;

//     if(debuglevel==3 && parameter_names[0]=="sigmaa")
//       {
// 	LOGOUT("-----------");
// 	if(samples[0].size()==50)
// 	  set_debuglevel(5);

// 	LogSingleton::getInstance().str() << "Old energy=" << calc_energy() << endl;

// 	if(samples[0].size()==50)
// 	  set_debuglevel(3);

//       }

    if((debuglevel==3 && parameter_names[0]=="DLPFC_mvn_mean_standard_space_1") || debuglevel==5)
       {
	LOGOUT("********************* start")
	for(int n=0; n<nparams; n++)
	  LogSingleton::getInstance().str() << parameter_names[n] << "_old=" << values[n] << endl;

	for(int n=0; n<nparams; n++)
	  {
	    LogSingleton::getInstance().str() << parameter_names[n] << "_max=" << param_max[n] << endl;	
	    LogSingleton::getInstance().str() << parameter_names[n] << "_min=" << param_min[n] << endl;	
	  }

	LOGOUT(in_bounds());
      }

    // propose new values      
    RowVector mn(nparams);
    for(int n=0; n<nparams; n++)
      mn(n+1)=values[n];

    RowVector newmn(nparams);  

    newmn=normrand.next(mn);

    for(int n=0; n<nparams; n++)
      values[n]=newmn(n+1);

    if((debuglevel==3 && parameter_names[0]=="DLPFC_mvn_mean_standard_space_1") || debuglevel==5)
      {
	for(int n=0; n<nparams; n++)
	  {
	    LogSingleton::getInstance().str() << parameter_names[n] << "_new=" << values[n] << endl;
	  }
	
      }

    set_values();
      
    bool accept=true;
    bool inbounds=true;
    double new_energy = 1e32;
    double new_lik_energy = 1e32;

    if(!in_bounds())
      {
	inbounds=false;
	accept=false;

	if((debuglevel==3 && parameter_names[0]=="DLPFC_mvn_mean_standard_space_1") || debuglevel==5)
 	  {
	    LOGOUT("Parameter out of range in Mcmc_Component::jump");
	  }
      }
    else
      {
	// calc new energy
	new_lik_energy = calc_energy();         	
	double new_prior_energy = calc_prior_energy(vector2ColumnVector(values),prior_mean,prior_precision);
	new_energy = new_lik_energy+new_prior_energy; 
      }
      
    if(accept)
      {
	// test acceptance
	double tmp = unifrnd().AsScalar();   
	accept = ((old_energy - new_energy) > std::log(tmp));
      }
      
    if((debuglevel==3 && parameter_names[0]=="DLPFC_mvn_mean_standard_space_1") || debuglevel==5)
       {
	for(int n=0; n<nparams; n++)
	  LogSingleton::getInstance().str() << parameter_names[n] << "=" << values[n] << endl;
	
	LogSingleton::getInstance().str() << "old_energy()=" << old_energy << endl;
	LogSingleton::getInstance().str() << "new_energy()=" << new_energy << endl;
	LogSingleton::getInstance().str() << "new_energy()==old_energy()=" << (new_energy==old_energy) << endl;
	LogSingleton::getInstance().str() << "new_energy()-old_energy()=" << (new_energy-old_energy) << endl;
	LogSingleton::getInstance().str() << "accept=" << accept << endl;
	LOGOUT(in_bounds());
	LOGOUT(samples[0].size());

	LOGOUT("********************* end")
      }
        
    if(accept)
      {		  	 
	nacceptedtotal++;
	naccepted++;
	naccepted_cov++;

	bool nochange=true;
	for(int n=0; n<nparams; n++)
	  {
	    double del=values[n]-old_values[n];
	    if(del!=0) nochange=false;

	    ss_values[n].push_back(del);
	  }

	if(nochange)
	  {
	    for(int n=0; n<nparams; n++)
	      {
		LOGOUT(parameter_names[n]);
		double del=values[n]-old_values[n];
		LOGOUT(del);
	      }

	    LOGOUT("no change, but the jump is being accepted");
	    
	    LOGOUT(proposal_std);
	    LOGOUT(proposal_cov);

	    LogSingleton::getInstance().str() << "old_energy()=" << old_energy << endl;
	    LogSingleton::getInstance().str() << "new_energy()=" << new_energy << endl;
	    LogSingleton::getInstance().str() << "new_energy()-old_energy()=" << (new_energy-old_energy) << endl;

// 	    LogSingleton::getInstance().str() << "new_energy()-old_energy_debug()=" << (new_energy-old_energy_debug) << endl;
// 	    LogSingleton::getInstance().str() << "old_energy()-old_energy_debug()=" << (old_energy-old_energy_debug) << endl;
// 	    	    LogSingleton::getInstance().str() << "old_energy_debug()=" << old_energy_debug << endl;

	    LOGOUT(in_bounds());
	  }      

	old_lik_energy = new_lik_energy;
      }
    else // reject      
      {
	//	if((inbounds && !isnan(new_energy)))// || (old_energy == new_energy))
	  {
	    nrejectedtotal++;
	    nrejected++;
	    nrejected_cov++;
	  }

	values=old_values;
	set_values();
	restore_old();

      }

    return old_lik_energy;
  }

  void Mcmc_Component::set_values_to_sample_mean()
  {
    Tracer_Plus trace("Mcmc_Component::set_values_to_sample_mean");

    for(int n=0; n<nparams; n++)
      {
	values[n]=mean(vector2ColumnVector(samples[n])).AsScalar();
      }

  }

  void Mcmc_Component::set_values_to_sample(int index)
  {
    Tracer_Plus trace("Mcmc_Component::set_values_to_sample_mean");

    for(int n=0; n<nparams; n++)
      {
	values[n]=samples[n][index-1];
      }

  }

  void copy_proposal_distributions(const vector<Mcmc_Component*>& from, vector<Mcmc_Component*> to)
  {
    Tracer_Plus trace("copy_proposal_distributions");
    
    // copies matching components (matched by name) proposal distributions
    for(unsigned int i=0; i<from.size(); i++)      
      {
	bool match=false;
	for(unsigned int j=0; j<to.size() && !match; j++)
	  {
	    if(from[i]->name==to[j]->name && from[i]->nparams==to[j]->nparams)
	      {
		match=true;
		
		to[j]->proposal_std=from[i]->proposal_std;
		to[j]->proposal_cov=from[i]->proposal_cov;
	      }
	  }
      }
      
  }
}
