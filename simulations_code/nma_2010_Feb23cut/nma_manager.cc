 /*  nma_manager.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT  */

#include "nma_manager.h"
#include "utils/log.h"
#include "libprob/libprob.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include "subject.h"
#include "model.h"
//#include "mcmc_hmc.h"
#include <algorithm>
#include <iostream>

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace std;
using namespace MISCPLOT;
using namespace MISCPIC;

namespace Nma {

  ReturnMatrix devar(const Matrix& mat, const int dim)
  { 
    Matrix res;
    if (dim == 1) {res=mat;}
    else {res=mat.t();}
    
    Matrix sd;
    sd = stdev(res);
    
    for (int ctr = 1; ctr <= res.Nrows(); ctr++) {
      for (int ctr2 =1; ctr2 <= res.Ncols(); ctr2++) {
	res(ctr,ctr2)/=sd(1,ctr2);
      }
    }
    if (dim != 1) {res=res.t();}
    res.Release();
    return res;
  }

  const double bvn_energy(float& x1, float& x2)   
  {
      Tracer_Plus trace("bvn_energy");

      float rho=0.998;
      float sig1=0.1;
      float sig2=1;

      float z=Sqr(x1/sig1) - 2*rho*x1*x2/(sig1*sig2) + Sqr(x2/sig2);

      return z/(2*(1-Sqr(rho)));      
  }

  void Nma_manager::setup()
  {
    Tracer_Plus trace("Nma_manager::setup");

    // setup model
    num_subjects=opts.subject_names.value().size();

    bool random_initialise=opts.rand_init.value();

    LOGOUT(random_initialise);
    bool sts=true; //single_timeseries

    double res=opts.tr.value()/opts.resfactor.value();
    model_sts=new Model(sts, opts.data_mode.value(), opts.haemodynamic_model.value(), opts.data_directory.value(), opts.model_directory.value(), opts.node_names.value(), opts.stimuli_names.value(), opts.stim_amp_mod.value(), opts.stim_ard.value(), opts.stimuli_are_single_col_format.value(), res, opts.tr.value(), opts.decode.value());

    // setup subjects
    for(int s=0; s<num_subjects; s++)
      {
	subjects_sts.push_back(new Subject(opts.subject_names.value()[s], opts.data_directory.value(), *model_sts, sts, opts.data_mode.value(), opts.haemodynamic_model.value(),random_initialise,opts.debuglevel.value(), opts.decode.value(), opts.phi_every_voxel.value()));  
      }

    if(opts.data_mode.value()==1 || opts.data_mode.value()==3)
      {
	// model for doing full inference on roi data
	sts=false;
	model_roi=new Model(sts, opts.data_mode.value(), opts.haemodynamic_model.value(), opts.data_directory.value(), opts.model_directory.value(), opts.node_names.value(), opts.stimuli_names.value(), opts.stim_amp_mod.value(), opts.stim_ard.value(), opts.stimuli_are_single_col_format.value(), res, opts.tr.value(), opts.decode.value());

	// setup subjects
	for(int s=0; s<num_subjects; s++)
	  {
	    subjects_roi.push_back(new Subject(opts.subject_names.value()[s], opts.data_directory.value(), *model_roi, sts, opts.data_mode.value(), opts.haemodynamic_model.value(),random_initialise,opts.debuglevel.value(), opts.decode.value(), opts.phi_every_voxel.value()));
	  }

	// model for doing no roi inference on initial runs using basis set
	sts=false;

	model_basis=new Model(sts, opts.data_mode.value(), opts.haemodynamic_model.value(), opts.data_directory.value(), opts.model_directory.value(), opts.node_names.value(), opts.stimuli_names.value(), opts.stim_amp_mod.value(), opts.stim_ard.value(), opts.stimuli_are_single_col_format.value(), res, opts.tr.value(), true);

	// setup subjects
	for(int s=0; s<num_subjects; s++)
	  {
	    subjects_basis.push_back(new Subject(opts.subject_names.value()[s], opts.data_directory.value(), *model_basis, sts, opts.data_mode.value(), opts.haemodynamic_model.value(),random_initialise,opts.debuglevel.value(), true, opts.phi_every_voxel.value()));
	  }

	model_basis_best=new Model(sts, opts.data_mode.value(), opts.haemodynamic_model.value(), opts.data_directory.value(), opts.model_directory.value(), opts.node_names.value(), opts.stimuli_names.value(), opts.stim_amp_mod.value(), opts.stim_ard.value(), opts.stimuli_are_single_col_format.value(), res, opts.tr.value(), true);

	// setup subjects
	for(int s=0; s<num_subjects; s++)
	  {
	    subjects_basis_best.push_back(new Subject(opts.subject_names.value()[s], opts.data_directory.value(), *model_basis_best, sts, opts.data_mode.value(), opts.haemodynamic_model.value(),random_initialise,opts.debuglevel.value(), true, opts.phi_every_voxel.value()));
	  }
      }

  }

  void Nma_manager::initialising_runs(int s, string dirname, const vector<Mcmc_Component*>& components_sts, vector<Mcmc_Component*>& components_basis_best)
  {
    Tracer_Plus trace("Nma_manager::initialising_runs");	
	    
    // working subjects
 
    int nnodes=subjects_basis[s]->get_subject_model().get_subject_nodes().size();    
    
    // loop through nodes 
    for(int n=1; n<=nnodes; n++)
      {
	// set up list of sub ROIs within this node
	Subject_Node& subject_node = (*subjects_basis[s]->get_subject_model().get_subject_nodes()[n-1]);
	
	subject_node.setup_mvn_basis_set();

	// loop through basis set running inference
	int best=1;
	double best_loglik=-1e64;

	for(int b=1; b<=subject_node.get_num_basis(); b++)
	  {
	    LogSingleton::getInstance().setthenmakeDir(dirname+string("/logdir_initial_run_")+num2str(n)+"_"+num2str(b),"logfile",true,true);
	    
	    LOGOUT(LogSingleton::getInstance().getDir());

	    subjects_basis[s]->get_subject_model().copy_neural_connectivity_values(subjects_sts[s]->get_subject_model());
	    subjects_basis[s]->get_subject_model().copy_hrf_values(subjects_sts[s]->get_subject_model());

	    subject_node.set_mvn_basis(b);
	    subjects_basis[s]->get_subject_model().evaluate_decoded_node_data();
	    
	    // working components
	    vector<Mcmc_Component*> components_basis;	 	
	    setup_mcmc_components(*(subjects_basis[s]),2,false,components_basis);
	    copy_proposal_distributions(components_basis_best,components_basis);	    

	    double loglik=run_mcmc(*(subjects_basis[s]),components_basis,120,0,1,-1,-1,false); 

	    OUT(b);
	    OUT(loglik);
	    if(loglik>best_loglik)
	      {
		OUT(best);
		OUT(best_loglik);

		best_loglik=loglik;
		best=b;

		copy_proposal_distributions(components_basis,components_basis_best);	    
	    
		subjects_basis_best[s]->get_subject_model().copy_neural_connectivity_values(subjects_basis[s]->get_subject_model());
		subjects_basis_best[s]->get_subject_model().copy_hrf_values(subjects_basis[s]->get_subject_model());
		subjects_basis_best[s]->get_subject_model().copy_roi_values(subjects_basis[s]->get_subject_model());	    
	      }

	    for(unsigned int r=0; r<components_basis.size(); r++) delete components_basis[r]; components_basis.clear();
	  }
	
	// copy best param values into working subject and components
	subjects_basis[s]->get_subject_model().copy_neural_connectivity_values(subjects_basis_best[s]->get_subject_model());
	subjects_basis[s]->get_subject_model().copy_hrf_values(subjects_basis_best[s]->get_subject_model());
	subjects_basis[s]->get_subject_model().copy_roi_values(subjects_basis_best[s]->get_subject_model());
		
      }

  }

  void Nma_manager::run()
  {
    Tracer_Plus trace("Nma_manager::run");

    for(int s=0; s<num_subjects; s++)
      {
	string dirname=LogSingleton::getInstance().getDir();

	if(opts.data_mode.value()==1)
	  {	   
	    ///////////////////////////      
	    // do non-roi inference with haemodynamics on roi averaged single timeseries data
	    LogSingleton::getInstance().setthenmakeDir(dirname+string("/logdir_single_timeseries"),"logfile",true,true);

 	    LOGOUT(LogSingleton::getInstance().getDir());

	    vector<Mcmc_Component*> components_sts;
	    setup_mcmc_components(*(subjects_sts[s]),1,true,components_sts);
	    run_mcmc(*(subjects_sts[s]),components_sts,100,0,1,30,30);

	    /////////////////////////// 
	    // iterate through spatial MVN basis functions to initialise ROI locations

	    vector<Mcmc_Component*> components_basis_best;
	    setup_mcmc_components(*(subjects_basis_best[s]),2,false,components_basis_best);
	    initialising_runs(s, dirname, components_sts, components_basis_best);

	    ///////////////////////////
	    // do non-roi inference only on roi data
	    LogSingleton::getInstance().setthenmakeDir(dirname+string("/logdir_roi_no_roi_infer"),"logfile",true,true);

 	    LOGOUT(LogSingleton::getInstance().getDir());	    

	    // transfer results from initialising run
	    subjects_roi[s]->get_subject_model().copy_neural_connectivity_values(subjects_basis_best[s]->get_subject_model());
	    subjects_roi[s]->get_subject_model().copy_hrf_values(subjects_basis_best[s]->get_subject_model());
	    subjects_roi[s]->get_subject_model().copy_roi_values(subjects_basis_best[s]->get_subject_model());	    
	    
	    vector<Mcmc_Component*> components_roi_no_roi_infer;
	    setup_mcmc_components(*(subjects_roi[s]),2,true,components_roi_no_roi_infer);
	    copy_proposal_distributions(components_basis_best,components_roi_no_roi_infer);	    

// 	    for(unsigned int i=0; i< components_roi_no_roi_infer.size(); i++)
// 	      components_roi_no_roi_infer[i]->set_debuglevel(5);

// 	    subjects_roi[s]->get_subject_model().set_debuglevel(5);

	    run_mcmc(*(subjects_roi[s]),components_roi_no_roi_infer,120,0,1,30,30);

	    //	    subjects_roi[s]->get_subject_model().output_roi_values(); 

	    ///////////////////////////
	    // do full inference on roi data
	    LogSingleton::getInstance().setDir(dirname,"logfile",true,true);
	    
	    LOGOUT(LogSingleton::getInstance().getDir());
	    
	    vector<Mcmc_Component*> components_roi;
	    setup_mcmc_components(*(subjects_roi[s]),3,true,components_roi);	    
	    copy_proposal_distributions(components_roi_no_roi_infer, components_roi);
	    //	    subjects_roi[s]->get_subject_model().output_roi_values(); 

// 	    for(unsigned int i=0; i< components_roi.size(); i++)
//  	      components_roi[i]->set_debuglevel(5);
//  	    subjects_roi[s]->get_subject_model().set_debuglevel(5);

	    run_mcmc(*(subjects_roi[s]),components_roi,opts.nsamps.value(),opts.burnin.value(),opts.sampleevery.value(),-1,30);

	    // delete components
	    for(unsigned int r=0; r<components_sts.size(); r++) delete components_sts[r]; components_sts.clear();
	    for(unsigned int r=0; r<components_basis_best.size(); r++) delete components_basis_best[r]; components_basis_best.clear();
	    for(unsigned int r=0; r<components_roi_no_roi_infer.size(); r++) delete components_roi_no_roi_infer[r]; components_roi_no_roi_infer.clear();
	    for(unsigned int r=0; r<components_roi.size(); r++) delete components_roi[r]; components_roi.clear();
	  }
	else if(opts.data_mode.value()==2 || opts.data_mode.value()==0)
	  {
	    ///////////////////////////
	    // do non-haemo inference 
	    LogSingleton::getInstance().setthenmakeDir(dirname+string("/logdir_no_haemodynamics"),"logfile",true,true);

 	    LOGOUT(LogSingleton::getInstance().getDir());
	    vector<Mcmc_Component*> components_no_haemodynamics;
	    setup_mcmc_components(*(subjects_sts[s]),1,false,components_no_haemodynamics);

	    //run_mcmc(Subject& subject, const vector<Mcmc_Component*>& components, int nsamps, int burnin, int se, int sample_vb_jump_every, int burnin_vb_jump_every, bool write_report)
	    run_mcmc(*(subjects_sts[s]),components_no_haemodynamics,200,0,1,-1,-1);

	    ///////////////////////////
	    // do full inference 
	    LogSingleton::getInstance().setDir(dirname,"logfile",true,true);
	    
	    LOGOUT(LogSingleton::getInstance().getDir());
	   
	    vector<Mcmc_Component*> components;
	    setup_mcmc_components(*(subjects_sts[s]),1,true,components);
  	    copy_proposal_distributions(components_no_haemodynamics, components);

// 	    for(unsigned int i=0; i< components.size(); i++)
//  	      components[i]->set_debuglevel(5);

//  	    subjects_sts[s]->get_subject_model().set_debuglevel(5);

	    run_mcmc(*(subjects_sts[s]),components,opts.nsamps.value(),opts.burnin.value(),opts.sampleevery.value(),-1,-1);

	    // delete components
  	    for(unsigned int r=0; r<components_no_haemodynamics.size(); r++) delete components_no_haemodynamics[r]; components_no_haemodynamics.clear();

	    for(unsigned int r=0; r<components.size(); r++) delete components[r]; components.clear();
	  }
	else if(opts.data_mode.value()==3)
	  {

	    ///////////////////////////
	    // do inference with fixed MVN ROI covariance on roi data
	    LogSingleton::getInstance().setthenmakeDir(dirname+string("/logdir_no_cov_infer"),"logfile",true,true);

 	    LOGOUT(LogSingleton::getInstance().getDir());

	    vector<Mcmc_Component*> components_no_cov_infer;
	    setup_mcmc_components(*(subjects_roi[s]),4,true,components_no_cov_infer);
	    run_mcmc(*(subjects_roi[s]),components_no_cov_infer,120,0,1,30,30);

// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[0]->get_mvn_sqrt_cov_value());
// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[1]->get_mvn_sqrt_cov_value());
// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[2]->get_mvn_sqrt_cov_value());

// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[0]->get_mvn_mean_func_space_value().t());
// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[1]->get_mvn_mean_func_space_value().t());
// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[2]->get_mvn_mean_func_space_value().t());

	    ///////////////////////////
	    // do full inference on roi data
	    LogSingleton::getInstance().setDir(dirname,"logfile",true,true);

	    vector<Mcmc_Component*> components;
	    setup_mcmc_components(*(subjects_roi[s]),3,true,components);

 // 	    subjects_roi[s]->get_subject_model().copy_neural_connectivity_values(subjects_roi[s]->get_subject_model());
//  	    subjects_roi[s]->get_subject_model().copy_hrf_values(subjects_roi[s]->get_subject_model());
//  	    subjects_roi[s]->get_subject_model().copy_roi_values(subjects_roi[s]->get_subject_model()); 

// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[0]->get_mvn_sqrt_cov_value());
// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[1]->get_mvn_sqrt_cov_value());
// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[2]->get_mvn_sqrt_cov_value());

// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[0]->get_mvn_mean_func_space_value().t());
// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[1]->get_mvn_mean_func_space_value().t());
// 	    OUT(subjects_roi[s]->get_subject_model().get_subject_nodes()[2]->get_mvn_mean_func_space_value().t());

	    run_mcmc(*(subjects_roi[s]),components,opts.nsamps.value(),opts.burnin.value(),opts.sampleevery.value(),-1,30);

	    ////////////////////
	    // delete components
	    for(unsigned int r=0; r<components_no_cov_infer.size(); r++) delete components_no_cov_infer[r]; components_no_cov_infer.clear();
	    for(unsigned int r=0; r<components.size(); r++) delete components[r]; components.clear();
	  }
      }
  }

  void Nma_manager::setup_mcmc_components(Subject& subject, int flag, bool do_haemodynamics, vector<Mcmc_Component*>& components)
  {
    Tracer_Plus trace("Nma_manager::setup_mcmc_components");

    for(unsigned int r=0; r<components.size(); r++) delete components[r]; components.clear();	

    // Add MCMC components to MCMC object    
    vector<string> names;       
    vector<float> values;
    vector<float> mins;
    vector<float> maxs;    
    vector<float> gradient_step_sizes;
    ColumnVector prior_mean;
    ColumnVector prior_isard;
    SymmetricMatrix prior_precision;

    // a component 
    names = subject.get_subject_model().get_names_a();   
    if(names.size()>0 && names.size()<4)
      {
	values = subject.get_subject_model().get_value_a_vec();
	mins.clear(); mins.resize(values.size(),-10), maxs.clear(); maxs.resize(values.size(),10);
	gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
	prior_mean.ReSize(values.size()); prior_mean=0;
	ColumnVector prior_precision_col(values.size());
	//	prior_precision_col=1e-6; prior_precision<<diag(prior_precision_col);
	// dcm eqn A9
	int N=subject.get_subject_model().get_subject_nodes().size(); 
	//	prior_precision_col=1.0/(N*(N-1)/chi2inv(1-1e-5,N*(N-1)));
	prior_precision_col=1.0/(N*(N-1)/chdtri(N*(N-1),1e-10));
	//	OUT(1.0/(N*(N-1)/chdtri(N*(N-1),1e-5)));
	prior_precision<<diag(prior_precision_col);

// 	LOGOUT(names);
// 	LOGOUT(values);
// 	LOGOUT(mins);
// 	LOGOUT(maxs);

	components.push_back(new A_Mcmc_Component(subject, "A", names, values, mins, maxs,gradient_step_sizes,opts.debuglevel.value(),opts.output_samples.value(),prior_mean,prior_precision));
      }
    else if (names.size()>0)
      {
	for(unsigned int n=1; n<=subject.get_subject_model().get_model().get_nodes().size(); n++)
	  {
	    names = subject.get_subject_model().get_names_a(n);
	    if(names.size()>0)
	      {
		values = subject.get_subject_model().get_value_a_vec(n);

		mins.clear(); mins.resize(values.size(),-10), maxs.clear(); maxs.resize(values.size(),10);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		//	prior_precision_col=1e-6; prior_precision<<diag(prior_precision_col);
		// dcm eqn A9
		int N=subject.get_subject_model().get_subject_nodes().size(); 
		//	prior_precision_col=1.0/(N*(N-1)/chi2inv(1-1e-5,N*(N-1)));
// 		OUT(N);
		prior_precision_col=1.0/(N*(N-1)/chdtri(N*(N-1),1e-10));
// 		OUT(1.0/(N*(N-1)/chdtri(N*(N-1),1e-5)));
// 		OUT(1.0/(N*(N-1)/chdtri(N*(N-1),1-1e-5)));
		prior_precision<<diag(prior_precision_col);

		components.push_back(new A_Node_Mcmc_Component(subject, "A_"+num2str(n), names, values, mins, maxs,gradient_step_sizes, n, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));
	      }
	  }
      }

    // b component
    names = subject.get_subject_model().get_names_b();	
    if(names.size()>0 && names.size()<4)
      {
	values = subject.get_subject_model().get_value_b_vec();
	mins.clear(); 
	mins.resize(values.size(),-10); 
	//mins.resize(values.size(),0); 
	maxs.clear(); maxs.resize(values.size(),10);
	gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
	prior_mean.ReSize(values.size()); prior_mean=0;
	ColumnVector prior_precision_col(values.size());
	prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
	int N=subject.get_subject_model().get_subject_nodes().size(); 
	prior_precision_col=(1.0/(N*(N-1)/chdtri(N*(N-1),1e-5)))/4.0;
	prior_precision<<diag(prior_precision_col);
// 	LOGOUT(names);
// 	LOGOUT(mins);
// 	LOGOUT(maxs);
	components.push_back(new B_Mcmc_Component(subject, "B", names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));
      }
    else if (names.size()>0)
      {
	for(unsigned int n=1; n<=subject.get_subject_model().get_model().get_nodes().size(); n++)
	  {
	    names = subject.get_subject_model().get_names_b(n);
	    if(names.size()>0)
	      {
		values = subject.get_subject_model().get_value_b_vec(n);

		mins.clear(); mins.resize(values.size(),-10);
		maxs.clear(); maxs.resize(values.size(),10);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		int N=subject.get_subject_model().get_subject_nodes().size(); 
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		prior_precision_col=(1.0/(N*(N-1)/chdtri(N*(N-1),1e-5)))/4.0;
		prior_precision<<diag(prior_precision_col);

		components.push_back(new B_Node_Mcmc_Component(subject, "B_"+num2str(n), names, values, mins, maxs,gradient_step_sizes, n, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));
	      }
	  }
      }

    // b amp mod components
    for(unsigned int n=1; n<=subject.get_subject_model().get_model().get_nodes().size(); n++)
      for(unsigned int i=1; i<=subject.get_subject_model().get_model().get_marker_b()[n-1].size(); i++)
	{
	  //	      LOGOUT(i);LOGOUT(n);LOGOUT(subject.get_subject_model().get_model().is_b_amp_mod(n,i));
	  if(subject.get_subject_model().get_model().is_b_amp_mod(n,i))
	    {
	      names = subject.get_subject_model().get_names_b_amp_mod(n,i);	
	      values = subject.get_subject_model().get_value_b_amp_mod_vec(n,i);
	      vector<float> mins(values.size(),-10.0); vector<float> maxs(values.size(),5.0);
	      gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);

	      components.push_back(new B_Amp_Mod_Mcmc_Component(subject, string("B_amp_mod_")+num2str(n)+string("_")+num2str(i), names, values, mins, maxs,gradient_step_sizes, n, i, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));
	    }
	}

    // c component 
    names = subject.get_subject_model().get_names_c();
    float max_c;
    if(subject.get_subject_model().get_haemodynamic_model()=="balloon" || subject.get_subject_model().get_haemodynamic_model()=="balloon_epsilon")      
      max_c=1;
    else
      max_c=10;

    if(names.size()>0 && names.size()<4)
      {
	values = subject.get_subject_model().get_value_c_vec();
	mins.clear(); 
	mins.resize(values.size(),-max_c);
	//mins.resize(values.size(),0);
	maxs.clear(); maxs.resize(values.size(),max_c);
	gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
	prior_mean.ReSize(values.size()); prior_mean=0;
	ColumnVector prior_precision_col(values.size());
	prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
	prior_isard=prior_mean;prior_isard=0;

// 	LOGOUT(names);
// 	LOGOUT(mins);
// 	LOGOUT(maxs);

	components.push_back(new C_Mcmc_Component(subject, "C", names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision,prior_isard));
      }
    else if (names.size()>0)
      {
	for(unsigned int n=1; n<=subject.get_subject_model().get_model().get_nodes().size(); n++)
	  {
	    names = subject.get_subject_model().get_names_c(n);
	    if(names.size()>0)
	      {
		values = subject.get_subject_model().get_value_c_vec(n);
		mins.clear(); mins.resize(values.size(),-max_c);
		maxs.clear(); maxs.resize(values.size(),max_c);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		prior_isard=prior_mean;prior_isard=0;

		components.push_back(new C_Node_Mcmc_Component(subject, "C_"+num2str(n), names, values, mins, maxs,gradient_step_sizes, n, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision,prior_isard));
	      }
	  }
      }

    // c amp mod components
    for(unsigned int n=1; n<=subject.get_subject_model().get_model().get_nodes().size(); n++)
      for(unsigned int i=1; i<=subject.get_subject_model().get_model().get_marker_c()[n-1].size(); i++)
	{
	  // 	      LOGOUT(i);LOGOUT(n);LOGOUT(subject.get_subject_model().get_model().is_c_amp_mod(n,i));
	  if(subject.get_subject_model().get_model().is_c_amp_mod(n,i))
	    {
	      names = subject.get_subject_model().get_names_c_amp_mod(n,i);	
	      values = subject.get_subject_model().get_value_c_amp_mod_vec(n,i);
	      vector<float> mins(values.size(),-max_c); vector<float> maxs(values.size(),max_c);
	      gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
	
	      components.push_back(new C_Amp_Mod_Mcmc_Component(subject, string("C_amp_mod_")+num2str(n)+string("_")+num2str(i), names, values, mins, maxs,gradient_step_sizes, n, i, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));
	    }  
	}

    // d component
    names = subject.get_subject_model().get_names_d();

    if(names.size()>0 && names.size()<4)
      {
	values = subject.get_subject_model().get_value_d_vec();
	mins.clear();
	mins.resize(values.size(),-10);
	//mins.resize(values.size(),0); 
	maxs.clear(); maxs.resize(values.size(),10);
	gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
	prior_mean.ReSize(values.size()); prior_mean=0;
	ColumnVector prior_precision_col(values.size());
	prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
	// dcm eqn A9
	int N=subject.get_subject_model().get_subject_nodes().size(); 
	//	prior_precision_col=1.0/(N*(N-1)/chi2inv(1-1e-5,N*(N-1)));
	// 		OUT(N);
	prior_precision_col=(1.0/(N*(N-1)/chdtri(N*(N-1),1e-5)))/2.0;
	prior_precision<<diag(prior_precision_col);


// 	LOGOUT(names);
// 	LOGOUT(mins);
// 	LOGOUT(maxs);

	components.push_back(new D_Mcmc_Component(subject, "D", names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));
      }
    else if (names.size()>0)
      {
	for(unsigned int n=1; n<=subject.get_subject_model().get_model().get_nodes().size(); n++)
	  {
	    names = subject.get_subject_model().get_names_d(n);
	    if(names.size()>0)
	      {
		values = subject.get_subject_model().get_value_d_vec(n);
		mins.clear(); mins.resize(values.size(),-10);
		maxs.clear(); maxs.resize(values.size(),10);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		// dcm eqn A9
		int N=subject.get_subject_model().get_subject_nodes().size(); 
		//	prior_precision_col=1.0/(N*(N-1)/chi2inv(1-1e-5,N*(N-1)));
		// 		OUT(N);
		prior_precision_col=(1.0/(N*(N-1)/chdtri(N*(N-1),1e-5)))/4.0;
		prior_precision<<diag(prior_precision_col);

		components.push_back(new D_Node_Mcmc_Component(subject, "D_"+num2str(n), names, values, mins, maxs,gradient_step_sizes, n, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));
	      }
	  }
      }
    
    //    if(flag==3) 
    {
      // logsigma_a
      names = subject.get_subject_model().get_names_logsigmaa();
      values = subject.get_subject_model().get_value_logsigmaa_vec();
      OUT(values);
      float min_logsigmaa=log(1e-10); float max_logsigmaa=log(1e10);
      mins.clear(); mins.resize(1); mins[0]=min_logsigmaa;
      maxs.clear(); maxs.resize(1); maxs[0]=max_logsigmaa;
      gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
      prior_mean.ReSize(values.size()); prior_mean=subject.get_subject_model().get_prior_mean_logsigmaa();
      ColumnVector prior_precision_col(values.size());
      prior_precision_col=subject.get_subject_model().get_prior_precision_logsigmaa(); prior_precision<<diag(prior_precision_col);

      components.push_back(new Log_Sigmaa_Mcmc_Component(subject, "logsigmaa", names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));
    }
    
    if(subject.is_phi_every_voxel())
      {
	if((flag==3 || flag==2) && !subject.is_decode()) 
	  {
	    const vector<ColumnVector>&  coords=subject.get_subject_model().get_voxel_coordinates();		
	    for(unsigned int n=0; n<coords.size(); n++) // indexes voxels
	      {
		////////////////////////////////
		// do voxel variance
		names = subject.get_subject_model().get_names_log_phi_every_voxel(n+1);
		values = subject.get_subject_model().get_value_log_phi_every_voxel_vec(n+1);
		mins.clear(); mins.resize(1,-30); 
		maxs.clear(); maxs.resize(1,30);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);

		//  LOGOUT(names);
		// 	    LOGOUT(values);
	    
		if(names.size()>0)
		  //		if(names.size()>0 && n<3)
		  components.push_back(new Log_Phi_Every_Voxel_Mcmc_Component(subject, string("log_phi_every_voxel"), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision,n+1));	
	      }
	  }
      }
    //    else


      {
	if((flag==3 || flag==2) && !subject.is_decode()) 
	  {
	    ////////////////////////////////
	    // do voxel variance
	    names = subject.get_subject_model().get_names_log_phi_voxel();
	    values = subject.get_subject_model().get_value_log_phi_voxel_vec();
	    mins.clear(); mins.resize(1,-30); 
	    maxs.clear(); maxs.resize(1,30);
	    gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
	    prior_mean.ReSize(values.size()); prior_mean=0;
	    ColumnVector prior_precision_col(values.size());
	    prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);

	    LOGOUT(names);
	    LOGOUT(values);

	    if(names.size()>0)
	      //		if(names.size()>0 && n<3)
	      components.push_back(new Log_Phi_Voxel_Mcmc_Component(subject, string("log_phi_voxel"), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),prior_mean,prior_precision));	
	  }
      }

    //    for(unsigned int n=1; n<=1; n++)
    for(unsigned int n=1; n<=subject.get_subject_model().get_subject_nodes().size(); n++)
      {
	if((flag==3 || flag==2) && !subject.is_decode()) 
	  {
	    ////////////////////////////////
	    // do node variance
	    names = subject.get_subject_model().get_names_phi_node(n);
	    values = subject.get_subject_model().get_value_phi_node_vec(n);
	    mins.clear(); mins.resize(1,-30); 
	    maxs.clear(); maxs.resize(1,30);
	    gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
	    prior_mean.ReSize(values.size()); prior_mean=0;
	    ColumnVector prior_precision_col(values.size());
	    prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);

	    LOGOUT(names);
	    LOGOUT(values);
// 	    LOGOUT(mins);
// 	    LOGOUT(maxs);
	    
	    if(names.size()>0)
	      //		if(names.size()>0 && n<3)
	      components.push_back(new Phi_Node_Mcmc_Component(subject, string("log_phi_node_")+num2str(n), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),n,prior_mean,prior_precision));	
	  }

	////////////////////////////////
	Subject_Node& subject_node=*(subject.get_subject_model().get_subject_nodes()[n-1]);

	if(do_haemodynamics) 
	  {
 	    if(subject.get_subject_model().get_haemodynamic_model()=="balloon" || subject.get_subject_model().get_haemodynamic_model()=="balloon_epsilon")      
 	      {
		names = subject_node.get_names_balloon_cbf();
		values = subject_node.get_value_balloon_cbf();
		mins = subject_node.get_mins_balloon_cbf();
		maxs = subject_node.get_maxs_balloon_cbf();
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		prior_precision=subject_node.get_prior_prec_balloon_cbf();
		prior_mean=subject_node.get_prior_mean_balloon_cbf();
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		for(unsigned int j=1; j<=gradient_step_sizes.size(); j++)      
		  gradient_step_sizes[j-1]=Min((1.0/std::sqrt(prior_precision(j,j))*1e-4),1e-4);

// 		LOGOUT(names);
// 		LOGOUT(mins);
// 		LOGOUT(maxs);

		if(names.size()>0)
		  //		if(names.size()>0 && n<3)
  		  components.push_back(new Balloon_Cbf_Mcmc_Component(subject, string("Balloon_Cbf_")+num2str(n), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),n,prior_mean,prior_precision));	

		names = subject_node.get_names_balloon();
		values = subject_node.get_value_balloon();
		mins = subject_node.get_mins_balloon();
		maxs = subject_node.get_maxs_balloon();				
		prior_mean.ReSize(values.size()); prior_mean=0;
		prior_precision_col.ReSize(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		prior_precision=subject_node.get_prior_prec_balloon();
		prior_mean=subject_node.get_prior_mean_balloon();
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		for(unsigned int j=1; j<=gradient_step_sizes.size(); j++)      
		  {
		    gradient_step_sizes[j-1]=Min((1.0/std::sqrt(prior_precision(j,j))*1e-4),1e-4);
		    LOGOUT(1.0/std::sqrt(prior_precision(j,j)));
		  }

	    LOGOUT(gradient_step_sizes);

// 		LOGOUT(names);
// 		LOGOUT(mins);
// 		LOGOUT(maxs);

     		if(names.size()>0)
   		  components.push_back(new Balloon_Bold_Mcmc_Component(subject, string("Balloon_Bold_")+num2str(n), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),n,prior_mean,prior_precision));	

		if(subject.get_subject_model().get_haemodynamic_model()=="balloon_epsilon")
		  {
		    names = subject_node.get_names_balloon2();
		    values = subject_node.get_value_balloon2();
		    mins = subject_node.get_mins_balloon2();
		    maxs = subject_node.get_maxs_balloon2();				
		    prior_mean.ReSize(values.size()); prior_mean=0;
		    ColumnVector prior_precision_col(values.size());
		    prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		    prior_precision=subject_node.get_prior_prec_balloon2();
		    prior_mean=subject_node.get_prior_mean_balloon2();
		    gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		    for(unsigned int j=1; j<=gradient_step_sizes.size(); j++)      
		      {
			gradient_step_sizes[j-1]=Min((1.0/std::sqrt(prior_precision(j,j))*1e-4),1e-4);	
			LOGOUT(1.0/std::sqrt(prior_precision(j,j)));
		      }

		    LOGOUT(gradient_step_sizes);

// 		    LOGOUT(names);
// 		    LOGOUT(mins);
// 		    LOGOUT(maxs);
		    
   		    if(names.size()>0)
   		      components.push_back(new Balloon_Bold2_Mcmc_Component(subject, string("Balloon_Bold2_")+num2str(n), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),n,prior_mean,prior_precision));
 		  }
	      }
	    else
	      {
		// hrf component for each node
		names = subject_node.get_names_hrf();
		values = subject_node.get_value_hrf();
		mins = subject_node.get_mins_hrf();
		maxs = subject_node.get_maxs_hrf();	
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		for(unsigned int j=1; j<=gradient_step_sizes.size(); j++)      
		  gradient_step_sizes[j-1]=Min((1.0/std::sqrt(prior_precision(j,j))*1e-4),1e-4);

		// 	float min_sigmaa=0.5; float max_sigmaa=1.5;
		// 	mins[values.size()-1]=min_sigmaa;
		// 	maxs[values.size()-1]=max_sigmaa;
		
		// 	LOGOUT(n);
		// 	LOGOUT(names);
		// 	LOGOUT(values);
		// 	LOGOUT(mins);
		// 	LOGOUT(maxs);	
		
 		if(names.size()>0)
		  {
		    components.push_back(new HRF_Mcmc_Component(subject, string("HRF_")+num2str(n), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),n,prior_mean,prior_precision));	    
		  }
	      }
	  }
	
	if(flag==3) 
	  {	    	    
	    // mvn mean component for each node
	    names = subject_node.get_names_mvn_mean_standard_space();
	    
	    if(names.size()>0)
	      {
		values = subject_node.get_value_mvn_mean_standard_space();
		mins = subject_node.get_mins_mvn_mean_standard_space();
		maxs = subject_node.get_maxs_mvn_mean_standard_space();
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		for(unsigned int j=1; j<=gradient_step_sizes.size(); j++)      
		  gradient_step_sizes[j-1]=Min((1.0/std::sqrt(prior_precision(j,j))*1e-4),1e-4);

		LOGOUT(n);
		LOGOUT(names);
		LOGOUT(values);
		LOGOUT(mins);
		LOGOUT(maxs);	
		
// 		if(n==2)
// 		  subject.get_subject_model().output_roi_values(); 

 		components.push_back(new MVN_Mean_Mcmc_Component(subject, string("MVN_Mean_")+num2str(n), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),n,prior_mean,prior_precision));
	
// 		if(n==2)
// 		  subject.get_subject_model().output_roi_values(); 

	      }
	  }
	
	if(flag==3)
	  {
	    // mvn cov component for each node
	    names = subject_node.get_names_mvn_sqrt_cov();
	    
	    if(names.size()>0)
	      {
		values = subject_node.get_value_mvn_sqrt_cov();
		mins = subject_node.get_mins_mvn_sqrt_cov();
		maxs = subject_node.get_maxs_mvn_sqrt_cov();
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		for(unsigned int j=1; j<=gradient_step_sizes.size(); j++)      
		  gradient_step_sizes[j-1]=Min((1.0/std::sqrt(prior_precision(j,j))*1e-4),1e-4);

		// 		LOGOUT(n);
		// 		LOGOUT(names);
		// 		LOGOUT(values);
		// 		LOGOUT(mins);
		// 		LOGOUT(maxs);	
		
		components.push_back(new MVN_Cov_Mcmc_Component(subject, string("MVN_Cov_")+num2str(n), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),n,prior_mean,prior_precision));
	      }
	  }
	
	if(flag==3 && opts.spatial_crosscorr.value())
	  {	
	    // mvn cov component for each node
	    names = subject_node.get_names_mvn_sqrt_cov_offdiag();
	    if(names.size()>0)
	      {
		values = subject_node.get_value_mvn_sqrt_cov_offdiag();
		mins = subject_node.get_mins_mvn_sqrt_cov_offdiag();
		maxs = subject_node.get_maxs_mvn_sqrt_cov_offdiag();
		prior_mean.ReSize(values.size()); prior_mean=0;
		ColumnVector prior_precision_col(values.size());
		prior_precision_col=1e-10; prior_precision<<diag(prior_precision_col);
		gradient_step_sizes.clear(); gradient_step_sizes.resize(values.size(),1e-4);
		for(unsigned int j=1; j<=gradient_step_sizes.size(); j++)      
		  gradient_step_sizes[j-1]=Min((1.0/std::sqrt(prior_precision(j,j))*1e-4),1e-4);

		// 		    		LOGOUT(n);
		// 		    		LOGOUT(names);
		// 		    		LOGOUT(values);
		// 		    		LOGOUT(mins);
		// 		    		LOGOUT(maxs);	

		components.push_back(new MVN_Cov_Offdiag_Mcmc_Component(subject, string("MVN_Cov_Offdiag")+num2str(n), names, values, mins, maxs,gradient_step_sizes, opts.debuglevel.value(), opts.output_samples.value(),n,prior_mean,prior_precision));
		
	      }
	  } 
	  
      }
    
  }

  double Nma_manager::run_mcmc(Subject& subject, const vector<Mcmc_Component*>& components, int nsamps, int burnin, int se, int sample_vb_jump_every, int burnin_vb_jump_every, bool write_report)
  {
    Tracer_Plus trace("Nma_manager::run_mcmc");

    // MCMC stuff    
    LogSingleton::getInstance().str() << "nsamps = " << nsamps << endl;
    LogSingleton::getInstance().str() << "burnin samps = " << burnin << endl;
    LogSingleton::getInstance().str() << "sampleevery = " << se << endl;
    LogSingleton::getInstance().str() << "sample_vb_jump_every = " << sample_vb_jump_every << endl;
    LogSingleton::getInstance().str() << "burnin_vb_jump_every = " << burnin_vb_jump_every << endl;
        
    // MCMC likelihood - used to calc model evidence
    Nma_Mcmc_Log_Likelihood nma_mcmc_log_likelihood(subject);
    Nma_Vb_Component vb_component(components,subject,opts.debuglevel.value());

    //    bool use_hmc=false;

    // setup appropriate mcmc object
    Mcmc* mcmc_ptr;
//     if(use_hmc)
//       mcmc_ptr=new Mcmc_hmc(components,nma_mcmc_log_likelihood,nsamps,burnin,se,opts.debuglevel.value(),false);
//     else
      {	
	mcmc_ptr=new Mcmc_Mh(components,vb_component,nma_mcmc_log_likelihood,nsamps,burnin,se,sample_vb_jump_every,burnin_vb_jump_every,opts.debuglevel.value(),opts.output_samples.value());
      }

    Mcmc& mcmc=*mcmc_ptr;
    
    if(false)
      {
	//////////
	// write data to log dir
	if(subject.is_single_timeseries())
	  {
	    const vector<ColumnVector>& node_data=subject.get_node_data();
	    for(unsigned int n=1; n<=node_data.size(); n++)      
	      {
		write_ascii_matrix(node_data[n-1],LogSingleton::getInstance().appendDir("node_data"+num2str(n)));    
	      }
	  }
	else
	  {
	    if(opts.data_mode.value()==1 || opts.data_mode.value()==3)
	      {
		volume4D<float> voxelwise_data_vol=subject.get_data();
		const Matrix& voxelwise_data=subject.get_voxelwise_data();
		const vector<ColumnVector>&  coords=subject.get_subject_model().get_voxel_coordinates();
		
		for(unsigned int n=0; n<coords.size(); n++) // indexes voxels
		  {
		    RowVector tmp=voxelwise_data.Row(n+1);
		    voxelwise_data_vol.setvoxelts(tmp.t(),int(coords[n](1)),int(coords[n](2)),int(coords[n](3)));		  
		  }
		
		save_volume4D(voxelwise_data_vol, LogSingleton::getInstance().appendDir("voxelwise_data"));
	      }
	  }
	
	/////////
      }
    //    if(opts.debuglevel.value()>0)
    
    if(write_report)
      {
// 	LogSingleton::getInstance().str()<< "Creating initialised data fit report..." << endl;
// 	create_report_fit(subject, "start");
// 	LogSingleton::getInstance().str()<< "done"<< endl;
      }
    
    mcmc.run();
    mcmc.save();

    ////////////////////
    // do plots 
    {

      if(write_report)
	{
	  //	  create_report_fit(subject, "end");

	  // reevaluate forward model
	  Subject_Model& subject_model=subject.get_subject_model();

	  subject_model.evaluate_neuronal_activity();
	  subject_model.evaluate_hrf();
	  subject_model.evaluate_node_bold();
	  
	  if(!subject.is_single_timeseries() && subject.is_decode())
	    subject_model.evaluate_decoded_node_data();
	  if(!subject.is_single_timeseries() && !subject.is_decode())
	    subject_model.evaluate_voxelwise_bold();
	 
	  LogSingleton::getInstance().str()<< "Creating data fit report"<< endl;
// 	  create_report_fit(subject, "end");
 
	  // set all parameters to sample MAP
	  mcmc.set_values_to_sample_map();	
	  
	  // reevaluate forward model
	  //Subject_Model& subject_model=subject.get_subject_model();
	  subject_model.evaluate_neuronal_activity();
	  subject_model.evaluate_hrf();
	  subject_model.evaluate_node_bold();
	  
	  if(!subject.is_single_timeseries() && subject.is_decode())
	    subject_model.evaluate_decoded_node_data();
	  if(!subject.is_single_timeseries() && !subject.is_decode())
	    subject_model.evaluate_voxelwise_bold();
	  
	  create_report_fit(subject, "map");
	  LogSingleton::getInstance().str()<< "Creating MCMC chains report"<< endl;
	  create_report_mcmc_chains(subject, mcmc, nma_mcmc_log_likelihood);
	  LogSingleton::getInstance().str()<< "Creating connectivity tables report"<< endl;
	  create_report_tables(subject, mcmc);
	  LogSingleton::getInstance().str()<< "Creating logfile report"<< endl;
	  create_logfile_report(subject);
	  LogSingleton::getInstance().str()<< "Creating report"<< endl;
	  create_report(subject);
	}
      
    }

    ////////////////////

    delete mcmc_ptr;

    ColumnVector cols=vector2ColumnVector(nma_mcmc_log_likelihood.get_log_likelihood_hist()); 
    double max_log_lik=Maximum(cols);

    return max_log_lik;

  }

  void Nma_manager::write_matrix_as_html_table(Log& htmllog, const Matrix& mean_mat, const Matrix& std_mat, const vector<vector<string> >& names, const string& table_name, const vector<string>& row_names, const vector<string>& col_names, const Matrix& is_amp_mod) const
  {
    Tracer_Plus trace("Nma_manager::write_matrix_as_html_table");

    write_ascii_matrix(mean_mat, LogSingleton::getInstance().appendDir(table_name+"_mean_connections.txt"));
    write_ascii_matrix(std_mat, LogSingleton::getInstance().appendDir(table_name+"_stddev_connections.txt"));

    htmllog << "<TABLE border=\"1\" frame=\"box\" rules=\"all\" summary=\"" << table_name << "\"cellspacing=\"5\" cellpadding=\"5%\">  <CAPTION>" << table_name << "</CAPTION>" << endl;
    // do column headings
    htmllog << "<THEAD align=\"center\"><TR><TD>";
    for(int c=1; c<=mean_mat.Ncols(); c++)
      {
	htmllog << "<TD><em> " << col_names[c-1] << "</em>";
      }
    htmllog << "</THEAD><TBODY align=\"center\">" << endl;

    // do table
    for(int r=1; r<=mean_mat.Nrows(); r++)
      {
	htmllog << "<TR>" << " <TD><em>" << row_names[r-1]<< "</em>";

	for(int c=1; c<=mean_mat.Ncols(); c++)
	  {

	    if(names[r-1][c-1]==string(""))
	      htmllog << "<TD> X<br>X" << endl;
	    else
	      {
		if(is_amp_mod(r,c)==0)		  
		  htmllog << "<TD> <A href=\"" << names[r-1][c-1]+string(".html") << "\">" << mean_mat(r,c) << " </A> <br> ("<< std_mat(r,c) << ")" << endl;
		else
		  {
		    OUT("herea");
		    OUT(names.size());
		    OUT(names[r-1].size());
		    OUT(names[r-1][c-1]);
		    htmllog << "<TD> <A href=\"" << names[r-1][c-1]+string(".html") << "\">" << mean_mat(r,c) << " </A> &nbsp <A href=\"" << names[r-1][c-1]+string("_amp_mod.html") << "\"> AM </A><br> ("<< std_mat(r,c) << ")" << endl;
		  }
	      }
	  }
	htmllog << endl;
	
      }

    htmllog << "</TBODY></TABLE>" << endl;
  }

void Nma_manager::write_amp_mod_report(const Matrix& samples, const string& name)
  {
    string logfilename2=LogSingleton::getInstance().getLogFileName();
    LogSingleton::getInstance().setLogFile(string(name+".html"));
    Log& htmllog = LogSingleton::getInstance();

    // output page for this parameter to link to from table
    miscplot boxplot;
    boxplot.add_xlabel("epochs");    
    boxplot.set_xysize(300,150);
    boxplot.boxplot(samples, LogSingleton::getInstance().appendDir(name+"_bp"), name);
    htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_bp") << ".png\"> <br>"<< endl;
			 
    LogSingleton::getInstance().setLogFile(logfilename2);
  }

  void Nma_manager::write_mcmc_chain_report(const ColumnVector& samples, const string& name)
  {
    string logfilename2=LogSingleton::getInstance().getLogFileName();


    LogSingleton::getInstance().setLogFile(string(name+".html"));
    Log& htmllog = LogSingleton::getInstance();

    // output page for this parameter to link to from table
    miscplot tsplot;
    tsplot.add_xlabel("samples");    
    tsplot.set_xysize(300,150);
    tsplot.set_minmaxscale(1);
    tsplot.timeseries(samples.t(), LogSingleton::getInstance().appendDir(name+"_ts"), name, 0,400,3,0,false);
    htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_ts") << ".png\">"<< endl;
    if(opts.output_samples.value())
      write_ascii_matrix(samples, LogSingleton::getInstance().appendDir(name+"_ts.txt"));

    miscplot histplot;
    histplot.add_xlabel("");    
    histplot.set_xysize(300,150);
    histplot.histogram(samples.t(), LogSingleton::getInstance().appendDir(name+"_hist"), name);
    htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_hist") << ".png\">"<< endl;    

    LogSingleton::getInstance().setLogFile(logfilename2);
  }

  void Nma_manager::create_report_tables(Subject& subject)
  {
  }

  void Nma_manager::create_report_tables(Subject& subject, const Mcmc& mcmc)
  {
    Tracer_Plus trace("Nma_manager::create_report_tables");

    // setup html report file
    string logfilename=LogSingleton::getInstance().getLogFileName();
    LogSingleton::getInstance().setLogFile(string("report_tables.html"));
    LogSingleton::getInstance().set_stream_to_cout(false);
    Log& htmllog = LogSingleton::getInstance();
	
    htmllog << "<HTML> " << endl
	    << "<TITLE>NMA connectivities for " << subject.get_name() << "</TITLE>" << endl
	    << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	    << "/doc/images/fsl-bg.jpg\">" << endl 
	    << "<hr><CENTER><H1>NMA connectivities for<br>" << subject.get_name() << " </H1></CENTER>"<< endl
	    << "<hr><p>" << endl;

    // output evidence
    htmllog << "MEN_1=" << mcmc.get_model_evidence_normalisation() << ";" << endl;
    htmllog << " UME_1=" << mcmc.get_unnormalised_model_evidence() << ";<br>" << endl;
    htmllog << "where UME is Unnormalised Model Evidence, and MEN is Model Evidence Normalisation," << endl;
    htmllog << "and model evidence=exp(MEN)*UME <br>" << endl;
    htmllog << "(so bayes factor between model 1 and model 2 is: exp(MEN_1-MEN_2)*UME_1/UME_2<br>" << endl;
    htmllog << "<br>AIC=" << mcmc.get_aic() << ";" << endl;
    htmllog << "<br>BIC=" << mcmc.get_bic() << ";" << endl;
    htmllog << "<br><br>Tables below show mean posterior values with standard deviations in brackets.<br><br>" << endl;

    // go through a,b,c,d params and create tables
    int nnodes=subject.get_subject_model().get_subject_nodes().size();    
    int nstim=subject.get_subject_model().get_stimuli().size();
    const Model& model=subject.get_subject_model().get_model();
    int component_number, parameter_number;

    // matrix A
    htmllog << "<br>A matrix is the influence of the column node on the activity in the row node.<br><br>" << endl;
    Matrix mean_matA(nnodes,nnodes); mean_matA=0;
    Matrix std_matA(nnodes,nnodes); std_matA=0;
    vector<vector<string> > namesA(nnodes);
    Matrix is_amp_modA(nnodes,nnodes); is_amp_modA=0;
    for(int n=1; n<=nnodes; n++)
      {
	namesA[n-1].resize(nnodes, string(""));
	for(unsigned int i=1; i<=model.get_marker_a()[n-1].size(); i++)
	  {
	    // identify which is the relevant parameter in the MCMC object
	    string name=string("a_"+num2str(n)+"_"+num2str(model.get_marker_a()[n-1][i-1]));	    
	    mcmc.find_parameter(name, component_number, parameter_number);
	    const Mcmc_Component& component=*(mcmc.get_components()[component_number]);

	    // now extract pertinent info about the parameter from the MCMC object
	    namesA[n-1][model.get_marker_a()[n-1][i-1]-1]=component.get_parameter_names()[parameter_number];
	    ColumnVector samples=vector2ColumnVector(component.get_samples()[parameter_number]);
	    mean_matA(n,model.get_marker_a()[n-1][i-1])=mean(samples).AsScalar();
	    std_matA(n,model.get_marker_a()[n-1][i-1])=sqrt(var(samples).AsScalar());

	    write_mcmc_chain_report(samples, name);	    
	  }
      }
    write_matrix_as_html_table(htmllog, mean_matA, std_matA, namesA, "A matrix", subject.get_subject_model().get_model().get_node_names(), subject.get_subject_model().get_model().get_node_names(),is_amp_modA);     
    htmllog << "<br>" << endl;

    // matrix B
    int count=0;
    for(int n=1; n<=nnodes; n++)
      count+=model.get_marker_b()[n-1].size()>0;
    if(count>0)
      {
	htmllog << "<br>B matrix for a node is the influence of the modulation of the row node by the column stimulus on the activity in that node.<br><br>" << endl;
	for(int n=1; n<=nnodes; n++)
	  {
	    Matrix mean_matB(nnodes,nstim); mean_matB=0;
	    Matrix std_matB(nnodes,nstim); std_matB=0;
	    vector<vector<string> > namesB(nnodes, vector<string>(nstim,""));

	    Matrix is_amp_modB(nnodes,nstim); is_amp_modB=0;

	    if(model.get_marker_b()[n-1].size()>0)
	      {
		for(unsigned int i=1; i<=model.get_marker_b()[n-1].size(); i++)
		  {    
		    // identify which is the relevant parameter in the MCMC object
		    int stim_index=model.get_marker_b()[n-1][i-1].second;  
		    int node_index=model.get_marker_b()[n-1][i-1].first;
		    string name=string("b_"+num2str(n)+"_"+num2str(node_index)+"_"+num2str(stim_index));
		    mcmc.find_parameter(name, component_number, parameter_number);
		    const Mcmc_Component& component=*(mcmc.get_components()[component_number]);
		    
		    // now extract pertinent info about the parameter from the MCMC object		    
		    namesB[node_index-1][stim_index-1]=component.get_parameter_names()[parameter_number];
		    ColumnVector samples=vector2ColumnVector(component.get_samples()[parameter_number]);
		    mean_matB(node_index,stim_index)=mean(samples).AsScalar();
		    std_matB(node_index,stim_index)=sqrt(var(samples).AsScalar());
		    
		    write_mcmc_chain_report(samples, name);

		    // see if it is amp mod
		    if(subject.get_subject_model().get_stimuli()[stim_index-1].size()>1)
		      {
			is_amp_modB(node_index,stim_index)=1;
			Matrix samples(mcmc.get_nsamps(), subject.get_subject_model().get_stimuli()[stim_index-1].size());
			for(unsigned int e=1; e<=subject.get_subject_model().get_stimuli()[stim_index-1].size(); e++)
			  {   	      			  
			    // identify which is the relevant parameter in the MCMC object
			    string name_amp_mod=string("b_"+num2str(n)+"_"+num2str(node_index)+"_"+num2str(stim_index)+"_"+num2str(e));
			    mcmc.find_parameter(name_amp_mod, component_number, parameter_number);
			    const Mcmc_Component& component=*(mcmc.get_components()[component_number]);
			    
			    // now extract pertinent info about the parameter from the MCMC object		    
			    samples.Column(e)=vector2ColumnVector(component.get_samples()[parameter_number]);		
			  }
			
			write_amp_mod_report(samples, name+"_amp_mod");
   
		      }		  
		  }
		
		write_matrix_as_html_table(htmllog, mean_matB, std_matB, namesB, string("B matrix into node ")+subject.get_subject_model().get_model().get_node_names()[n-1], subject.get_subject_model().get_model().get_node_names(), subject.get_subject_model().get_model().get_stimuli_names(),is_amp_modB); 
		htmllog << "<br>" << endl;     
	      }
	  }
      }

    // matrix C
    htmllog << "<br>C matrix is the influence of the column stimulus on the activity in the row node.<br><br>" << endl;
    Matrix mean_matC(nnodes,nstim); mean_matC=0;
    Matrix std_matC(nnodes,nstim); std_matC=0;
    vector<vector<string> > namesC(nnodes);
    Matrix is_amp_modC(nnodes,nstim); is_amp_modC=0;
    for(int n=1; n<=nnodes; n++)
      {
	namesC[n-1].resize(nstim, string(""));
	for(unsigned int i=1; i<=model.get_marker_c()[n-1].size(); i++)
	  {
	    int stim_index=model.get_marker_c()[n-1][i-1]; 

	    // identify which is the relevant parameter in the MCMC object
	    string name=string("c_"+num2str(n)+"_"+num2str(stim_index));	    
	    mcmc.find_parameter(name, component_number, parameter_number);
	    const Mcmc_Component& component=*(mcmc.get_components()[component_number]);

	    // now extract pertinent info about the parameter from the MCMC object
	    namesC[n-1][stim_index-1]=component.get_parameter_names()[parameter_number];
	    ColumnVector samples=vector2ColumnVector(component.get_samples()[parameter_number]);
	    mean_matC(n,stim_index)=mean(samples).AsScalar();
	    std_matC(n,stim_index)=sqrt(var(samples).AsScalar());

	    write_mcmc_chain_report(samples, name);

	    // see if it is amp mod
	    if(subject.get_subject_model().get_stimuli()[stim_index-1].size()>1)
	      {
		int node_index=n;
		is_amp_modC(node_index,stim_index)=1;
		Matrix samples(mcmc.get_nsamps(), subject.get_subject_model().get_stimuli()[stim_index-1].size());
		for(unsigned int e=1; e<=subject.get_subject_model().get_stimuli()[stim_index-1].size(); e++)
		  {   	      			  
		    // identify which is the relevant parameter in the MCMC object
		    string name_amp_mod=string("c_"+num2str(node_index)+"_"+num2str(stim_index)+"_"+num2str(e));
		    mcmc.find_parameter(name_amp_mod, component_number, parameter_number);
		    const Mcmc_Component& component=*(mcmc.get_components()[component_number]);
		    
		    // now extract pertinent info about the parameter from the MCMC object		    
		    samples.Column(e)=vector2ColumnVector(component.get_samples()[parameter_number]);		
		  }
		
		write_amp_mod_report(samples, name+"_amp_mod");		
	      }	
	  }
      }
    write_matrix_as_html_table(htmllog, mean_matC, std_matC, namesC, "C matrix", subject.get_subject_model().get_model().get_node_names(), subject.get_subject_model().get_model().get_stimuli_names(),is_amp_modC);     
    htmllog << "<br>" << endl;

    // matrix D
    count=0;
    for(int n=1; n<=nnodes; n++)
      count+=model.get_marker_d()[n-1].size()>0;
    if(count>0)
      {
	htmllog << "<br>D matrix for a node is the influence of the modulation of the row node by the column node on the activity in that node.<br><br>" << endl;
	for(int n=1; n<=nnodes; n++)
	  {
	    Matrix mean_matD(nnodes,nnodes); mean_matD=0;
	    Matrix std_matD(nnodes,nnodes); std_matD=0;
	    vector<vector<string> > namesD(nnodes, vector<string>(nnodes,""));
	    Matrix is_amp_modD(nnodes,nnodes); is_amp_modD=0;
	    if(model.get_marker_d()[n-1].size()>0)
	      {
		for(unsigned int i=1; i<=model.get_marker_d()[n-1].size(); i++)
		  {    
		    // identify which is the relevant parameter in the MCMC object	
		    int node_index=model.get_marker_d()[n-1][i-1].first;
		    int node_index2=model.get_marker_d()[n-1][i-1].second;  
		    string name=string("d_"+num2str(n)+"_"+num2str(node_index)+"_"+num2str(node_index2));
		    mcmc.find_parameter(name, component_number, parameter_number);
		    const Mcmc_Component& component=*(mcmc.get_components()[component_number]);
		    
		    // now extract pertinent info about the parameter from the MCMC object	
		    
		    namesD[node_index-1][node_index2-1]=component.get_parameter_names()[parameter_number];
		    ColumnVector samples=vector2ColumnVector(component.get_samples()[parameter_number]);
		    mean_matD(node_index,node_index2)=mean(samples).AsScalar();
		    std_matD(node_index,node_index2)=sqrt(var(samples).AsScalar());
		    
		    write_mcmc_chain_report(samples, name);
		  }
		write_matrix_as_html_table(htmllog, mean_matD, std_matD, namesD, string("D matrix into node ")+subject.get_subject_model().get_model().get_node_names()[n-1], subject.get_subject_model().get_model().get_node_names(), subject.get_subject_model().get_model().get_node_names(),is_amp_modD); 
		htmllog << "<br>" << endl;
		
	      }
	  }
      }

    ///////////////////
    // plot stimuli
    htmllog << "<br><br> <H2>Stimuli <br><br> " << endl;
    const vector<vector<ColumnVector> >& stimuli=subject.get_subject_model().get_stimuli(); // nstim*nevents*ntpts
    for(unsigned int s=1; s<=stimuli.size(); s++)
      {
	Matrix mat(stimuli[s-1].size(),stimuli[s-1][0].Nrows());
	for(unsigned int e=1; e<=stimuli[s-1].size(); e++)
	  {
	    mat.Row(e)=stimuli[s-1][e-1].t();
	  }
	string name=subject.get_subject_model().get_model().get_stimuli_names()[s-1];

	miscplot tsplot;
	tsplot.add_xlabel("");    
	tsplot.set_xysize(300,150);
	tsplot.set_minmaxscale(1);
	tsplot.timeseries(mat, LogSingleton::getInstance().appendDir(name+"_ts"),name, 0,400,3,0,false);
	htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_ts") << ".png\">"<< endl;	 	
	write_ascii_matrix(mat, LogSingleton::getInstance().appendDir(name+"_ts.txt"));

	htmllog << "<br> " << endl;
      }

    ///////////////////
    // close log file
    htmllog<< "<HR><FONT SIZE=1>This page produced automatically by "
	   << "nma </A>" 
	   << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	   << "FMRIB Software Library</A>.</FONT>" << endl
	   << "</BODY></HTML>" << endl;    
     
    LogSingleton::getInstance().setLogFile(logfilename);
    LogSingleton::getInstance().set_stream_to_cout(true);
     
  }
  
  void Nma_manager::create_report_mcmc_chains(Subject& subject, const Mcmc& mcmc, const Nma_Mcmc_Log_Likelihood& nma_mcmc_log_likelihood)
  {
    Tracer_Plus trace("Nma_manager::create_report_mcmc_chains");

    ///////////////////
    // setup html report file
    string logfilename=LogSingleton::getInstance().getLogFileName();
    LogSingleton::getInstance().setLogFile(string("report_mcmc_chains.html"));
    LogSingleton::getInstance().set_stream_to_cout(false);
    Log& htmllog = LogSingleton::getInstance();
	
    htmllog << "<HTML> " << endl
	    << "<TITLE>NMA MCMC results for " << subject.get_name() << "</TITLE>" << endl
	    << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	    << "/doc/images/fsl-bg.jpg\">" << endl 
	    << "<hr><CENTER><H1>NMA MCMC results for<br>" << subject.get_name() << " </H1>"<< endl
	    << "<hr><LEFT><p>" << endl;

    ColumnVector cols=vector2ColumnVector(mcmc.get_energy_hist());
    bool energy_ok=true;

    for(int i=1; i<=cols.Nrows(); i++)      
      if(isnan(cols(i)) && cols(i)>1e-32 && cols(i)<1e32) energy_ok=false;

    if((Maximum(cols)-Minimum(cols))<1e-16) energy_ok=false;
   
    if(!isnan(mcmc.get_unnormalised_model_evidence()) && energy_ok)
      //if(energy_ok)  
    {
 	///////////////////////////
 	// output energy plot 
 	ColumnVector cols=vector2ColumnVector(mcmc.get_energy_hist());

	miscplot tsplot;
	tsplot.add_xlabel("samples");    
	tsplot.set_xysize(300,150);
	tsplot.set_minmaxscale(1);
	tsplot.timeseries(cols.t(), LogSingleton::getInstance().appendDir("energy_hist_ts"), "Energy", 0,400,3,0,false);
	htmllog << "<img BORDER=0 SRC=\"energy_hist_ts.png\"> "<< endl;	

	write_ascii_matrix(cols, LogSingleton::getInstance().appendDir(string("energy_hist_ts.txt")));
	
	///////////////////////////
	// output log likelihood plot 
	cols=vector2ColumnVector(nma_mcmc_log_likelihood.get_log_likelihood_hist());    
	{
	  miscplot tsplot;
	  tsplot.add_xlabel("samples");   
	  tsplot.set_xysize(300,150);
	  tsplot.set_minmaxscale(1);
	  tsplot.timeseries(cols.t(), LogSingleton::getInstance().appendDir("log_likelihood_hist_ts"), "log(likelihood)", 0,400,3,0,false);
	}
	htmllog << "<img BORDER=0 SRC=\"log_likelihood_hist_ts.png\"> <br><br><br>"<< endl;
	write_ascii_matrix(cols, LogSingleton::getInstance().appendDir(string("log_likelihood_hist_ts.txt")));
	
	///////////////////////////
	// output log likelihood fit for each ROI plot 
	int nnodes=subject.get_subject_model().get_subject_nodes().size();
	for(int n=1; n<=nnodes; n++)      
	  {
	    ColumnVector cols=nma_mcmc_log_likelihood.get_roi_log_likelihood_hist(n);
	    miscplot tsplot;
	    tsplot.add_xlabel("samples");    
	    tsplot.set_xysize(300,150);
	    tsplot.set_minmaxscale(1);
	    if(Minimum(cols)!=Maximum(cols))
	      tsplot.timeseries(cols.t(), LogSingleton::getInstance().appendDir("log_likelihood_"+num2str(n)+"_hist_ts"), string("Log Likelihood ")+subject.get_subject_model().get_subject_nodes()[n-1]->get_node().get_name(), 0,400,3,0,false);
	    htmllog << "<img BORDER=0 SRC=\"log_likelihood_" << num2str(n) <<"_hist_ts.png\">"<< endl;
	  }
	htmllog << "<br><br><br>" << endl;
      }
    else
      {
	ColumnVector cols=vector2ColumnVector(mcmc.get_energy_hist());	
	write_ascii_matrix(cols, LogSingleton::getInstance().appendDir(string("energy_hist_woo.txt")));
	
      }

    ///////////////////////////
    // loop through components output plots of MCMC chains:        
    const vector<Mcmc_Component*>& components=mcmc.get_components();
    for(unsigned int i=0; i<components.size(); i++)
      {	
	if(components[i]->get_save_out())
	  {	 
	    const vector<string>& parameter_names=components[i]->get_parameter_names();
	    const vector<vector<float> >& samples=components[i]->get_samples();// num_params*num_samps	
      
	for(unsigned int n=0; n< parameter_names.size(); n++)
	  {
	    ColumnVector cols=vector2ColumnVector(samples[n]);

	    miscplot tsplot;
	    tsplot.add_xlabel("samples");    
	    tsplot.set_xysize(300,150);
	    tsplot.set_minmaxscale(1);
	    tsplot.timeseries(cols.t(), LogSingleton::getInstance().appendDir(parameter_names[n]+"_ts"), parameter_names[n], 0,400,3,0,false);
	    htmllog << "<img BORDER=0 SRC=\"" <<  parameter_names[n]+string("_ts") << ".png\">"<< endl;
	    if(opts.output_samples.value())
	      write_ascii_matrix(cols, LogSingleton::getInstance().appendDir(parameter_names[n]+string("_ts.txt")));

	    miscplot histplot;
	    histplot.add_xlabel("");    
	    histplot.set_xysize(300,150);
	    histplot.histogram(cols.t(), LogSingleton::getInstance().appendDir(parameter_names[n]+"_hist"), parameter_names[n]);
	    htmllog << "<img BORDER=0 SRC=\"" <<  parameter_names[n]+string("_hist") << ".png\">"<< endl;
	    
	    htmllog << "<br>" << endl;
 	  }
	  }
      }

    ///////////////////
    // close log file
    htmllog<< "<HR><FONT SIZE=1>This page produced automatically by "
	   << "nma </A>" 
	   << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	   << "FMRIB Software Library</A>.</FONT>" << endl
	   << "</BODY></HTML>" << endl;    
    
    LogSingleton::getInstance().setLogFile(logfilename);
    LogSingleton::getInstance().set_stream_to_cout(true);
  }

  void Nma_manager::create_report(Subject& subject)
  {
    Tracer_Plus trace("Nma_manager::create_report");
    
    ///////////////////
    // setup html report file
    string logfilename=LogSingleton::getInstance().getLogFileName();
    LogSingleton::getInstance().setLogFile(string("report.html"));
    LogSingleton::getInstance().set_stream_to_cout(false);
    Log& htmllog = LogSingleton::getInstance();
	
    htmllog << "<HTML> " << endl
	    << "<TITLE>NMA results for " << subject.get_name() << "</TITLE>" << endl
	    << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	    << "/doc/images/fsl-bg.jpg\">" << endl 
	    << "<hr><CENTER><H1>NMA results for<br>" << subject.get_name() << " </H1>"<< endl
	    << "<hr><p>" << endl;

    htmllog << "<TD> <A href=\"report_tables.html\"> Inferred Connectivities <br> " << endl;
    htmllog << "<TD> <A href=\"report_map.html\"> Maximum a posterior model fit<br>" << endl;
    htmllog << "<TD> <A href=\"report_mcmc_chains.html\"> MCMC chains<br>" << endl;
    htmllog << "<TD> <A href=\"report_logfile.html\"> Log file<br>" << endl;
    
    ///////////////////
    // close log file
    htmllog<< "<HR><FONT SIZE=1>This page produced automatically by "
	   << "nma </A>" 
	   << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	   << "FMRIB Software Library</A>.</FONT>" << endl
	   << "</BODY></HTML>" << endl;    


    LogSingleton::getInstance().setLogFile(logfilename);
    LogSingleton::getInstance().set_stream_to_cout(true);
    
  }
  
  void Nma_manager::create_logfile_report(Subject& subject)
  {
    Tracer_Plus trace("Nma_manager::create_logfile_report");
    
    ///////////////////
    // setup html report file
    string logfilename=LogSingleton::getInstance().getLogFileName();
    LogSingleton::getInstance().setLogFile(string("report_logfile.html"));
    LogSingleton::getInstance().set_stream_to_cout(false);
    Log& htmllog = LogSingleton::getInstance();
	
    htmllog << "<HTML> " << endl
	    << "<TITLE>NMA logfile for " << subject.get_name() << "</TITLE>" << endl
	    << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	    << "/doc/images/fsl-bg.jpg\">" << endl 
	    << "<hr><CENTER><H1>NMA logfile for<br>" << subject.get_name() << " </H1></CENTER>"<< endl
	    << "<hr><p>" << endl;

    // stream logfile into html file
    ifstream logfilein;
    logfilein.open((LogSingleton::getInstance().appendDir(logfilename)).c_str(), ios::in);
    if(logfilein.bad())
      {
	throw Exception((string("Unable to read logfile ")+logfilename).c_str());
      }

    string line;
    while (!logfilein.eof())
    {
      getline (logfilein,line);
      htmllog << line << "<br>" << endl;
    }
   
    logfilein.close();

    ///////////////////
    // close log file
    htmllog<< "<HR><FONT SIZE=1>This page produced automatically by "
	   << "nma </A>" 
	   << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	   << "FMRIB Software Library</A>.</FONT>" << endl
	   << "</BODY></HTML>" << endl;    


    LogSingleton::getInstance().setLogFile(logfilename);
    LogSingleton::getInstance().set_stream_to_cout(true);
    
  }  

  ///////////////////////////////////////////////////

  void fit_model_using_c_matrix_only(Subject& subject, int debuglevel)
  {
    Tracer_Plus trace("fit_model_using_c_matrix_only");

    const Model& model=subject.get_subject_model().get_model();
    Subject_Model& subject_model=subject.get_subject_model();
    unsigned int nnodes=model.get_nodes().size();

    for(unsigned int nod=0; nod<nnodes; nod++)	  
      {
	//	LOGOUT(nod);
	// construct design matrix for this node
	vector<ColumnVector> dm_vec;
	for(unsigned int i=1; i<=subject_model.value_c_amp_mod[nod].size(); i++)
	  {
	    int stim_index=model.marker_c[nod][i-1];
	    //	    LOGOUT(stim_index);
	    for(unsigned int e=1; e<=subject_model.stimuli[stim_index-1].size(); e++)
	      {
		//		LOGOUT(e);
		ColumnVector ev=subject_model.stimuli[stim_index-1][e-1];
		ColumnVector ev_real;
		ColumnVector ev_imag;

		double tmp;

		subject_model.halfcos_hrf.fft_stimulus(ev,ev_real,ev_imag,tmp);

		ColumnVector ffthrf_real;
		ColumnVector ffthrf_imag;
		subject_model.halfcos_hrf.fft_hrf(1, 5, 7, 8, 0.0001, 0.2, ffthrf_real, ffthrf_imag);

		write_ascii_matrix(ffthrf_real,LogSingleton::getInstance().appendDir("hrfreal"));

		//		subject_model.halfcos_hrf.convolve_hrf_fft(ev_real,ev_imag,tmp,subject_model.hrf_fft_real[nod],subject_model.hrf_fft_imag[nod],ev);
		subject_model.halfcos_hrf.convolve_hrf_fft(ev_real,ev_imag,tmp,ffthrf_real,ffthrf_imag,ev);

		// regress out confound EVs
		ev = subject_model.residual_forming_confound_evs*ev;

		// add ev to design matrix
		dm_vec.push_back(ev);
	      }
	  } 
	
	if(dm_vec.size()>0) // node has c inputs
	  {
	    // convert dm_ vec to Matrix:
	    Matrix dm(dm_vec[0].Nrows(),dm_vec.size());
	    for(unsigned int e=1; e<=dm_vec.size(); e++)
	      {
		dm.Column(e)=devar(dm_vec[e-1],1);
	      }

	    if(subject.is_single_timeseries())    
	      {
		const ColumnVector& data=subject.get_node_data()[nod];
		ColumnVector betas = pinv(dm)*devar(data,1);  
		
// 		LOGOUT(LogSingleton::getInstance().appendDir("dm.txt"));

//  		write_ascii_matrix(dm,LogSingleton::getInstance().appendDir("dm.txt"));
//  		write_ascii_matrix(dm*betas,LogSingleton::getInstance().appendDir("fitc_dm.txt"));
//  		write_ascii_matrix(data,LogSingleton::getInstance().appendDir("fitc_data.txt"));				
		betas/=10;	

		int beta_index=1;
		for(unsigned int i=1; i<=subject_model.value_c_amp_mod[nod].size(); i++)
		  {
		    int stim_index=model.marker_c[nod][i-1];
		    
		    vector<float> amp_mod_betas;
		    for(unsigned int e=1; e<=subject_model.stimuli[stim_index-1].size(); e++)
		      {
			amp_mod_betas.push_back(betas(beta_index++));
		      }
		    
		    float beta_mean=mean(vector2ColumnVector(amp_mod_betas)).AsScalar();
		    subject_model.value_c[nod][i-1]=beta_mean;
		
		    if(subject_model.stimuli[stim_index-1].size()!=amp_mod_betas.size())
		      {
			LogSingleton::getInstance().str()<< "subject_model.stimuli[stim_index-1].size()!=amp_mod_betas.size()" << endl;
			exit(1);
		      }
		    for(unsigned int e=1; e<=subject_model.stimuli[stim_index-1].size(); e++)
		      {
			subject_model.value_c_amp_mod[nod][i-1][e-1]=(amp_mod_betas[e-1]-beta_mean);
		      }
		  }
	      }

// 	    else // voxelwise data
// 	      {
// 		const Matrix& voxelwise_data=subject.get_voxelwise_data();
// 		Matrix pinvdm=pinv(dm);

// 		float max_sum_beta_mean=0.0;
// 		vector<vector<float> > max_amp_mod_betas; // num_c_values * num_c_amp_mod_values
// 		vector<float> max_beta_mean; // num_c_values
// 		int max_vox=0;

// 		// loop through voxels finding best fit to model
// 		for(int vox=0; vox<voxelwise_data.Nrows(); vox++)
// 		  {
// 		    if(subject_model.is_voxel_in_roi(vox+1,nod+1))
// 		      {
// 			ColumnVector data=voxelwise_data.Row(vox+1).t();
// 			ColumnVector betas=pinvdm*data;  
			
// 			int beta_index=1;
// 			vector<vector<float> > amp_mod_betas(subject_model.value_c_amp_mod[nod].size()); // num_c_values * num_c_amp_mod_values
// 			vector<float> beta_mean(subject_model.value_c_amp_mod[nod].size()); // num_c_values
			
// 			float sum_beta_mean=0.0;
			
// 			for(unsigned int i=1; i<=subject_model.value_c_amp_mod[nod].size(); i++) // iter through c_values
// 			  {
// 			    int stim_index=model.marker_c[nod][i-1];
			    
// 			    for(unsigned int e=1; e<=subject_model.stimuli[stim_index-1].size(); e++) // iter through c_amp_mod_values
// 			      {
// 				amp_mod_betas[i-1].push_back(betas(beta_index++));
// 			      }
// 			    beta_mean[i-1]=mean(vector2ColumnVector(amp_mod_betas[i-1])).AsScalar();
// 			    sum_beta_mean+=abs(beta_mean[i-1]);			
// 			  }
			
// 			if(abs(sum_beta_mean)>abs(max_sum_beta_mean))
// 			  {
// 			    max_sum_beta_mean=sum_beta_mean;
// 			    max_amp_mod_betas=amp_mod_betas;
// 			    max_beta_mean=beta_mean;
// 			    max_vox=vox;
// 			  }
// 		      }
// 		  }
		
// // 		LOGOUT(max_vox);
// // 		LOGOUT(max_sum_beta_mean);
// // 		LOGOUT(max_beta_mean);

// 		// set c values in model
// 		for(unsigned int i=1; i<=subject_model.value_c_amp_mod[nod].size(); i++)
// 		  {
// 		    subject_model.value_c[nod][i-1]=max_beta_mean[i-1];
// 		    int stim_index=model.marker_c[nod][i-1];
// 		    for(unsigned int e=1; e<=subject_model.stimuli[stim_index-1].size(); e++)
// 		      {
// 			subject_model.value_c_amp_mod[nod][i-1][e-1]=(max_amp_mod_betas[i-1][e-1]-max_beta_mean[i-1]);
// 		      }
// 		  }
		
// // 		// set node mean coordinates at max voxel
// // 		const ColumnVector& max_coords=subject_model.get_subject_nodes()[nod]->get_voxel_coordinates()[max_vox];
// // 		//		LOGOUT(max_coords);
// // 		subject_model.get_subject_nodes()[nod]->set_full_value_mvn_mean_func_space(max_coords);

// // 		SymmetricMatrix sqrt_cov_mat;
// // 		sqrt_cov_mat<<2*IdentityMatrix(3);
// // 		subject_model.get_subject_nodes()[nod]->set_full_value_mvn_sqrt_cov(sqrt_cov_mat);	

// 	      } // end of voxelwise data


	  } // end of if node has any c inputs
      } // end of node

    // reinitialise forward model
    subject_model.evaluate_neuronal_activity();
    subject_model.evaluate_hrf();
    subject_model.evaluate_node_bold();

    if(!subject.is_single_timeseries() && subject.is_decode())
      subject_model.evaluate_decoded_node_data();
    if(!subject.is_single_timeseries() && !subject.is_decode())
	subject_model.evaluate_voxelwise_bold();
  }

  void find_best_vox(Subject& subject, int debuglevel)
  {
    Tracer_Plus trace("find_best_vox");

    bool usemode=true;

    const Model& model=subject.get_subject_model().get_model();
    Subject_Model& subject_model=subject.get_subject_model();
    unsigned int nnodes=model.get_nodes().size();

    //    LOGOUT(nnodes);

    for(unsigned int nod=0; nod<nnodes; nod++)	  
      {
	
	const Subject_Node& subject_node=*(subject_model.get_subject_nodes()[nod]);

	//	LOGOUT(nod);
	// construct design matrix for this node
	vector<ColumnVector> dm_vec;

	ColumnVector node_bold = subject.get_subject_model().get_node_bold()[nod];
	Matrix dm=node_bold;

	dm=dm/sqrt(var(dm).AsScalar());

	const Matrix& voxelwise_data=subject.get_voxelwise_data();
	Matrix pinvdm=pinv(dm);

	// We are going to fit a diagonal 3D MVN distribution to the parameter estimates from a fit of the current model at each voxel

	//////////////////////////////////
	// loop through storing voxelwise fits to current model
	ColumnVector betas(voxelwise_data.Nrows());
	betas=0;
	const volume4D<float>& sub_data=subject.get_data();
	volume<float> betas_newimage(sub_data.xsize(),sub_data.ysize(),sub_data.zsize());
	betas_newimage=0;

	Matrix coords(voxelwise_data.Nrows(),3);
	ColumnVector sumx(int(1+subject_node.get_max_func_space_voxel_coordinate()(1)-subject_node.get_min_func_space_voxel_coordinate()(1))); sumx=0;
	ColumnVector sumy(int(1+subject_node.get_max_func_space_voxel_coordinate()(2)-subject_node.get_min_func_space_voxel_coordinate()(2))); sumy=0;
	ColumnVector sumz(int(1+subject_node.get_max_func_space_voxel_coordinate()(3)-subject_node.get_min_func_space_voxel_coordinate()(3))); sumz=0;

	for(int vox=0; vox<voxelwise_data.Nrows(); vox++)
	  if(subject_model.is_voxel_in_roi(vox+1,nod+1))
	    {
	      ColumnVector data=voxelwise_data.Row(vox+1).t();
	      // 	    data=data/sqrt(var(data).AsScalar());
	      betas(vox+1)=(pinvdm*data).AsScalar();

// 	      if(isnan(betas(vox+1)))
// 		{
// 		  write_ascii_matrix(data,LogSingleton::getInstance().appendDir("data"));
// 		  write_ascii_matrix(dm,LogSingleton::getInstance().appendDir("dm"));
// 		  OUT(vox);
// 		  exit(1);
// 		}

	      coords.Row(vox+1)=subject_model.get_voxel_coordinates()[vox].t();
	      betas_newimage(int(coords(vox+1,1)),int(coords(vox+1,2)),int(coords(vox+1,3)))=betas(vox+1);
	    }

	save_volume(betas_newimage,LogSingleton::getInstance().appendDir("betas"+num2str(nod+1)));

	// end of voxelwise fits to current model
	//////////////////////////////////

	// Now fit a diagonal 3D MVN distribution to the parameter estimates 

	betas=(betas-Minimum(betas))/(sum(betas).AsScalar()); // need all PEs to be positive and looking like a PDF
	// could spatially smooth the betas here before finding the mode

	int best_vox=0;
	float best_beta=0;
	for(int vox=0; vox<voxelwise_data.Nrows(); vox++)
	  if(subject_model.is_voxel_in_roi(vox+1,nod+1))
	    {
	      if(betas(vox+1)>best_beta)
		{
		  best_beta=betas(vox+1);
		  best_vox=vox;
		}
	      betas_newimage(int(coords(vox+1,1)),int(coords(vox+1,2)),int(coords(vox+1,3)))=betas(vox+1);
	    }

	for(int vox=0; vox<voxelwise_data.Nrows(); vox++)
	  if(subject_model.is_voxel_in_roi(vox+1,nod+1))	   
	    {
	      sumx(int(1+coords(vox+1,1)-subject_node.get_min_func_space_voxel_coordinate()(1)))+=betas(vox+1);
	      sumy(int(1+coords(vox+1,2)-subject_node.get_min_func_space_voxel_coordinate()(2)))+=betas(vox+1);
	      sumz(int(1+coords(vox+1,3)-subject_node.get_min_func_space_voxel_coordinate()(3)))+=betas(vox+1);
	    }

// 	OUT(subject_node.get_min_voxel_coordinate());
// 	OUT(subject_node.get_max_voxel_coordinate());

	ColumnVector xs(sumx.Nrows());
	for(int i=1; i<=sumx.Nrows(); i++)
	  xs(i)=i+subject_node.get_min_func_space_voxel_coordinate()(1)-1;
	ColumnVector ys(sumy.Nrows());
	for(int i=1; i<=sumy.Nrows(); i++)
	  ys(i)=i+subject_node.get_min_func_space_voxel_coordinate()(2)-1;
	ColumnVector zs(sumz.Nrows());
	for(int i=1; i<=sumz.Nrows(); i++)
	  zs(i)=i+subject_node.get_min_func_space_voxel_coordinate()(3)-1;

	// save_volume(betas_newimage,LogSingleton::getInstance().appendDir("betas"+num2str(nod+1)));

//  	write_ascii_matrix(betas,LogSingleton::getInstance().appendDir("betas"+num2str(nod+1)));
// 	write_ascii_matrix(xs,LogSingleton::getInstance().appendDir("xs"+num2str(nod+1)));
//  	write_ascii_matrix(ys,LogSingleton::getInstance().appendDir("ys"+num2str(nod+1)));
//  	write_ascii_matrix(zs,LogSingleton::getInstance().appendDir("zs"+num2str(nod+1)));
// 	write_ascii_matrix(sumx,LogSingleton::getInstance().appendDir("sumx"+num2str(nod+1)));
//  	write_ascii_matrix(sumy,LogSingleton::getInstance().appendDir("sumy"+num2str(nod+1)));
//  	write_ascii_matrix(sumz,LogSingleton::getInstance().appendDir("sumz"+num2str(nod+1)));
// 	write_ascii_matrix(coords,LogSingleton::getInstance().appendDir("coords"+num2str(nod+1)));

	// calc node mean coordinates 
	ColumnVector mean_coords(3);
	SymmetricMatrix stddev_coords_mat(3);
	stddev_coords_mat=0;

	// set to mode for now
	mean_coords=coords.Row(best_vox+1).t();

	// do x
	sumx=sumx/(sum(sumx).AsScalar());
	if(!usemode) mean_coords(1)=sum(SP(sumx,xs)).AsScalar();
	
	float tmpx=sqrt(sum(SP(sumx,SP(xs-mean_coords(1),xs-mean_coords(1)))).AsScalar());
	if(tmpx<=0)tmpx=1;	
	stddev_coords_mat(1,1)=tmpx;
	// do y
	sumy=sumy/(sum(sumy).AsScalar());
	if(!usemode) mean_coords(2)=sum(SP(sumy,ys)).AsScalar();
	float tmpy=sqrt(sum(SP(sumy,SP(ys-mean_coords(2),ys-mean_coords(2)))).AsScalar());
	if(tmpy<=0)tmpy=1;
	stddev_coords_mat(2,2)=tmpy;
	// do z
	sumz=sumz/(sum(sumz).AsScalar());
	if(!usemode) mean_coords(3)=sum(SP(sumz,zs)).AsScalar();
	float tmpz=sqrt(sum(SP(sumz,SP(zs-mean_coords(3),zs-mean_coords(3)))).AsScalar());
	if(tmpz<=0)tmpz=1;	
	stddev_coords_mat(3,3)=tmpz;

// 	stddev_coords_mat(1,1)=1.5;
// 	stddev_coords_mat(2,2)=1.5;
// 	stddev_coords_mat(3,3)=1.5;

	stddev_coords_mat(1,2)=0.001;	stddev_coords_mat(1,3)=0.001;	stddev_coords_mat(2,3)=0.001;

	subject_model.get_subject_nodes()[nod]->set_full_value_mvn_mean_func_space(mean_coords);
	
	// set node ROI shape covariance 
  	subject_model.get_subject_nodes()[nod]->set_full_value_mvn_sqrt_cov(stddev_coords_mat);

//  	LOGOUT(nod);
//   	LOGOUT(mean_coords);
// 	//	LOGOUT(betas);
// 	LOGOUT(subject_model.get_subject_nodes()[nod]->get_node().get_name());
// 	LOGOUT(stddev_coords_mat);

// 	LOGOUT(sumx);
// 	LOGOUT(xs);
// 	LOGOUT(tmpx);


      } // end of node

    // reinitialise forward model
    subject_model.evaluate_neuronal_activity();
    subject_model.evaluate_hrf();
    subject_model.evaluate_node_bold();

    if(!subject.is_single_timeseries() && subject.is_decode())
      subject_model.evaluate_decoded_node_data();
    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold();

//     // initialise voxelwise phi's if doing decoding
//     if(subject.is_decode())
//       {
// 	const vector<ColumnVector>& node_bold = subject_model.get_node_bold();
// 	const vector<ColumnVector>& node_data = subject.get_subject_model().get_decoded_node_data();
// 	const vector<ColumnVector>& voxelwise_pvf = subject_model.get_voxelwise_pvf();
// 	const Matrix& voxelwise_data = subject.get_voxelwise_data();

// 	for(unsigned int r=0; r<node_bold.size(); r++)
// 	  {

// 	    double ss = SumSquare(node_bold[r]-node_data[r]);
// 	    double est_node_phi=(node_data[r].Nrows())/ss;

// 	    float sumpvf=0;
// 	    for(int i=1; i<=voxelwise_data.Nrows(); i++)
// 	      {
// 		if(subject_model.is_voxel_in_roi(i,r+1))
// 		  {
// 		    sumpvf+=voxelwise_pvf[r](i);
// 		  }
// 	      }

// 	    for(int i=1; i<=voxelwise_data.Nrows(); i++)
// 	      {
// 		if(subject_model.is_voxel_in_roi(i,r+1))
// 		  {
// 		    vector<float> phivox(1);
// 		    phivox[0]=1.0/( (1.0/est_node_phi)/(voxelwise_data.Nrows()*Sqr(voxelwise_pvf[r](i)/sumpvf)));

// 		    subject_model.set_value_log_phi_every_voxel(phivox,i);
// 		  }
// 	      }
// 	  }
//       }

  }

  
  ///////////////////////////////////////////////////

  Nma_Mcmc_Log_Likelihood::Nma_Mcmc_Log_Likelihood(Subject& psubject)
    : Mcmc_Log_Likelihood(),
      subject(psubject),
      subject_model(psubject.get_subject_model()),
      nnodes(subject.get_subject_model().get_subject_nodes().size()),
      roi_log_likelihood(nnodes)
    {
    }

  const double Nma_Mcmc_Log_Likelihood::calc_energy()
  {
    Tracer_Plus trace("Nma_Mcmc_Log_Likelihood::calc_energy");

    subject_model.evaluate_neuronal_activity();

    if(subject_model.get_haemodynamic_model()=="halfcosine")
      subject_model.evaluate_hrf();      
    
    subject_model.evaluate_node_bold();

    if(!subject.is_single_timeseries() && subject.is_decode())
      subject_model.evaluate_decoded_node_data();
    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold();
    
    double en = calculate_model_energy(subject, 0);
    
    return en;  
  }

  const double Nma_Mcmc_Log_Likelihood::evaluate()
  {
    Tracer_Plus trace("Nma_Mcmc_Log_Likelihood::evaluate");

    double like=0;

    if(subject.is_single_timeseries())
      {
	const vector<ColumnVector>& node_bold = subject.get_subject_model().get_node_bold();
	const vector<ColumnVector>& node_data = subject.get_node_data();

	for(unsigned int r=0; r<node_bold.size(); r++)
	  {
	    double ss = SumSquare(node_bold[r]-node_data[r]);
	    like += -node_bold[r].Nrows()/2.0*std::log(ss);
	  }

      }
    else if(subject.is_decode())
      {
	const vector<ColumnVector>& node_bold = subject.get_subject_model().get_node_bold();
	const vector<ColumnVector>& node_data = subject.get_subject_model().get_decoded_node_data();
	
	for(unsigned int r=0; r<node_bold.size(); r++)
	  {
	    double ss = SumSquare(node_bold[r]-node_data[r]);
	    like += -node_bold[r].Nrows()/2.0*std::log(ss);
	  }

// 	const vector<ColumnVector>& node_bold = subject_model.get_node_bold();
// 	const vector<ColumnVector>& node_data = subject.get_subject_model().get_decoded_node_data();
// 	const vector<ColumnVector>& voxelwise_pvf = subject_model.get_voxelwise_pvf();
// 	const Matrix& voxelwise_data = subject.get_voxelwise_data();
// 	for(unsigned int r=0; r<node_bold.size(); r++)
// 	  {
// 	float sumvar=0;
// 	float sumpvf=0;
	
// 	for(int i=1; i<=voxelwise_data.Nrows(); i++)
// 	  {
// 	    if(subject_model.is_voxel_in_roi(i,r+1))
// 	      {
// 		sumvar+=Sqr(voxelwise_pvf[r](i))/subject_model.get_log_phi_every_voxel()(i);
// 		sumpvf+=voxelwise_pvf[r](i);
// 	      }		
// 	  }
// 	double variance = sumvar/Sqr(sumpvf);
// 	double ss = SumSquare(node_bold[r]-node_data[r]);
// 	like +=  -0.5*node_bold[r].Ncols()*std::log(variance)-0.5*ss/variance;	  
//       }

      }
    else
      {
	const Matrix& voxelwise_bold = subject_model.get_voxelwise_bold();
	const Matrix& voxelwise_data = subject.get_voxelwise_data();
	const vector<ColumnVector>& voxelwise_pvf = subject_model.get_voxelwise_pvf();

	for(int i=1; i<=voxelwise_bold.Nrows(); i++) // loop through voxels
	  {    
	    RowVector tmp=voxelwise_bold.Row(i);
	    double ss = SumSquare(tmp-voxelwise_data.Row(i));
	    
	    double variance;
	    if(subject.is_phi_every_voxel())
	      variance = 1.0/std::exp(subject_model.get_log_phi_voxel()(1)+subject_model.get_log_phi_every_voxel()(i));
	    else
	      variance = 1.0/std::exp(subject_model.get_log_phi_voxel()(1));
	    

	    for(unsigned int n=0; n<voxelwise_pvf.size(); n++)	      
	      {
		variance+=1.0/subject_model.get_phi_node()[n]*Sqr(voxelwise_pvf[n](i));
	      }
	    
	    like += -0.5*tmp.Ncols()*std::log(variance) - 0.5*ss/variance;	    

// 	    LOGOUT(0.5*ss/variance);
// 	    LOGOUT(0.5*tmp.Ncols()*std::log(variance));
	  }


// 	const Matrix& voxelwise_bold = subject.get_subject_model().get_voxelwise_bold();
// 	const Matrix& voxelwise_data = subject.get_voxelwise_data();
		  
// 	for(int i=0; i<voxelwise_bold.Nrows(); i++) //indexes voxels

// 	      RowVector tmp=voxelwise_bold.Row(i);
// 	      double ss = SumSquare(tmp-voxelwise_data.Row(i));
// 	      like += -tmp.Ncols()/2.0*std::log(ss);
      }
    
    // call roi ones as well
    for(int r=0; r<nnodes; r++)
      {
	evaluate(r+1);
      }

    log_likelihood=like;

    return like;  
  }

  const double Nma_Mcmc_Log_Likelihood::evaluate(int node_index)
  {
    Tracer_Plus trace("Nma_Mcmc_Log_Likelihood::evaluate");
    
    double like=0;
    int r=node_index-1;

    if(subject.is_single_timeseries())
      {
	const vector<ColumnVector>& node_bold = subject.get_subject_model().get_node_bold();
	const vector<ColumnVector>& node_data = subject.get_node_data();

	double ss = SumSquare(node_bold[r]-node_data[r]);
	like += -node_bold[r].Nrows()/2.0*std::log(ss);
      }
    else if(subject.is_decode())
      {
	const vector<ColumnVector>& node_bold = subject.get_subject_model().get_node_bold();
	const vector<ColumnVector>& node_data = subject.get_subject_model().get_decoded_node_data();
	
	double ss = SumSquare(node_bold[r]-node_data[r]);
	like += -node_bold[r].Nrows()/2.0*std::log(ss);

// 	const vector<ColumnVector>& node_bold = subject_model.get_node_bold();
// 	const vector<ColumnVector>& node_data = subject.get_subject_model().get_decoded_node_data();
// 	const vector<ColumnVector>& voxelwise_pvf = subject_model.get_voxelwise_pvf();
// 	const Matrix& voxelwise_data = subject.get_voxelwise_data();

// 	float sumvar=0;
// 	float sumpvf=0;
	
// 	for(int i=1; i<=voxelwise_data.Nrows(); i++)
// 	  {
// 	    if(subject_model.is_voxel_in_roi(i,r+1))
// 	      {
// 		sumvar+=Sqr(voxelwise_pvf[r](i))/std::exp(subject_model.get_log_phi_every_voxel()(i));
// 		sumpvf+=voxelwise_pvf[r](i);
// 	      }		
// 	  }
// 	double variance = sumvar/Sqr(sumpvf);
// 	double ss = SumSquare(node_bold[r]-node_data[r]);
// 	like +=  0.5*node_bold[r].Ncols()*std::log(variance)+0.5*ss/variance;	  

      }
    else
      {
	const Matrix& voxelwise_bold = subject_model.get_voxelwise_bold();
	const Matrix& voxelwise_data = subject.get_voxelwise_data();
	const vector<ColumnVector>& voxelwise_pvf = subject_model.get_voxelwise_pvf();

	for(int i=1; i<=voxelwise_bold.Nrows(); i++) // loop through voxels
	  {	
	    if(subject_model.is_voxel_in_roi(i,r+1))
	      {      
		RowVector tmp=voxelwise_bold.Row(i);
		double ss = SumSquare(tmp-voxelwise_data.Row(i));
		
		double variance;
		if(subject.is_phi_every_voxel())
		  variance =1.0/std::exp(subject_model.get_log_phi_voxel()(1)+subject_model.get_log_phi_every_voxel()(i));
		else
		  variance =1.0/std::exp(subject_model.get_log_phi_voxel()(1));
		
		for(unsigned int n=0; n<voxelwise_pvf.size(); n++)	      
		  {
		    variance+=1.0/subject_model.get_phi_node()[n]*Sqr(voxelwise_pvf[n](i));
		  }
	      
		like += -0.5*tmp.Ncols()*std::log(variance)-0.5*ss/variance;
	      }
	  }

// 	const Matrix& voxelwise_bold = subject.get_subject_model().get_voxelwise_bold();
// 	const Matrix& voxelwise_data = subject.get_voxelwise_data();
		  
// 	for(int i=0; i<voxelwise_bold.Nrows(); i++) //indexes voxels
// 	  if(subject_model.is_voxel_in_roi(i,r+1))
// 	    {
// 	      RowVector tmp=voxelwise_bold.Row(i);
// 	      double ss = SumSquare(tmp-voxelwise_data.Row(i));
// 	      like += -tmp.Ncols()/2.0*std::log(ss);
// 	    }
      }
    
    roi_log_likelihood[r]=like;

    return like;  
  }
  
  void Nma_Mcmc_Log_Likelihood::sample()
  {   
   Tracer_Plus trace("Nma_Mcmc_Log_Likelihood::sample");

   Mcmc_Log_Likelihood::sample();
     
   roi_log_likelihood_hist.push_back(roi_log_likelihood);
   
  }

  ReturnMatrix Nma_Mcmc_Log_Likelihood::get_roi_log_likelihood_hist(int node_index) const 
  {   
    Tracer_Plus trace("Nma_Mcmc_Log_Likelihood::get_roi_log_likelihood_hist");

    ColumnVector ret(roi_log_likelihood_hist.size());
    for(unsigned int r=0; r<roi_log_likelihood_hist.size(); r++) //indexes samples
      { 
	ret(r+1)=roi_log_likelihood_hist[r][node_index-1];
      }

    ret.Release();
    return ret;
  }

  ///////////////////////////////////////////////////

  void Connection_Mcmc_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("Connection_Mcmc_Component::calc_noise_precision");

    calculate_noise_precision(subject,noise_precision);

    
  }

  void Connection_Mcmc_Component::calc_forward_model(ColumnVector& forward_model)
  {
    Tracer_Plus trace("Connection_Mcmc_Component::calc_forward_model");

    subject_model.evaluate_neuronal_activity();
    
    subject_model.evaluate_node_bold();

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold();

    calculate_forward_model(subject,forward_model);

  }

  const ColumnVector& Connection_Mcmc_Component::get_vb_data() 
  {
    Tracer_Plus trace("Connection_Mcmc_Component::get_vb_data");

     return subject.get_vb_data();    
  }

  double Connection_Mcmc_Component::calc_energy()
  {
    Tracer_Plus trace("Connection_Mcmc_Component::calc_energy");

    subject_model.evaluate_neuronal_activity();
    
    subject_model.evaluate_node_bold();

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold();

    double energy = calculate_model_energy(subject, debuglevel);

    return energy;
  }

  void Connection_Mcmc_Component::store_old()
  {
    Tracer_Plus trace("Connection_Mcmc_Component::store_old");

    if(subject_model.get_haemodynamic_model()=="halfcosine")
      {
	zfft_real_old=subject_model.get_zfft_real(); 
	zfft_imag_old=subject_model.get_zfft_imag(); 
	z_mean_old=subject_model.get_z_mean(); 
      }
    else
      {
	z_old=subject_model.get_z(); 
	node_cbf_old=subject_model.get_node_cbf();
      }

    node_bold_old=subject_model.get_node_bold();

    if(!subject.is_single_timeseries())
      {
	voxelwise_bold_old=subject_model.get_voxelwise_bold();
	voxelwise_pvf_old = subject_model.get_voxelwise_pvf();
	voxelwise_energy_old=subject_model.get_voxelwise_energy();
      }
  }

  void Connection_Mcmc_Component::restore_old()
  {
    Tracer_Plus trace("Connection_Mcmc_Component::restore_old");

    if(subject_model.get_haemodynamic_model()=="halfcosine")
      {
	subject_model.set_zfft_real(zfft_real_old); 
	subject_model.set_zfft_imag(zfft_imag_old); 
	subject_model.set_z_mean(z_mean_old); 
      }
    else
      {
	subject_model.set_z(z_old); 
	subject_model.set_node_cbf(node_cbf_old);
      }

    subject_model.set_node_bold(node_bold_old);

    if(!subject.is_single_timeseries())
      {
	subject_model.set_voxelwise_bold(voxelwise_bold_old);
	subject_model.set_voxelwise_pvf(voxelwise_pvf_old);
	subject_model.set_voxelwise_energy(voxelwise_energy_old);
      }
  }

  ///////////////////////////////////////////////////
   
  void HRF_Mcmc_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("HRF_Mcmc_Component::calc_noise_precision");

    calculate_noise_precision(subject,noise_precision);
    
  }

  void HRF_Mcmc_Component::calc_forward_model(ColumnVector& forward_model)
  {
    Tracer_Plus trace("HRF_Mcmc_Component::calc_forward_model");
    
    subject_model.evaluate_hrf(node_index);
    
    subject_model.evaluate_node_bold();
    
    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold();

    calculate_forward_model(subject,forward_model);
    
  }

  const ColumnVector& HRF_Mcmc_Component::get_vb_data()
  {
    Tracer_Plus trace("HRF_Mcmc_Component::get_vb_data");

     return subject.get_vb_data();    
  }

  double HRF_Mcmc_Component::calc_energy()
  {
    Tracer_Plus trace("HRF_Mcmc_Component::calc_energy");

    subject_model.evaluate_hrf(node_index);
    
    subject_model.evaluate_node_bold(node_index);

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold(node_index);

    double energy = calculate_model_energy(subject, debuglevel);

    return energy;
  }      
      
  void HRF_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("HRF_Mcmc_Component::set_values");
    
    const vector<Subject_Node*>& subject_nodes = subject.get_subject_model().get_subject_nodes();
    
    subject_nodes[node_index-1]->set_value_hrf(values);

  }

  void HRF_Mcmc_Component::store_old()
  {
    Tracer_Plus trace("HRF_Mcmc_Component::store_old");

    hrf_fft_real_old=subject_model.get_hrf_fft_real(); 
    hrf_fft_imag_old=subject_model.get_hrf_fft_imag();

    node_bold_old=subject_model.get_node_bold();

    if(!subject.is_single_timeseries())
      {
	voxelwise_bold_old=subject_model.get_voxelwise_bold();
	voxelwise_pvf_old = subject_model.get_voxelwise_pvf();
	voxelwise_energy_old=subject_model.get_voxelwise_energy();
      }
  }

  void HRF_Mcmc_Component::restore_old()
  {
    Tracer_Plus trace("HRF_Mcmc_Component::restore_old");

    subject_model.set_hrf_fft_real(hrf_fft_real_old); 
    subject_model.set_hrf_fft_imag(hrf_fft_imag_old); 

    subject_model.set_node_bold(node_bold_old);

     if(!subject.is_single_timeseries())
      {
	subject_model.set_voxelwise_bold(voxelwise_bold_old);
	subject_model.set_voxelwise_pvf(voxelwise_pvf_old);
	subject_model.set_voxelwise_energy(voxelwise_energy_old);
      }  
  }

  ///////////////////////////////////////////////////
   
  void Balloon_Mcmc_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("Balloon_Mcmc_Component::calc_noise_precision");

    calculate_noise_precision(subject,noise_precision);

    
  }

  void Balloon_Mcmc_Component::calc_forward_model(ColumnVector& forward_model)
  {
    Tracer_Plus trace("Balloon_Mcmc_Component::calc_forward_model");
    
    subject_model.evaluate_node_bold_balloon(node_index);

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold(node_index);

    calculate_forward_model(subject,forward_model);

    
  }

  const ColumnVector& Balloon_Mcmc_Component::get_vb_data()
  {
    Tracer_Plus trace("Balloon_Mcmc_Component::get_vb_data");

     return subject.get_vb_data();    
  }

  double Balloon_Mcmc_Component::calc_energy()
  {
    Tracer_Plus trace("Balloon_Mcmc_Component::calc_energy");

    subject_model.evaluate_node_bold_balloon(node_index);

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold(node_index);

    double energy = calculate_model_energy(subject, debuglevel);

    return energy;
  }            

  void Balloon_Mcmc_Component::store_old()
  {
    Tracer_Plus trace("Balloon_Mcmc_Component::store_old");

    node_bold_old=subject_model.get_node_bold();

    if(!subject.is_single_timeseries())
      {
	voxelwise_bold_old=subject_model.get_voxelwise_bold();
	voxelwise_pvf_old = subject_model.get_voxelwise_pvf();
	voxelwise_energy_old=subject_model.get_voxelwise_energy();
      } 
  }

  void Balloon_Mcmc_Component::restore_old()
  {
    Tracer_Plus trace("Balloon_Mcmc_Component::restore_old");

    subject_model.set_node_bold(node_bold_old);

    if(!subject.is_single_timeseries())
      {
	subject_model.set_voxelwise_bold(voxelwise_bold_old);
	subject_model.set_voxelwise_pvf(voxelwise_pvf_old);
	subject_model.set_voxelwise_energy(voxelwise_energy_old);
      }  
  }

  ///////////////////////////////////////////////////
   
  ///////////////////////////////////////////////////
   
  void Balloon_Cbf_Mcmc_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("Balloon_Cbf_Mcmc_Component::calc_noise_precision");

    calculate_noise_precision(subject,noise_precision);

    
  }

  void Balloon_Cbf_Mcmc_Component::calc_forward_model(ColumnVector& forward_model)
  {
    Tracer_Plus trace("Balloon_Cbf_Mcmc_Component::calc_forward_model");
    
    subject_model.evaluate_node_cbf_balloon(node_index);
    subject_model.evaluate_node_bold_balloon(node_index);

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold(node_index);

    calculate_forward_model(subject,forward_model);

    
  }

  const ColumnVector& Balloon_Cbf_Mcmc_Component::get_vb_data()
  {
    Tracer_Plus trace("Balloon_Cbf_Mcmc_Component::get_vb_data");

     return subject.get_vb_data();    
  }

  double Balloon_Cbf_Mcmc_Component::calc_energy()
  {
    Tracer_Plus trace("Balloon_Cbf_Mcmc_Component::calc_energy");

    subject_model.evaluate_node_cbf_balloon(node_index);
    subject_model.evaluate_node_bold_balloon(node_index);

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold(node_index);

    double energy = calculate_model_energy(subject, debuglevel);

    return energy;
  }            

  void Balloon_Cbf_Mcmc_Component::store_old()
  {
    Tracer_Plus trace("Balloon_Cbf_Mcmc_Component::store_old");

    node_cbf_old=subject_model.get_node_cbf();
    node_bold_old=subject_model.get_node_bold();

     if(!subject.is_single_timeseries())
      {
	voxelwise_bold_old=subject_model.get_voxelwise_bold();
	voxelwise_pvf_old = subject_model.get_voxelwise_pvf();
	voxelwise_energy_old=subject_model.get_voxelwise_energy();
      }  
  }

  void Balloon_Cbf_Mcmc_Component::restore_old()
  {
    Tracer_Plus trace("Balloon_Cbf_Mcmc_Component::restore_old");

    subject_model.set_node_cbf(node_cbf_old);
    subject_model.set_node_bold(node_bold_old);

    if(!subject.is_single_timeseries())
      {
	subject_model.set_voxelwise_bold(voxelwise_bold_old);
	subject_model.set_voxelwise_pvf(voxelwise_pvf_old);
	subject_model.set_voxelwise_energy(voxelwise_energy_old);
      } 
  }

  void Balloon_Cbf_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("Balloon_Cbf_Mcmc_Component::set_values");
    
    const vector<Subject_Node*>& subject_nodes = subject.get_subject_model().get_subject_nodes();
    
    subject_nodes[node_index-1]->set_value_balloon_cbf(values);

  }


  ///////////////////////////////////////////////////
   
      
  void Balloon_Bold_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("Balloon_Bold_Mcmc_Component::set_values");
    
    const vector<Subject_Node*>& subject_nodes = subject.get_subject_model().get_subject_nodes();
    
    subject_nodes[node_index-1]->set_value_balloon(values);

  }

  ///////////////////////////////////////////////////
      
  void Balloon_Bold2_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("Balloon_Bold2_Mcmc_Component::set_values");
    
    const vector<Subject_Node*>& subject_nodes = subject.get_subject_model().get_subject_nodes();
    
    subject_nodes[node_index-1]->set_value_balloon2(values);

  }

  ///////////////////////////////////////////////////
   
  void MVN_Mcmc_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("MVN_Mcmc_Component::calc_noise_precision");

    calculate_noise_precision(subject,noise_precision);

    
  }

  void MVN_Mcmc_Component::calc_forward_model(ColumnVector& forward_model)
  {
    Tracer_Plus trace("MVN_Mcmc_Component::calc_forward_model");
    
    if(subject.is_single_timeseries())
      throw Exception("Should not be in this function, MVN_Mcmc_Component::calc_fm");
    
    if(subject.is_decode())
      subject_model.evaluate_decoded_node_data(node_index);
    else
      subject_model.evaluate_voxelwise_bold(node_index);

    calculate_forward_model(subject,forward_model);

    
  }

  const ColumnVector& MVN_Mcmc_Component::get_vb_data()
  {
    Tracer_Plus trace("MVN_Mcmc_Component::get_vb_data");

     return subject.get_vb_data();    
  }

  double MVN_Mcmc_Component::calc_energy()
  {
    Tracer_Plus trace("MVN_Mcmc_Component::calc_energy");

    if(subject.is_single_timeseries())
      throw Exception("Should not be in this function, MVN_Mcmc_Component::calc_energy");
    
    if(subject.is_decode())
      subject_model.evaluate_decoded_node_data(node_index);
    else
      subject_model.evaluate_voxelwise_bold(node_index);

    double energy = calculate_model_energy(subject, debuglevel);

    return energy;
  }            

  void MVN_Mcmc_Component::store_old()
  {
    Tracer_Plus trace("MVN_Mcmc_Component::store_old");

    if(!subject.is_single_timeseries())
      {
	voxelwise_bold_old=subject_model.get_voxelwise_bold();
	voxelwise_pvf_old = subject_model.get_voxelwise_pvf();
	voxelwise_energy_old=subject_model.get_voxelwise_energy();
      }  
  }

  void MVN_Mcmc_Component::restore_old()
  {
    Tracer_Plus trace("MVN_Mcmc_Component::restore_old");

    if(!subject.is_single_timeseries())
      {
	subject_model.set_voxelwise_bold(voxelwise_bold_old);
	subject_model.set_voxelwise_pvf(voxelwise_pvf_old);
	subject_model.set_voxelwise_energy(voxelwise_energy_old);
      } 
  }

  ///////////////////////////////////////////////////    

  void MVN_Mean_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("MVN_Mean_Mcmc_Component::set_values");
    
    const vector<Subject_Node*>& subject_nodes = subject.get_subject_model().get_subject_nodes();
    
    subject_nodes[node_index-1]->set_value_mvn_mean_standard_space(values);      
  }

  ///////////////////////////////////////////////////    
      
  void MVN_Cov_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("MVN_Cov_Mcmc_Component::set_values");
    
    const vector<Subject_Node*>& subject_nodes = subject.get_subject_model().get_subject_nodes();
    
    subject_nodes[node_index-1]->set_value_mvn_sqrt_cov(values);      
  }

 ///////////////////////////////////////////////////
      
  void MVN_Cov_Offdiag_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("MVN_Cov_Offdiag_Mcmc_Component::set_values");
    
    const vector<Subject_Node*>& subject_nodes = subject.get_subject_model().get_subject_nodes();
    
    subject_nodes[node_index-1]->set_value_mvn_sqrt_cov_offdiag(values);      
  }

  ///////////////////////////////////////////////////
   
  void Log_Phi_Voxel_Mcmc_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("Log_Phi_Voxel_Mcmc_Component::calc_noise_precision");

    throw Exception("Log_Phi_Voxel_Mcmc_Component::calc_noise_precision");
    
  }

  void Log_Phi_Voxel_Mcmc_Component::calc_forward_model(ColumnVector&)
  {
    Tracer_Plus trace("Log_Phi_Voxel_Mcmc_Component::calc_forward_model");
    
    throw Exception("Log_Phi_Voxel_Mcmc_Component::calc_forward_model");
    
  }

  const ColumnVector& Log_Phi_Voxel_Mcmc_Component::get_vb_data()
  {
    Tracer_Plus trace("Log_Phi_Voxel_Mcmc_Component::get_vb_data");

    throw Exception("Log_Phi_Voxel_Mcmc_Component::get_vb_data");
  }

  double Log_Phi_Voxel_Mcmc_Component::calc_energy()
  {
    Tracer_Plus trace("Log_Phi_Voxel_Mcmc_Component::calc_energy");

    double energy = calculate_model_energy(subject, debuglevel);

    return energy;
  }      
      
  void Log_Phi_Voxel_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("Log_Phi_Voxel_Mcmc_Component::set_values");
    
    subject_model.set_value_log_phi_voxel(values);

  }

  void Log_Phi_Voxel_Mcmc_Component::store_old()
  {
    Tracer_Plus trace("Log_Phi_Voxel_Mcmc_Component::store_old");

  }

  void Log_Phi_Voxel_Mcmc_Component::restore_old()
  {
    Tracer_Plus trace("Log_Phi_Voxel_Mcmc_Component::restore_old");

  }

///////////////////////////////////////////////////
   
  void Log_Phi_Every_Voxel_Mcmc_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("Log_Phi_Every_Voxel_Mcmc_Component::calc_noise_precision");

    throw Exception("Log_Phi_Every_Voxel_Mcmc_Component::calc_noise_precision");

    
  }

  void Log_Phi_Every_Voxel_Mcmc_Component::calc_forward_model(ColumnVector&)
  {
    Tracer_Plus trace("Log_Phi_Every_Voxel_Mcmc_Component::calc_forward_model");
    
    throw Exception("Log_Phi_Every_Voxel_Mcmc_Component::calc_forward_model");

    
  }

  const ColumnVector& Log_Phi_Every_Voxel_Mcmc_Component::get_vb_data()
  {
    Tracer_Plus trace("Log_Phi_Every_Voxel_Mcmc_Component::get_vb_data");

    throw Exception("Log_Phi_Every_Voxel_Mcmc_Component::get_vb_data");
  }

  double Log_Phi_Every_Voxel_Mcmc_Component::calc_energy()
  {
    Tracer_Plus trace("Log_Phi_Every_Voxel_Mcmc_Component::calc_energy");

    evaluate_model_voxelwise_energy(subject, vox, debuglevel);
    double energy=addup_model_voxelwise_energy(subject, debuglevel);
    
    //double energy = calculate_model_energy(subject, debuglevel);

    //    OUT(energy2-energy);

    return energy;
  }      
      
  void Log_Phi_Every_Voxel_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("Log_Phi_Every_Voxel_Mcmc_Component::set_values");
    
    subject_model.set_value_log_phi_every_voxel(values, vox);
  }

  void Log_Phi_Every_Voxel_Mcmc_Component::store_old()
  {
    Tracer_Plus trace("Log_Phi_Every_Voxel_Mcmc_Component::store_old");

    if(!subject.is_single_timeseries())
      voxelwise_energy_old=subject_model.get_voxelwise_energy(vox);
  }

  void Log_Phi_Every_Voxel_Mcmc_Component::restore_old()
  {
    Tracer_Plus trace("Log_Phi_Every_Voxel_Mcmc_Component::restore_old");

    if(!subject.is_single_timeseries())      
      subject_model.set_voxelwise_energy(vox, voxelwise_energy_old);
  }

  ///////////////////////////////////////////////////
   
  void Phi_Node_Mcmc_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("Phi_Node_Mcmc_Component::calc_noise_precision");

    throw Exception("Phi_Node_Mcmc_Component::calc_noise_precision");

    
  }

  void Phi_Node_Mcmc_Component::calc_forward_model(ColumnVector&)
  {
    Tracer_Plus trace("Phi_Node_Mcmc_Component::calc_forward_model");
    
    throw Exception("Phi_Node_Mcmc_Component::calc_forward_model");
  }

  const ColumnVector& Phi_Node_Mcmc_Component::get_vb_data()
  {
    Tracer_Plus trace("Phi_Node_Mcmc_Component::get_vb_data");

    throw Exception("Phi_Node_Mcmc_Component::get_vb_data");
  }

  double Phi_Node_Mcmc_Component::calc_energy()
  {
    Tracer_Plus trace("Phi_Node_Mcmc_Component::calc_energy");

    double energy = calculate_model_energy(subject, debuglevel);

    return energy;
  }      
      
  void Phi_Node_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("Phi_Node_Mcmc_Component::set_values");
    
    subject_model.set_value_phi_node(values,node_index);
  }

  void Phi_Node_Mcmc_Component::store_old()
  {
    Tracer_Plus trace("Phi_Node_Mcmc_Component::store_old");

    voxelwise_energy_old=subject_model.get_voxelwise_energy();
  }

  void Phi_Node_Mcmc_Component::restore_old()
  {
    Tracer_Plus trace("Phi_Node_Mcmc_Component::restore_old");
    
    subject_model.set_voxelwise_energy(voxelwise_energy_old);
  }

  ///////////////////////////////////////////////////


  ///////////////////////////////////////////////////

  void Nma_Vb_Component::calc_noise_precision(ColumnVector& noise_precision)
  {
    Tracer_Plus trace("Nma_Vb_Component::calc_noise_precision");

    calculate_noise_precision(subject,noise_precision);
    
  }

  void Nma_Vb_Component::calc_forward_model(ColumnVector& forward_model)
  {
    Tracer_Plus trace("Nma_Vb_Component::calc_forward_model");

    subject_model.evaluate_neuronal_activity();
    
    subject_model.evaluate_node_bold();

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold();

    calculate_forward_model(subject,forward_model);

  }

  const ColumnVector& Nma_Vb_Component::get_vb_data() 
  {
    Tracer_Plus trace("Nma_Vb_Component::get_vb_data");

     return subject.get_vb_data();    
  }

  double Nma_Vb_Component::calc_energy()
  {
    Tracer_Plus trace("Nma_Vb_Component::calc_energy");

    subject_model.evaluate_neuronal_activity();
    
    subject_model.evaluate_node_bold();

    if(!subject.is_single_timeseries() && !subject.is_decode())
      subject_model.evaluate_voxelwise_bold();

    double energy = calculate_model_energy(subject, debuglevel);

    return energy;
  }

//   void Nma_Vb_Component::store_old()
//   {
//     Tracer_Plus trace("Nma_Vb_Component::store_old");

//     if(subject_model.get_haemodynamic_model()=="halfcosine")
//       {
// 	zfft_real_old=subject_model.get_zfft_real(); 
// 	zfft_imag_old=subject_model.get_zfft_imag(); 
// 	z_mean_old=subject_model.get_z_mean(); 
//       }
//     else
//       {
// 	z_old=subject_model.get_z(); 
// 	node_cbf_old=subject_model.get_node_cbf();
//       }

//     node_bold_old=subject_model.get_node_bold();

//     if(!subject.is_single_timeseries())
//       {
// 	voxelwise_bold_old=subject_model.get_voxelwise_bold();
// 	voxelwise_pvf_old = subject_model.get_voxelwise_pvf();
// 	voxelwise_energy_old=subject_model.get_voxelwise_energy();
//       }
//   }

//   void Nma_Vb_Component::restore_old()
//   {
//     Tracer_Plus trace("Nma_Vb_Component::restore_old");

//     if(subject_model.get_haemodynamic_model()=="halfcosine")
//       {
// 	subject_model.set_zfft_real(zfft_real_old); 
// 	subject_model.set_zfft_imag(zfft_imag_old); 
// 	subject_model.set_z_mean(z_mean_old); 
//       }
//     else
//       {
// 	subject_model.set_z(z_old); 
// 	subject_model.set_node_cbf(node_cbf_old);
//       }

//     subject_model.set_node_bold(node_bold_old);

//     if(!subject.is_single_timeseries())
//       {
// 	subject_model.set_voxelwise_bold(voxelwise_bold_old);
// 	subject_model.set_voxelwise_pvf(voxelwise_pvf_old);
// 	subject_model.set_voxelwise_energy(voxelwise_energy_old);
//       }
//   }

}
