/*  nmamanager.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(nmamanager_h)
#define nmamanager_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "nmaoptions.h"
#include "newimage/newimageall.h"
#include "subject.h"
#include "model.h"
#include "mcmc_mh.h"
#include "vb_component.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

namespace Nma {

  class Nma_Mcmc_Log_Likelihood;

  // Give this class a file containing
  class Nma_manager
    {
    public:

      // constructor
      Nma_manager() : 
	opts(NmaOptions::getInstance())
	{
	}

      void setup();

      void run();

      // Destructor
      virtual ~Nma_manager() {
	//LOGOUT("~Nma_manager start");
	for(int r=0; r<num_subjects; r++) {
	  delete subjects_sts[r]; 
	}
	subjects_sts.clear(); delete model_sts; 
	if(opts.data_mode.value()==1)
	  {
	    for(int r=0; r<num_subjects; r++) {
	      delete subjects_roi[r]; 
	    }
	    subjects_roi.clear(); delete model_roi; 
	    
	  }
	//LOGOUT("~Nma_manager end");
      }
 
    protected:
      void setup_mcmc_components(Subject& subject, int flag, bool do_haemodynamics, vector<Mcmc_Component*>& components);
      double run_mcmc(Subject& subject, const vector<Mcmc_Component*>& components, int nsamps, int burnin, int se, int sample_vb_jump_every, int burnin_vb_jump_every, bool write_report=true);
      void initialising_runs(int s, string dirname, const vector<Mcmc_Component*>& components_sts, vector<Mcmc_Component*>& components_basis_best);

      void create_report(Subject& subject);
      void create_report_tables(Subject& subject);
      void create_report_tables(Subject& subject, const Mcmc& mcmc);
      void create_report_mcmc_chains(Subject& subject, const Mcmc& mcmc, const Nma_Mcmc_Log_Likelihood& nma_mcmc_log_likelihood);
      void create_logfile_report(Subject& subject);

      void write_matrix_as_html_table(Log& htmllog, const Matrix& mean_matA, const Matrix& std_matA, const vector<vector<string> >& namesA, const string& name, const vector<string>& row_names, const vector<string>& col_names, const Matrix& is_amp_mod) const; 
      void write_mcmc_chain_report(const ColumnVector& samples, const string& name);
      void write_amp_mod_report(const Matrix& samples, const string& name);
 
    private:

      const Nma_manager& operator=(Nma_manager& par);     
      Nma_manager(Nma_manager& des);      
     
      NmaOptions& opts;

      vector<Subject*> subjects_sts;

      Model* model_sts;

      vector<Subject*> subjects_roi;

      Model* model_roi;

      vector<Subject*> subjects_basis;

      Model* model_basis;

      vector<Subject*> subjects_basis_best;

      Model* model_basis_best;

      int num_subjects;

    };  


  //////////////////////////////////////////////

  class Nma_Mcmc_Log_Likelihood : public Mcmc_Log_Likelihood
  {   
  public:
    
    Nma_Mcmc_Log_Likelihood(Subject& psubject);

    const double evaluate();
    const double evaluate(int node_index);
    const double calc_energy();

    virtual void sample();

    const vector<vector<double> >& get_roi_log_likelihood_hist() const {return roi_log_likelihood_hist;}
    vector<double> get_roi_log_likelihood() const {return roi_log_likelihood;}

    ReturnMatrix get_roi_log_likelihood_hist(int node_index) const; // returns samples for log likelihood 

    virtual const int get_num_data_points() const
    {
      int num_data_points=0;
      if(subject.is_single_timeseries())
	{
	  const vector<ColumnVector>& node_data=subject.get_node_data();
	  num_data_points=node_data.size()*node_data[0].Nrows();	 
	}
      else if(subject.is_decode())
	{
	  const Matrix& voxelwise_data = subject.get_voxelwise_data();
	  num_data_points=voxelwise_data.Nrows()*voxelwise_data.Ncols();
	}
      else
	{	
	  const Matrix& voxelwise_data = subject.get_voxelwise_data();
	  num_data_points=voxelwise_data.Nrows()*voxelwise_data.Ncols();
	}

      return num_data_points;
    }

  protected:    

    Subject& subject;
    Subject_Model& subject_model;
    
    int nnodes;

    vector<vector<double> > roi_log_likelihood_hist; // nsamples*nnodes
    vector<double> roi_log_likelihood; // nodes

  };

  //////////////////////////////////////////////

 //////////////////////////////////////////////

  const double bvn_energy(float& x1, float& x2);

  class Bvn_Mcmc_Log_Likelihood : public Mcmc_Log_Likelihood
  {   
  public:
    
    Bvn_Mcmc_Log_Likelihood(float& px1, float& px2) :   Mcmc_Log_Likelihood(), x1(px1),x2(px2)
    {}

    virtual const double evaluate() {
      return 1;
    }

    virtual const double calc_energy() {
      Tracer_Plus trace("Bvn_Mcmc_Log_Likelihood::calc_energy");

      return bvn_energy(x1,x2);      
    }

    virtual const int get_num_data_points() const
    {
      
      return 1;
    }

  protected:    

    float& x1;
    float& x2;
  };

  //////////////////////////////////////////////

  class Connection_Mcmc_Component : public Mcmc_Component
  {
    
  public:
    
    Connection_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples, pprior_mean, pprior_precision),
      subject(psubject),
      subject_model(psubject.get_subject_model())
    {
    }
    
    Connection_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision, const ColumnVector& pprior_isard)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples, pprior_mean, pprior_precision, pprior_isard),
      subject(psubject),
      subject_model(psubject.get_subject_model())
    {
    }
    
    virtual void calc_forward_model(ColumnVector& forward_model);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector& get_vb_data();

    virtual double calc_energy();

  protected:    
        
    Subject& subject;
    Subject_Model& subject_model;       

    // store calculations for old values when new values are proposed
    virtual void store_old();
    // restore calculations for old values when proposal is rejected
    virtual void restore_old();

    // stored old stuff when a proposal is made
    vector<ColumnVector> node_cbf_old;
    vector<ColumnVector> node_bold_old;
    Matrix voxelwise_bold_old;
    vector<ColumnVector> voxelwise_pvf_old;
    ColumnVector voxelwise_energy_old;

    vector<ColumnVector> z_old;
    vector<ColumnVector> zfft_real_old; 
    vector<ColumnVector> zfft_imag_old; 
    ColumnVector z_mean_old; 
  };

  //////////////////////////////////////////////

  class A_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    A_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }
    
    
  virtual bool in_bounds() const {
      Tracer_Plus trace("MCMC_manager::in_bounds");

 	bool inbounds=true;

	// sum squares of all a_{ij} pairs needs to be less than nnodes/(nnodes-1) (see appendix a.3 in dcm paper)
	float ss=0.0;
	for(int n=0; n<nparams && inbounds; n++)
	  {
	    inbounds=(values[n]<param_max[n] && values[n]>param_min[n]);
	    ss+=Sqr(values[n]);
	  }

	int nnodes=subject.get_subject_model().get_subject_nodes().size();    
	if(inbounds)
	  inbounds=ss<nnodes/float(nnodes-1);

	return inbounds;
    }

  protected:
    
    virtual void set_values();    
    
  };

  //////////////

  class A_Node_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    A_Node_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pnode_index, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
	node_index(pnode_index)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }
    
    virtual bool in_bounds() const {
      Tracer_Plus trace("MCMC_manager::in_bounds");

 	bool inbounds=true;
	
	for(int n=0; n<nparams && inbounds; n++)
	  {
	    inbounds=(values[n]<param_max[n] && values[n]>param_min[n]);	    
	  }

	if(inbounds)
	  {
	    int nnodes=subject.get_subject_model().get_subject_nodes().size();    

//  	    OUT("++++++++++++++");

// 	    float ss2=0.0;    
// 	    for(int n=0; n<nnodes; n++)
// 	      {
// 		ColumnVector ab(nnodes);
// 		ab=0;
// 		for(unsigned int i=0; i<subject.get_subject_model().get_model().get_marker_a()[n].size(); i++)
// 		  {
// 		    int node_index2=subject.get_subject_model().get_model().get_marker_a()[n][i];
// 		    ab(node_index2)+=subject.get_subject_model().get_value_a()[n][i];
// 		    OUT(node_index2);
// 		    OUT(subject.get_subject_model().get_value_a()[n][i]);
// 		  }
// 		ss2+=SumSquare(ab);
// 	      }
// 	    OUT(ss2);
// 	    OUT("==============");

	    const vector<float> all_values=subject.get_subject_model().get_value_a_vec();
	    // sum squares of all a_{ij} pairs needs to be less than nnodes/(nnodes-1) (see appendix a.3 in dcm paper)
	    float ss=0.0;

// 	    OUT(node_index);

	    for(unsigned int n=0; n<all_values.size() && inbounds; n++)
	      {
		ss+=Sqr(all_values[n]);
		
// 		OUT(n);
// 		OUT(all_values[n]);
	      }
	    
	    inbounds=ss<(nnodes/float(nnodes-1));

// 	    OUT(ss);
// 	    OUT(nnodes/float(nnodes-1));
// 	    OUT(inbounds);
// 	    OUT("++++++++++++++");

	  }

	return inbounds;
    }
    
  protected:

    virtual void set_values();    

    int node_index;
  };

  //////////////////////////////////////////////

  class B_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    B_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }        
    
    virtual bool in_bounds() const {
      Tracer_Plus trace("MCMC_manager::in_bounds");

 	bool inbounds=true;
		
	for(int n=0; n<nparams && inbounds; n++)
	  {
	    inbounds=(values[n]<param_max[n] && values[n]>param_min[n]);
	    
	  }

// 	if(inbounds)
// 	  {
// 	    int nnodes=subject.get_subject_model().get_subject_nodes().size();    
// 	    int nstim=subject.get_subject_model().get_stimuli().size();
	    
// 	    const vector<vector<float> >& all_value_b=subject.get_subject_model().get_value_b(); // nnodes*(inplay nnodes*nstim)
// 	    const vector<vector<pair<int,int> > >& marker_b=subject.get_subject_model().get_model().get_marker_b();

// 	    // sum squares of all b_{ij} pairs from B matrix for each stimulus needs to be less than nnode/(nnodes-1) (see appendix a.3 in dcm paper)

// 	    vector<float> ss(nstim,0.0);
	    
// 	    for(int n=1; n<=nnodes; n++)
// 	      for(unsigned int i=1; i<=all_value_b[n-1].size(); i++)	
// 		{
// 		  int stim_index=marker_b[n-1][i-1].second;	
// 		  ss[stim_index-1]+=Sqr(all_value_b[n-1][i-1]);
// 		}
	    
// 	    for(int e=1; e<=nstim && inbounds; e++)
// 	      inbounds=ss[e-1]<(nnodes/float(nnodes-1));
// 	  }

	return inbounds;
    }

  protected:

    virtual void set_values();    
    
  };

  //////////////////////////////////////////////

  class B_Node_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    B_Node_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pnode_index, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      node_index(pnode_index)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }        
    
    virtual void set_values();    
    
    virtual bool in_bounds() const {
      Tracer_Plus trace("MCMC_manager::in_bounds");
      
 	bool inbounds=true;
		
	for(int n=0; n<nparams && inbounds; n++)
	  {
	    inbounds=(values[n]<param_max[n] && values[n]>param_min[n]);
	    
	  }

	if(inbounds)
	  {
	    int nnodes=subject.get_subject_model().get_subject_nodes().size();    
	    int nstim=subject.get_subject_model().get_stimuli().size();
	    
	    const vector<vector<float> >& all_value_b=subject.get_subject_model().get_value_b(); // nnodes*(inplay nnodes*nstim)
	    const vector<vector<pair<int,int> > >& marker_b=subject.get_subject_model().get_model().get_marker_b();

	    // sum squares of all b_{ij} pairs from B matrix for each stimulus needs to be less than nnodes/(nnodes-1) (see appendix a.3 in dcm paper)

	    vector<float> ss(nstim,0.0);
	    
	    for(int n=1; n<=nnodes; n++)
	      for(unsigned int i=1; i<=all_value_b[n-1].size(); i++)	
		{
		  int stim_index=marker_b[n-1][i-1].second;	
		  ss[stim_index-1]+=Sqr(all_value_b[n-1][i-1]);
		}
	    
	    for(int e=1; e<=nstim && inbounds; e++)
	      inbounds=ss[e-1]<(nnodes/float(nnodes-1));
	  }

	return inbounds;
    }

  protected:

    int node_index;
  };

  //////////////////////////////////////////////

  class B_Amp_Mod_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    B_Amp_Mod_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pnode_index, int pinplay_index, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      node_index(pnode_index),
      inplay_index(pinplay_index)
    {
      //     if(subject.is_decode())
	do_vb=false;

      setup();
    }    
    
  protected:
    
    virtual void set_values();    
    
    int node_index;
    int inplay_index;

  };

  //////////////////////////////////////////////

  class Log_Sigmaa_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    Log_Sigmaa_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }
    
  protected:
    
    virtual void set_values();    
    
  };

  //////////////////////////////////////////////

  class C_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    C_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision, const ColumnVector& pprior_isard)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision, pprior_isard)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }
    
  protected:
    
    virtual void set_values();    
    
  };
  
  //////////////////////////////////////////////

  class C_Node_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    C_Node_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pnode_index, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision, const ColumnVector& pprior_isard)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision, pprior_isard),
	node_index(pnode_index)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }
    
  protected:
    
    virtual void set_values();    

    int node_index;
  };


  //////////////////////////////////////////////

  class C_Amp_Mod_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    C_Amp_Mod_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pnode_index, int pinplay_index, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      node_index(pnode_index),
      inplay_index(pinplay_index)
    {      
      do_vb=false;
      
      setup();
    }
    
    
  protected:
    
    virtual void set_values();
    
    int node_index;
    int inplay_index;
  };

  //////////////////////////////////////////////

  class D_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    D_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }    
    
    
    virtual bool in_bounds() const {
      Tracer_Plus trace("MCMC_manager::in_bounds");
      
 	bool inbounds=true;
		
	for(int n=0; n<nparams && inbounds; n++)
	  {
	    inbounds=(values[n]<param_max[n] && values[n]>param_min[n]);
	    
	  }

	if(inbounds)
	  {
	    int nnodes=subject.get_subject_model().get_subject_nodes().size();    
	    int nstim=subject.get_subject_model().get_stimuli().size();
	    
	    const vector<vector<float> >& all_value_d=subject.get_subject_model().get_value_d(); // nnodes*(inplay nnodes*nstim)
	    const vector<vector<pair<int,int> > >& marker_d=subject.get_subject_model().get_model().get_marker_d();

	    // sum squares of all d_{ij} pairs from D matrix for each stimulus needs to be less than nnodes/(nnodes-1) (see appendix a.3 in dcm paper)

	    vector<float> ss(nstim,0.0);
	    
	    for(int n=1; n<=nnodes; n++)
	      for(unsigned int i=1; i<=all_value_d[n-1].size(); i++)	
		{
		  int stim_index=marker_d[n-1][i-1].second;	
		  ss[stim_index-1]+=Sqr(all_value_d[n-1][i-1]);
		}
	    
	    for(int e=1; e<=nstim && inbounds; e++)
	      inbounds=ss[e-1]<2*nnodes/float(nnodes-1);
	  }

	return inbounds;
    }

  protected:

    virtual void set_values();      

  };

  //////////////////////////////////////////////

  class D_Node_Mcmc_Component : public Connection_Mcmc_Component
  {
    
  public:
    
    D_Node_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pnode_index, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Connection_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      node_index(pnode_index)
    {
      if(subject.is_decode())
	do_vb=false;
      setup();
    }    
    
    virtual bool in_bounds() const {
      Tracer_Plus trace("MCMC_manager::in_bounds");

 	bool inbounds=true;
		
	for(int n=0; n<nparams && inbounds; n++)
	  {
	    inbounds=(values[n]<param_max[n] && values[n]>param_min[n]);
	    
	  }

	if(inbounds)
	  {
	    int nnodes=subject.get_subject_model().get_subject_nodes().size();    
	    int nstim=subject.get_subject_model().get_stimuli().size();
	    
	    const vector<vector<float> >& all_value_d=subject.get_subject_model().get_value_d(); // nnodes*(inplay nnodes*nstim)
	    const vector<vector<pair<int,int> > >& marker_d=subject.get_subject_model().get_model().get_marker_d();

	    // sum squares of all d_{ij} pairs from D matrix for each stimulus needs to be less than nnodes/(nnodes-1) (see appendix a.3 in dcm paper)

	    vector<float> ss(nstim,0.0);
	    
	    for(int n=1; n<=nnodes; n++)
	      for(unsigned int i=1; i<=all_value_d[n-1].size(); i++)	
		{
		  int stim_index=marker_d[n-1][i-1].second;	
		  ss[stim_index-1]+=Sqr(all_value_d[n-1][i-1]);
		}
	    
	    for(int e=1; e<=nstim && inbounds; e++)
	      inbounds=ss[e-1]<nnodes/float(nnodes-1);
	  }

	return inbounds;
    }

  protected:
    
    virtual void set_values();      
    int node_index;
  };

  //////////////////////////////////////////////

  class HRF_Mcmc_Component : public Mcmc_Component
  {
    
  public:
    
    HRF_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      subject(psubject),
      subject_model(psubject.get_subject_model()),
      node_index(pnode_index)
    {
      if(subject.is_decode())
	do_vb=false;
      setup();
    }
    
    virtual void calc_forward_model(ColumnVector& fm);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector& get_vb_data();

    virtual double calc_energy();
    
  protected:
    
    virtual void set_values();
    
    Subject& subject;
    Subject_Model& subject_model;
    
    int node_index;

    // store calculations for old values when new values are proposed
    virtual void store_old();
    // restore calculations for old values when proposal is rejected
    virtual void restore_old();

    // stored old stuff when a proposal is made
    //    vector<ColumnVector> node_cbf_old;
    vector<ColumnVector> node_bold_old;
    Matrix voxelwise_bold_old;
    vector<ColumnVector> voxelwise_pvf_old;
    ColumnVector voxelwise_energy_old;

    vector<ColumnVector> hrf_fft_real_old; 
    vector<ColumnVector> hrf_fft_imag_old; 

  };

//////////////////////////////////////////////

  class Balloon_Mcmc_Component : public Mcmc_Component
  {
    
  public:
    
    Balloon_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      subject(psubject),
      subject_model(psubject.get_subject_model()),
      node_index(pnode_index)
    {
    }
    
    virtual void calc_forward_model(ColumnVector& fm);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector& get_vb_data();

    virtual double calc_energy();
    
  protected:
    
    virtual void set_values() = 0;
    
    Subject& subject;
    Subject_Model& subject_model;
    
    int node_index;

    // store calculations for old values when new values are proposed
    virtual void store_old();
    // restore calculations for old values when proposal is rejected
    virtual void restore_old();

    // stored old stuff when a proposal is made
    vector<ColumnVector> node_bold_old;
    Matrix voxelwise_bold_old;
    vector<ColumnVector> voxelwise_pvf_old;
    ColumnVector voxelwise_energy_old;
  };

  //////////////////////////////////////////////

  class Balloon_Cbf_Mcmc_Component : public Mcmc_Component
  {
    
  public:
    
    Balloon_Cbf_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      subject(psubject),
      subject_model(psubject.get_subject_model()),
      node_index(pnode_index)
    {
      if(subject.is_decode())
	do_vb=false;
      setup();
    }
    
    virtual void calc_forward_model(ColumnVector& fm);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector& get_vb_data();

    virtual double calc_energy();

  protected:
        
    Subject& subject;
    Subject_Model& subject_model;
    
    int node_index;

    // store calculations for old values when new values are proposed
    virtual void store_old();
    // restore calculations for old values when proposal is rejected
    virtual void restore_old();

    virtual void set_values();
 
    // stored old stuff when a proposal is made
    vector<ColumnVector> node_bold_old;
    Matrix voxelwise_bold_old;
    vector<ColumnVector> voxelwise_pvf_old;
    ColumnVector voxelwise_energy_old;
    vector<ColumnVector> node_cbf_old;
  };

  //////////////////////////////////////////////

  class Balloon_Bold_Mcmc_Component : public Balloon_Mcmc_Component
  {
    
  public:
    
    Balloon_Bold_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Balloon_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples, pnode_index, pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;
      setup();
    }
    
  protected:
    
    virtual void set_values();
    
  };

  //////////////////////////////////////////////

  class Balloon_Bold2_Mcmc_Component : public Balloon_Mcmc_Component
  {
    
  public:
    
    Balloon_Bold2_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Balloon_Mcmc_Component(psubject, pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples, pnode_index, pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;
      setup();
    }
    
  protected:
    
    virtual void set_values();
    
  };

 //////////////////////////////////////////////

  class MVN_Mcmc_Component : public Mcmc_Component
  {
    
  public:
    
    MVN_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples, pprior_mean, pprior_precision),
      subject(psubject),
      subject_model(psubject.get_subject_model()),
      node_index(pnode_index)
    {
    }
    
    virtual void calc_forward_model(ColumnVector& fm);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector& get_vb_data();

    virtual double calc_energy();
    
  protected:    
    
    Subject& subject;
    Subject_Model& subject_model;
     
    int node_index;

    // store calculations for old values when new values are proposed
    virtual void store_old();
    // restore calculations for old values when proposal is rejected
    virtual void restore_old();

    // stored old stuff when a proposal is made
    Matrix voxelwise_bold_old;
    vector<ColumnVector> voxelwise_pvf_old;
    ColumnVector voxelwise_energy_old;
  };


  //////////////////////////////////////////////

  class MVN_Mean_Mcmc_Component : public MVN_Mcmc_Component
  {
    
  public:
    
    MVN_Mean_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : MVN_Mcmc_Component(psubject,pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pnode_index,pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;
      setup();
    }
    
    virtual bool in_bounds() const {
      Tracer_Plus trace("MVN_Mean_Mcmc_Component::in_bounds");

      bool inbounds=true;
      
      for(int n=0; n<nparams && inbounds; n++)
	{
	  inbounds=(values[n]<param_max[n] && values[n]>param_min[n]);
	}

      if(inbounds)
	{
	  const Subject_Node& subject_node=*(subject.get_subject_model().get_subject_nodes()[node_index-1]);    
	  const vector<float>& mvn_mean_standard_space_value=subject_node.get_mvn_mean_standard_space_value();
	  const volume<float>& std_mask=subject_node.get_std_mask();
	  
	  inbounds=(std_mask(round(mvn_mean_standard_space_value[0]),round(mvn_mean_standard_space_value[1]),round(mvn_mean_standard_space_value[2]))>0);
 
	}

      return inbounds;
    }

 
  protected:
    virtual void set_values();
    
  };

  //////////////////////////////////////////////

  class MVN_Cov_Mcmc_Component : public MVN_Mcmc_Component
  {
    
  public:
    
    MVN_Cov_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : MVN_Mcmc_Component(psubject,pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pnode_index,pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;

      setup();
    }
    
  protected:
    
    virtual void set_values();
    
  };
 
  //////////////////////////////////////////////

  class MVN_Cov_Offdiag_Mcmc_Component : public MVN_Mcmc_Component
  {
    
  public:
    
    MVN_Cov_Offdiag_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : MVN_Mcmc_Component(psubject,pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, pdebuglevel, poutput_samples,pnode_index,pprior_mean, pprior_precision)
    {
      if(subject.is_decode())
	do_vb=false;
      setup();
    }
    
  protected:
    
    virtual void set_values();
    
  };

  //////////////////////////////////////////////

  class Log_Phi_Voxel_Mcmc_Component : public Mcmc_Component
  {
    
  public:
    
    Log_Phi_Voxel_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      subject(psubject),
      subject_model(psubject.get_subject_model())
    {
      do_vb=false;
      setup();
    }
    
    virtual void calc_forward_model(ColumnVector& fm);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector&  get_vb_data();

    virtual double calc_energy();
    
  protected:
    
    virtual void set_values();
    
    Subject& subject;
    Subject_Model& subject_model;
 
    // store calculations for old values when new values are proposed
    virtual void store_old();
    // restore calculations for old values when proposal is rejected
    virtual void restore_old();

    // stored old stuff when a proposal is made
    // none

  };

  //////////////////////////////////////////////

  class Log_Phi_Every_Voxel_Mcmc_Component : public Mcmc_Component
  {
    
  public:
    
    Log_Phi_Every_Voxel_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision, int pvox)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      subject(psubject),
      subject_model(psubject.get_subject_model()),
      vox(pvox)
    {
      do_vb=false;
      save_out=false;

      setup();
    }
    
    virtual void calc_forward_model(ColumnVector& fm);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector& get_vb_data();

    virtual double calc_energy();
    
  protected:
    
    virtual void set_values();
    
    Subject& subject;
    Subject_Model& subject_model;
 
    int vox;

    // store calculations for old values when new values are proposed
    virtual void store_old();
    // restore calculations for old values when proposal is rejected
    virtual void restore_old();

    // stored old stuff when a proposal is made
    float voxelwise_energy_old;

  };

  //////////////////////////////////////////////

  class Phi_Node_Mcmc_Component : public Mcmc_Component
  {
    
  public:
    
    Phi_Node_Mcmc_Component(Subject& psubject, const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, int pnode_index, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      subject(psubject),
      subject_model(psubject.get_subject_model()),
      node_index(pnode_index)
    {
      do_vb=false;
      setup();
    }
    
    virtual void calc_forward_model(ColumnVector& fm);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector& get_vb_data();

    virtual double calc_energy();
    
  protected:
    
    virtual void set_values();
    
    Subject& subject;
    Subject_Model& subject_model;
    
    int node_index;

    // store calculations for old values when new values are proposed
    virtual void store_old();
    // restore calculations for old values when proposal is rejected
    virtual void restore_old();

    // stored old stuff when a proposal is made
    ColumnVector voxelwise_energy_old;
  };

  ///////////////////////////////////////

  void fit_model_using_c_matrix_only(Subject& subject, int debuglevel);
  void find_best_vox(Subject& subject, int debuglevel);

  ///////////////////////////////////////////////////
      
  inline void A_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("A_Mcmc_Component::set_values");

    subject.get_subject_model().set_value_a(values);
  }
 
  ///////////////////////////////////////////////////
      
  inline void A_Node_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("A_Node_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_a(node_index, values);
  }

  ///////////////////////////////////////////////////
      
  inline void B_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("B_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_b(values);
  }

  ///////////////////////////////////////////////////
      
  inline void B_Amp_Mod_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("B_Amp_Mod_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_b_amp_mod(node_index,inplay_index,values);
  }

  ///////////////////////////////////////////////////

  inline void B_Node_Mcmc_Component::set_values() 
  {
    //   Tracer_Plus trace("B_Node_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_b(node_index, values);
  }

  ///////////////////////////////////////////////////
      
  inline void C_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("C_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_c(values);
  }

  ///////////////////////////////////////////////////
      
  inline void Log_Sigmaa_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("Log_Sigmaa_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_logsigmaa(values);
  }

  ///////////////////////////////////////////////////

  inline void C_Amp_Mod_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("C_Amp_Mod_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_c_amp_mod(node_index,inplay_index,values);
  }

  ///////////////////////////////////////////////////
     
  inline void C_Node_Mcmc_Component::set_values() 
  {
    Tracer_Plus trace("C_Node_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_c(node_index, values);
  }

  ///////////////////////////////////////////////////
   
  inline void D_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("D_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_d(values);
  }

  ///////////////////////////////////////////////////

  inline void D_Node_Mcmc_Component::set_values() 
  {
    //Tracer_Plus trace("D_Node_Mcmc_Component::set_values");
    
    subject.get_subject_model().set_value_d(node_index, values);
  }  

  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////

  class Bvn_Mcmc_Component : public Mcmc_Component
  {
    
  public:
  
    Bvn_Mcmc_Component(float& px1,    float& px2,const string& pname, const vector<string>& pparameter_names, const vector<float>& pinitial_values, const vector<float>& pparam_min, const vector<float>& pparam_max, const vector<float>& pgradient_step_sizes, int pdebuglevel, bool poutput_samples, const ColumnVector& pprior_mean, const SymmetricMatrix& pprior_precision)
      : Mcmc_Component(pname, pparameter_names, pinitial_values, pparam_min, pparam_max, pgradient_step_sizes, 0, pdebuglevel, poutput_samples,pprior_mean, pprior_precision),
      x1(px1),x2(px2)
    {
      do_vb=false;
      setup();
    }
    
    virtual void calc_forward_model(ColumnVector& forward_model){}
    void calc_noise_precision(ColumnVector& prec){}
    virtual const ColumnVector& get_vb_data(){return tmp;}

    virtual double calc_energy() {
   
      return bvn_energy(x1,x2);
    }
    
//     virtual void append_gradient(vector<float>& grad_vec, double energy) {
//       float rho=0.998;
//       float sig1=0.1;
//       float sig2=1;

//       float z1=2*(x1*sig2-rho*x2*sig1)/Sqr(sig1)/sig2;
//       float z2=2*(x2*sig1-rho*x1*sig2)/Sqr(sig2)/sig1;

//       grad_vec.push_back(z1/(2*(1-Sqr(rho))));
//       grad_vec.push_back(z2/(2*(1-Sqr(rho))));
//     }
  
  protected:
    
    virtual void set_values(){
      x1=values[0];
      x2=values[1];
    }  

    // store calculations for old values when new values are proposed
    virtual void store_old() {}
    // restore calculations for old values when proposal is rejected
    virtual void restore_old() {}
  
    float& x1;
    float& x2;

    ColumnVector tmp;

  };


  //////////////////////////////////////////////

  class Nma_Vb_Component : public Vb_Component
  {
    
  public:
    
    Nma_Vb_Component(const vector<Mcmc_Component*>& pmcmc_comps, Subject& psubject,  int pdebuglevel)
      : Vb_Component(pmcmc_comps, pdebuglevel),
      subject(psubject),
      subject_model(psubject.get_subject_model())
    {
    }
    
    virtual void calc_forward_model(ColumnVector& forward_model);
    void calc_noise_precision(ColumnVector& prec); 
    virtual const ColumnVector& get_vb_data();

    virtual double calc_energy();

  protected:    
        
    Subject& subject;
    Subject_Model& subject_model;       

//     // store calculations for old values when new values are proposed
//     virtual void store_old();
//     // restore calculations for old values when proposal is rejected
//     virtual void restore_old();

//     // stored old stuff when a proposal is made
//     vector<ColumnVector> node_cbf_old;
//     vector<ColumnVector> node_bold_old;
//     Matrix voxelwise_bold_old;
//     vector<ColumnVector> voxelwise_pvf_old;
//     ColumnVector voxelwise_energy_old;

//     vector<ColumnVector> z_old;
//     vector<ColumnVector> zfft_real_old; 
//     vector<ColumnVector> zfft_imag_old; 
//     ColumnVector z_mean_old; 
  };


}

#endif







