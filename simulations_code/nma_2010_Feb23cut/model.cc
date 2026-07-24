/*  model.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT  */

#include "model.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include "node.h"
#include <algorithm>
#include <exception>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace std;

namespace Nma {

  ReturnMatrix threecol2ev(const Matrix& threecol, double res, double nsecs)
  {
    // threecol is nevents * 3
    ColumnVector tmp(int(::round(nsecs/res)));
    tmp=0;

    // loop through events
    for(int e=1; e<=threecol.Nrows(); e++)
      {
	int from=int(::round(threecol(e,1)/res));
	int to=int(::round((threecol(e,1)+threecol(e,2))/res));

	if(to>tmp.Nrows()) to=tmp.Nrows();
	if(from<1) from=1;

	tmp.Rows(from,to)=threecol(e,3);
      }

    tmp.Release();
    return tmp;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  Model::Model(bool psingle_timeseries, int pdata_mode, string phaemodynamic_model, const string& pdata_dir, const string& pmodel_dir,const vector<string>& pnode_names, const vector<string>& pstimuli_names, const vector<int>& pstim_amp_mod, const vector<int>& pstim_ard, bool pstim_single_col_format, double pres, float ptr, bool pdecode)
    : single_timeseries(psingle_timeseries),      
      haemodynamic_model(phaemodynamic_model),
      decode(pdecode),
      data_mode(pdata_mode),
      data_dir(pdata_dir),
      model_dir(pmodel_dir),
      node_names(pnode_names),
      stimuli_names(pstimuli_names),
      stim_amp_mod(pstim_amp_mod),
      stim_ard(pstim_ard),
      stim_single_col_format(pstim_single_col_format),
      res(pres),
      tr(ptr)
  {
    setup();
  }

  void Model::setup()
  {
    Tracer_Plus trace("Model::setup");      


    // setup nodes   
    nnodes=node_names.size();
    for(int r=0; r<nnodes; r++)
      {
	bool do_load_mask=(data_mode!=0);
	nodes.push_back(new Node(node_names[r], data_dir, single_timeseries, do_load_mask)); 
      }


   if(haemodynamic_model=="halfcosine")
      {
	// read HRF files
	ColumnVector hrf_prior_mean;
	ColumnVector hrf_prior_min;
	ColumnVector hrf_prior_max;
	
	hrf_prior_mean=read_ascii_matrix(data_dir+"/hrf_prior_mean.txt");
	hrf_prior_min=read_ascii_matrix(data_dir+"/hrf_prior_min.txt");
	hrf_prior_max=read_ascii_matrix(data_dir+"/hrf_prior_max.txt");
      
	// check input compatabilities
	if(hrf_prior_mean.Nrows() != 6 || hrf_prior_min.Nrows()!=6 || hrf_prior_max.Nrows() != 6)
	  {
	    LogSingleton::getInstance().str()<<"HRF prior inputs have wrong dimension"<<endl;
	    throw Exception(string("HRF prior inputs have wrong dimension").data());
	  }

	for(int r=0; r<nnodes; r++)
	  {
	    nodes[r]->set_hrf_prior(hrf_prior_mean,hrf_prior_min,hrf_prior_max);
	  }
      }

    // load matrices
    matA=read_ascii_matrix(model_dir+"/matA.txt");
    matC=read_ascii_matrix(model_dir+"/matC.txt");
    matB.resize(nnodes);
    matD.resize(nnodes);
    for(int n=0; n<nnodes; n++)
      {
	matB[n]=read_ascii_matrix(model_dir+"/matB_into_"+node_names[n]+".txt");
	matD[n]=read_ascii_matrix(model_dir+"/matD_into_"+node_names[n]+".txt");
      }

//     LOGOUT(matA);
//     LOGOUT(matC);

//     LOGOUT(matB[0]);
//     LOGOUT(matB[1]);
//     LOGOUT(matB[2]);

//     LOGOUT(size(matA));
//     LOGOUT(size(matC));
//     LOGOUT(size(matB[0]));
 
    // check input compatabilities
    if(nnodes != matA.Ncols() || nnodes != matA.Nrows())
      {
	LOGOUT(nnodes);LOGOUT(size(matA));
	LogSingleton::getInstance().str()<<"Number of nodes different to A matrix"<<endl;
	throw Exception(string("Number of nodes different to A matrix").data());
      }
    if(nnodes != int(matB.size()) || nnodes != matB[0].Nrows())
      {
	LOGOUT(nnodes);LOGOUT(matB.size());LOGOUT(matB[0].Nrows());
	LogSingleton::getInstance().str()<<"Number of nodes different to B matrix"<< endl;
	throw Exception(string("Number of nodes different to B matrix").data());
      }
    if(nnodes != matC.Nrows())
      {
	LOGOUT(nnodes);LOGOUT(matC.Nrows());
	LogSingleton::getInstance().str() << "Number of nodes different to C matrix" << endl;
	throw Exception(string("Number of nodes different to C matrix").data());
      }
    if(matB[0].Ncols() != matC.Ncols())
      {
	LOGOUT(matB[0].Ncols());LOGOUT(matC.Ncols());
	LogSingleton::getInstance().str() << "Number of stimuli in B matrix different to C matrix" << endl;
	throw Exception(string("Number of stimuli in B matrix different to C matrix").data());
      }
    if(nnodes != int(matD.size()) || nnodes != matD[0].Nrows() || nnodes != matD[0].Ncols())
      {
	LOGOUT(nnodes);LOGOUT(matD.size());LOGOUT(matD[0].Nrows());LOGOUT(matD[0].Ncols());
	LogSingleton::getInstance().str() << "Number of nodes different to D matrix" << endl;
	throw Exception(string("Number of nodes different to D matrix").data());
      }

    // setup markers
    marker_a.resize(nnodes);
    marker_b.resize(nnodes);
    marker_c.resize(nnodes);
    marker_d.resize(nnodes);
    num_a_params=0;
    num_b_params=0;
    num_c_params=0;
    num_d_params=0;

    for(int n=1; n<=nnodes; n++)
      {
	// matrix A
	for(int n2=1; n2<=nnodes; n2++)
	  if(abs(matA(n,n2))>1e-8 && n!=n2)
	    {
	      marker_a[n-1].push_back(n2);
	      num_a_params++;
	    }
	  
	// matrix B
	for(int n2=1; n2<=nnodes; n2++)	  
	  for(int s=1; s<=matB[n-1].Ncols(); s++)   	    
	    if(abs(matB[n-1](n2,s))>1e-8)
	      {
		marker_b[n-1].push_back(make_pair(n2,s));
		num_b_params++;
	      }	      

	// matrix C
	for(int s=1; s<=matC.Ncols(); s++)
	  if(abs(matC(n,s))>1e-8)
	    {
	      marker_c[n-1].push_back(s);
	      num_c_params++;
	    }

	// matrix D
	for(int n2=1; n2<=nnodes; n2++)
	  for(int n3=1; n3<=nnodes; n3++)   	    
	    if(abs(matD[n-1](n2,n3))>1e-8)
	      {
		marker_d[n-1].push_back(make_pair(n2,n3));
		num_d_params++;
	      }	      
      }
  }
 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  Subject_Model::Subject_Model(string psubject_name, const Model& pmodel, string phaemodynamic_model, int pdebuglevel) : model(pmodel), subject_name(psubject_name), halfcos_hrf(model.res,model.tr), haemodynamic_model(phaemodynamic_model), debuglevel(pdebuglevel)
    {
      // setup subject specific nodes
      nnodes=model.nnodes;

      for(int r=0; r<nnodes; r++)
	{
	  subject_nodes.push_back(new Subject_Node(subject_name, *(model.nodes[r])));
	}
    }

  void Subject_Model::copy_neural_connectivity_values(const Subject_Model& sub_model)
  {
    Tracer_Plus trace("Subject_Model::copy_neural_connectivity_values");
    
    value_a=sub_model.value_a;
    value_b=sub_model.value_b;
    value_c=sub_model.value_c;
    value_d=sub_model.value_d;
    value_b_amp_mod=sub_model.value_b_amp_mod;
    value_c_amp_mod=sub_model.value_c_amp_mod;
    value_logsigmaa=sub_model.value_logsigmaa;

    evaluate_neuronal_activity();
    evaluate_hrf();
    evaluate_node_bold();
       
    if(!model.single_timeseries)
      evaluate_voxelwise_bold();

  }

  void Subject_Model::copy_hrf_values(const Subject_Model& sub_model)
  {
    Tracer_Plus trace("Subject_Model::copy_hrf_values");

    for(int r=0; r<nnodes; r++)
      {
	subject_nodes[r]->copy_hrf_values(*(sub_model.get_subject_nodes()[r]));
      }

    if(haemodynamic_model=="halfcosine")
      evaluate_hrf();

    evaluate_node_bold();
   
    if(!model.single_timeseries)
      evaluate_voxelwise_bold();

  }

  void Subject_Model::output_roi_values()
  {
    Tracer_Plus trace("Subject_Model::output_roi_values");

//     for(int r=0; r<nnodes; r++)
//       {
// 	subject_nodes[r]->output_roi_values();
//       }
//    subject_nodes[1]->output_roi_values();
  }

  void Subject_Model::copy_roi_values(const Subject_Model& sub_model)
  {
    Tracer_Plus trace("Subject_Model::copy_roi_values");

    for(int r=0; r<nnodes; r++)
      {
	subject_nodes[r]->copy_roi_values(*(sub_model.get_subject_nodes()[r]));

	if(!model.single_timeseries)
	  {
	    if(model.is_decode())
	      evaluate_decoded_node_data(r+1);
	    else
	      evaluate_voxelwise_bold(r+1);
	  }
      }
  }

  void Subject_Model::setup()
  {
    Tracer_Plus trace("Subject_Model::setup");

    debugcount1=0;
    debugcount2=0;
    halfcos_hrf.set_nsecs(nsecs);

    LOGOUT(ntpts);

    // sigmaa
    ColumnVector tmp=read_ascii_matrix(model.model_dir+"/sigmaa_prior_mean.txt");
    prior_mean_logsigmaa=std::log(tmp(1));
    tmp=read_ascii_matrix(model.model_dir+"/sigmaa_prior_var.txt");
    prior_precision_logsigmaa=1.0/Sqr(std::log(tmp(1)));
    value_logsigmaa=prior_mean_logsigmaa;

    // balloon priors
    for(int r=0; r<nnodes; r++)
      {
	ColumnVector balloon_cbf_means;
	balloon_cbf_means=read_ascii_matrix(model.model_dir+"/"+model.node_names[r]+"_balloon_cbf_mean.txt");
	SymmetricMatrix balloon_cbf_covs;	   
	balloon_cbf_covs << read_ascii_matrix(model.model_dir+"/"+model.node_names[r]+"_balloon_cbf_cov.txt");       

	ColumnVector balloon_means;
	balloon_means=read_ascii_matrix(model.model_dir+"/"+model.node_names[r]+"_balloon_mean.txt");
	SymmetricMatrix balloon_covs;	   
	balloon_covs << read_ascii_matrix(model.model_dir+"/"+model.node_names[r]+"_balloon_cov.txt");
	
	ColumnVector balloon2_means;
	balloon2_means=read_ascii_matrix(model.model_dir+"/"+model.node_names[r]+"_balloon2_mean.txt");
	SymmetricMatrix balloon2_covs;	   
	balloon2_covs << read_ascii_matrix(model.model_dir+"/"+model.node_names[r]+"_balloon2_cov.txt");
	
	subject_nodes[r]->set_balloon_prior(balloon_cbf_means, balloon_cbf_covs,balloon_means, balloon_covs,balloon2_means, balloon2_covs);

      }
    
    // load subject-wise EVs
    int num_evs=model.stimuli_names.size();
    
    // check input compatabilities
    if(num_evs != model.matC.Ncols())
      {
	LOGOUT(num_evs); LOGOUT(model.matC.Ncols());
	LogSingleton::getInstance().str() << "Number of stimulus EVs different to C matrix" << endl;
	throw Exception("Number of stimulus EVs different to C matrix");
      }

    stimuli.resize(num_evs);
    OUT(num_evs);

    for(int r=0; r<num_evs; r++)
      {

	if(model.stim_single_col_format)
	  {
	    if(model.stim_amp_mod[r])
	      throw Exception("stim_single_col_format incompatible with amp modulation");
	    else
	      {
		OUT(model.data_dir+"/" + subject_name + "/" + model.stimuli_names[r] + ".txt");

		ColumnVector tmp_singlecol=read_ascii_matrix(model.data_dir+"/" + subject_name + "/" + model.stimuli_names[r] + ".txt");
	       
// 		tmp_singlecol

		if(tmp_singlecol.Nrows()!=ntpts)
		  {
		    OUT(tmp_singlecol.Nrows());
		    OUT(ntpts);
		    throw Exception("tmp_singlecol.Nrows!=halfcos_hrf.get_nscans()");
		  }

		stimuli[r].push_back(tmp_singlecol);
		
	      }
	    
	  }
	else
	  {
	    Matrix tmp_threecol=read_ascii_matrix(model.data_dir+"/" + subject_name + "/" + model.stimuli_names[r] + ".txt");
	    
	    if(model.stim_amp_mod[r])
	      {		
		for(int e=1; e<=tmp_threecol.Nrows(); e++)
		  {
		    stimuli[r].push_back(threecol2ev(tmp_threecol.Row(e),model.res,nsecs));

		    stimuli[r][e-1]=stimuli[r][e-1]/(Maximum(stimuli[r][e-1])-Minimum(stimuli[r][e-1]));

		    write_ascii_matrix(stimuli[r][e-1],LogSingleton::getInstance().appendDir("stim"+num2str(r+1)+"_"+num2str(e)));
		  }
	      }
	    else
	      {
		stimuli[r].push_back(threecol2ev(tmp_threecol,model.res,nsecs));

		int e=1;
		stimuli[r][e-1]=stimuli[r][e-1]/(Maximum(stimuli[r][e-1])-Minimum(stimuli[r][e-1]));
		

		// 	    ColumnVector tmp=read_ascii_matrix("/Users/woolrich/homedir/matlab/nma_cpp/data_for_spmdcm/stim"+num2str(r+1));
		// 	    stimuli[r].push_back(tmp);
		
		// 	    for (int t=1; t<=stimuli[r][0].Nrows(); t++)
		// 	      if(stimuli[r][0](t)>1) stimuli[r][0](t)=1;
		
	      }	    
	  }
	
	write_ascii_matrix(stimuli[r][0],LogSingleton::getInstance().appendDir("stim"+num2str(r+1)));
			
      }

    confound_evs=read_ascii_matrix(model.data_dir+string("/") + subject_name + "/confound_evs.txt");
    residual_forming_confound_evs=(IdentityMatrix(confound_evs.Nrows())-confound_evs*pinv(confound_evs));

    // check input compatabilities
    if(confound_evs.Nrows() != halfcos_hrf.get_nscans())
      {
	LOGOUT(confound_evs.Nrows()); LOGOUT(halfcos_hrf.get_nscans());
	LogSingleton::getInstance().str() << "Confound matrix has wrong number of time points"<< endl;
	throw Exception("Confound matrix has wrong number of time points");
      }

  }

  void Subject_Model::initialise(bool random_initialise)
  {
    Tracer_Plus trace("Subject_Model::initialise");

    // initialise parameter values

    value_a.resize(nnodes);
    value_b.resize(nnodes);
    value_c.resize(nnodes);
    value_d.resize(nnodes);
    value_b_amp_mod.resize(nnodes);
    value_c_amp_mod.resize(nnodes);

    num_b_amp_mod_params.resize(nnodes);
    num_c_amp_mod_params.resize(nnodes);
 
    num_a_node_params.resize(nnodes);
    num_b_node_params.resize(nnodes);
    num_c_node_params.resize(nnodes);
    num_d_node_params.resize(nnodes);

    LOGOUT(random_initialise);
    for(int n=1; n<=nnodes; n++)
      {
	// matrix A
	for(int n2=1; n2<=nnodes; n2++)
	  if(abs(model.matA(n,n2))>1e-8 && n!=n2)
	    {	   
	      if(random_initialise)
		//		value_a[n-1].push_back(normrnd(1,1,0,0.01));
		value_a[n-1].push_back(0.01);
	      else
		value_a[n-1].push_back(model.matA(n,n2));
	    }
	num_a_node_params[n-1]=value_a[n-1].size();

	// matrix B
	for(int n2=1; n2<=nnodes; n2++)	  
	  for(unsigned int s=1; s<=stimuli.size(); s++)   	    
	    if(abs(model.matB[n-1](n2,s))>1e-8)
	      {
		if(random_initialise)
		  {
		    value_b[n-1].push_back(0.01);   
		    if(stimuli[s-1].size()>1)
		      value_b_amp_mod[n-1].push_back(ColumnVector2vector(normrnd(stimuli[s-1].size(),1,0,0.001)));
		    else
		      {
			vector<float> tmp(1); tmp[0]=0;
			value_b_amp_mod[n-1].push_back(tmp);
		      }
		  }
		else
		  {
		    value_b[n-1].push_back(model.matB[n-1](n2,s));  
		    if(stimuli[s-1].size()>1)
		      value_b_amp_mod[n-1].push_back(ColumnVector2vector(normrnd(stimuli[s-1].size(),1,0,model.matB[n-1](n2,s)*0.001)));
		  else
		    {
		      vector<float> tmp(1); tmp[0]=0;
		      value_b_amp_mod[n-1].push_back(tmp);
		    }
		  }

		if(stimuli[s-1].size()>1)
		  num_b_amp_mod_params[n-1].push_back(stimuli[s-1].size());
		else
		  num_b_amp_mod_params[n-1].push_back(0);
	      }	      
	num_b_node_params[n-1]=value_b[n-1].size();

	// matrix C
	for(unsigned int s=1; s<=stimuli.size(); s++)
	  if(abs(model.matC(n,s))>1e-8)
	    {
	      if(random_initialise)
		{
		  value_c[n-1].push_back(0.01);	      
		  if(stimuli[s-1].size()>1)
		    value_c_amp_mod[n-1].push_back(ColumnVector2vector(normrnd(stimuli[s-1].size(),1,0,0.001)));
		  else
		    {
		      vector<float> tmp(1); tmp[0]=0;
		      value_c_amp_mod[n-1].push_back(tmp);
		    }
		}
	      else
		{
		  value_c[n-1].push_back(model.matC(n,s));	      
		  if(stimuli[s-1].size()>1)
		    value_c_amp_mod[n-1].push_back(ColumnVector2vector(normrnd(stimuli[s-1].size(),1,0,model.matC(n,s)*0.001)));
		  else
		    {
		      vector<float> tmp(1); tmp[0]=0;
		      value_c_amp_mod[n-1].push_back(tmp);
		    }
		}

	      if(stimuli[s-1].size()>1)
		num_c_amp_mod_params[n-1].push_back(stimuli[s-1].size());
	      else
		num_c_amp_mod_params[n-1].push_back(0);
	    }
	num_c_node_params[n-1]=value_c[n-1].size();

	// matrix D
	for(int n2=1; n2<=nnodes; n2++)
	  for(int n3=1; n3<=nnodes; n3++)   	    
	    if(abs(model.matD[n-1](n2,n3))>1e-8)
	      {
		if(random_initialise)
		  {	
		    value_d[n-1].push_back(0.01);
		  }
		else
		  {
		    value_d[n-1].push_back(model.matD[n-1](n2,n3));
		  }
	      }	  
    	num_d_node_params[n-1]=value_d[n-1].size();
		
      }

    // setup containers    
    node_cbf.resize(nnodes);
    node_bold.resize(nnodes);
    node_bold_high_res.resize(nnodes);

   for(int r=0; r<nnodes; r++)
     {
       node_cbf[r].ReSize(ntpts);
       node_bold[r].ReSize(halfcos_hrf.get_nscans());
       node_bold_high_res[r].ReSize(ntpts);
     }

    if(!model.single_timeseries)
      {
	voxelwise_bold.ReSize(get_voxel_coordinates().size(),halfcos_hrf.get_nscans());
	voxelwise_pvf.resize(nnodes);

	voxelwise_energy.ReSize(get_voxel_coordinates().size());

	log_phi_voxel.ReSize(get_voxel_coordinates().size());
	log_phi_voxel=0;
	log_phi_every_voxel.ReSize(get_voxel_coordinates().size());
	log_phi_every_voxel=0.1;
	
	phi_node.resize(nnodes);
	for(int r=0; r<nnodes; r++)
	  {
	    voxelwise_pvf[r].ReSize(get_voxel_coordinates().size());
	    phi_node[r]=1;
	  } 
      }

    //    set_debuglevel(5);

    // initialise forward model
    evaluate_neuronal_activity();
    evaluate_hrf();
    evaluate_node_bold();
    
    if(!model.single_timeseries)
      {	
	if(model.is_decode())
	  evaluate_decoded_node_data();
	else
	  evaluate_voxelwise_bold();	 
      }

    // test to make sure that initialisation is valid:
    for(int n=1; n<=nnodes; n++)
      for(int t=1; t<=node_bold[n-1].Nrows(); t++)
	if(isnan(node_bold[n-1](t)))
	  {
	    LOGOUT("Invalid initialisation");
	    throw Exception("Invalid initialisation");
	  }
  }

  //////
  void Subject_Model::set_value_a(const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_a");

    int index=0;
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=value_a[n-1].size(); i++)
	value_a[n-1][i-1]=pvalues[index++];
   
  }

  const vector<float> Subject_Model::get_value_a_vec() const
  {
    Tracer_Plus trace("Subject_Model::get_value_a_vec");
    
    vector<float> ret(model.get_num_a_params());
    
    int index=0;
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=value_a[n-1].size(); i++)
	ret[index++]=value_a[n-1][i-1];

    return ret;
  }

  /////

  void Subject_Model::set_value_a(int n, const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_a");

    int index=0;
    for(unsigned int i=1; i<=value_a[n-1].size(); i++)
      value_a[n-1][i-1]=pvalues[index++];
   
  }

  const vector<float> Subject_Model::get_value_a_vec(int n) const
  {
    Tracer_Plus trace("Subject_Model::get_value_a_vec");
    
    vector<float> ret(get_num_a_node_params(n));
    
    int index=0;
    for(unsigned int i=1; i<=value_a[n-1].size(); i++)
      ret[index++]=value_a[n-1][i-1];

    return ret;
  }

  ////


  void Subject_Model::set_value_b(const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_b");

    int index=0;
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=value_b[n-1].size(); i++)
	{	
	  value_b[n-1][i-1]=pvalues[index++];
	  int stim_index=model.marker_b[n-1][i-1].second;

	  if(!model.stim_amp_mod[stim_index-1])	   
	      value_b_amp_mod[n-1][i-1][0]=0;
	}
  }

  const vector<float> Subject_Model::get_value_b_vec() const
  {
    Tracer_Plus trace("Subject_Model::get_value_b_vec");
    
    vector<float> ret(model.get_num_b_params());
    
    int index=0;
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=value_b[n-1].size(); i++)
	{	
	  ret[index++]=value_b[n-1][i-1];
	}

    return ret;
  }

  //////

  void Subject_Model::set_value_b_amp_mod(int n, int i, const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_b_amp_mod");

    int index=0;

    int stim_index=model.marker_b[n-1][i-1].second;
    for(unsigned int e=1; e<=stimuli[stim_index-1].size(); e++)
      {
	value_b_amp_mod[n-1][i-1][e-1]=pvalues[index++];	   
      }  
  }

  const vector<float> Subject_Model::get_value_b_amp_mod_vec(int n, int i) const
  {
    Tracer_Plus trace("Subject_Model::get_value_b_amp_mod_vec");
    
    vector<float> ret(get_num_b_amp_mod_params(n,i));
    
    int index=0;

    int stim_index=model.marker_b[n-1][i-1].second;
    for(unsigned int e=1; e<=stimuli[stim_index-1].size(); e++)
      { 
	ret[index++]=value_b_amp_mod[n-1][i-1][e-1];
      }	

    return ret;
  }

  /////

  void Subject_Model::set_value_b(int n, const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_b");

    int index=0;
    for(unsigned int i=1; i<=value_b[n-1].size(); i++)
      {	
	value_b[n-1][i-1]=pvalues[index++];
	int stim_index=model.marker_b[n-1][i-1].second;
	
	if(!model.stim_amp_mod[stim_index-1])	   
	  value_b_amp_mod[n-1][i-1][0]=0;
      }
  }

  const vector<float> Subject_Model::get_value_b_vec(int n) const
  {
    Tracer_Plus trace("Subject_Model::get_value_b_vec");
    
    vector<float> ret(get_num_b_node_params(n));
    
    int index=0;
    for(unsigned int i=1; i<=value_b[n-1].size(); i++)
      {	
	ret[index++]=value_b[n-1][i-1];
      }
    
    return ret;
  }

  ////

  void Subject_Model::set_value_c(const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_c");

    int index=0;
    for(int n=1; n<=nnodes; n++)
	for(unsigned int i=1; i<=value_c[n-1].size(); i++)
	  {	     
	    value_c[n-1][i-1]=pvalues[index++];	   

	    int stim_index=model.marker_c[n-1][i-1];
	    if(!model.stim_amp_mod[stim_index-1])
	      value_c_amp_mod[n-1][i-1][0]=0;
	  }

  }

  const vector<float> Subject_Model::get_value_c_vec() const
  {
    Tracer_Plus trace("Subject_Model::get_value_c_vec");
    
    vector<float> ret(model.get_num_c_params());
    
    int index=0;
    for(int n=1; n<=nnodes; n++)
	for(unsigned int i=1; i<=value_c[n-1].size(); i++)
	  { 
	    ret[index++]=value_c[n-1][i-1];	      
	  }
    
    return ret;
  }

  /////////

  void Subject_Model::set_value_c_amp_mod(int n, int i, const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_c_amp_mod_vec");

    int index=0;

    int stim_index=model.marker_c[n-1][i-1];
    for(unsigned int e=1; e<=stimuli[stim_index-1].size(); e++)  
      {    	      
	value_c_amp_mod[n-1][i-1][e-1]=pvalues[index++];	   
      }
  
  }

  const vector<float> Subject_Model::get_value_c_amp_mod_vec(int n, int i) const
  {
    Tracer_Plus trace("Subject_Model::get_value_c_amp_mod_vec");
    
    vector<float> ret(get_num_c_amp_mod_params(n,i));
    
    int index=0;

    int stim_index=model.marker_c[n-1][i-1];
    for(unsigned int e=1; e<=stimuli[stim_index-1].size(); e++)  
      {  
	ret[index++]=value_c_amp_mod[n-1][i-1][e-1];
      }	  
    
    return ret;
  }

  ///////

  void Subject_Model::set_value_c(int n, const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_c");

    int index=0;
    for(unsigned int i=1; i<=value_c[n-1].size(); i++)
      {	     
	value_c[n-1][i-1]=pvalues[index++];	   
	
	int stim_index=model.marker_c[n-1][i-1];
	if(!model.stim_amp_mod[stim_index-1])
	  value_c_amp_mod[n-1][i-1][0]=0;
      }
    
  }

  const vector<float> Subject_Model::get_value_c_vec(int n) const
  {
    Tracer_Plus trace("Subject_Model::get_value_c_vec");
    
    vector<float> ret(get_num_c_node_params(n));
    
    int index=0;
    for(unsigned int i=1; i<=value_c[n-1].size(); i++)
      { 
	ret[index++]=value_c[n-1][i-1];	      
      }    

    return ret;
  }

  /////////

  void Subject_Model::set_value_d(const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_d");

    int index=0;
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=value_d[n-1].size(); i++)
	value_d[n-1][i-1]=pvalues[index++];	   
    
  }

  const vector<float> Subject_Model::get_value_d_vec() const
  {
    Tracer_Plus trace("Subject_Model::get_value_d_vec");
    
    vector<float> ret(model.get_num_d_params());
    
    int index=0;
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=value_d[n-1].size(); i++)
	ret[index++]=value_d[n-1][i-1];
    
    return ret;
  }

  /////////


  void Subject_Model::set_value_d(int n, const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_d");

    int index=0;
      for(unsigned int i=1; i<=value_d[n-1].size(); i++)
	value_d[n-1][i-1]=pvalues[index++];	   
    
  }

  const vector<float> Subject_Model::get_value_d_vec(int n) const
  {
    Tracer_Plus trace("Subject_Model::get_value_d_vec");
    
    vector<float> ret(get_num_d_node_params(n));
    
    int index=0;
    for(unsigned int i=1; i<=value_d[n-1].size(); i++)
      ret[index++]=value_d[n-1][i-1];
    
    return ret;
  }

  ////////

  void Subject_Model::set_value_logsigmaa(const vector<float>& pvalues) 
  {
    Tracer_Plus trace("Subject_Model::set_value_logsigmaa");

    value_logsigmaa=pvalues[0];
  }

  const vector<float> Subject_Model::get_value_logsigmaa_vec() const
  {
    Tracer_Plus trace("Subject_Model::get_value_logsigmaa_vec");
    
    vector<float> ret(1);    
    
    ret[0]=value_logsigmaa;

    return ret;
  }

  const vector<string> Subject_Model::get_names_a() const
  {
    Tracer_Plus trace("Subject_Model::get_names_a");
    
    vector<string> ret(model.get_num_a_params());
        
    int index=0;
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=model.marker_a[n-1].size(); i++)
	ret[index++]=string("a_"+num2str(n)+"_"+num2str(model.marker_a[n-1][i-1]));
    
    return ret;
  }
 

  const vector<string> Subject_Model::get_names_a(int n) const
  {
    Tracer_Plus trace("Subject_Model::get_names_a");
    
    vector<string> ret(get_num_a_node_params(n));
        
    int index=0;
    for(unsigned int i=1; i<=model.marker_a[n-1].size(); i++)
      ret[index++]=string("a_"+num2str(n)+"_"+num2str(model.marker_a[n-1][i-1]));
    
    return ret;
  }

  const vector<string> Subject_Model::get_names_b() const
  {
    Tracer_Plus trace("Subject_Model::get_names_b");
    
    vector<string> ret(model.get_num_b_params());
    
    int index=0;

    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=value_b[n-1].size(); i++)
	{    
	  int stim_index=model.marker_b[n-1][i-1].second;  
	  int node_index=model.marker_b[n-1][i-1].first;
	  ret[index++]=string("b_"+num2str(n)+"_"+num2str(node_index)+"_"+num2str(stim_index));
	}
    
    return ret;
  }  

  const vector<string> Subject_Model::get_names_b_amp_mod(int n, int i) const
  {
    Tracer_Plus trace("Subject_Model::get_names_b_amp_mod");
    
    vector<string> ret(get_num_b_amp_mod_params(n,i));
    
    if(ret.size()>0)
      {
	int index=0;
	
	int stim_index=model.marker_b[n-1][i-1].second;
	int node_index=model.marker_b[n-1][i-1].first;

	if(stimuli[stim_index-1].size()>1)
	  for(unsigned int e=1; e<=stimuli[stim_index-1].size(); e++)
	    {   	      
	      ret[index++]=string("b_"+num2str(n)+"_"+num2str(node_index)+"_"+num2str(stim_index)+"_"+num2str(e));
	    }
      }
    
    return ret;
  } 

  const vector<string> Subject_Model::get_names_b(int n) const
  {
    Tracer_Plus trace("Subject_Model::get_names_b");
    
    vector<string> ret(get_num_b_node_params(n));
    
    int index=0;
    
    for(unsigned int i=1; i<=value_b[n-1].size(); i++)
      {    
	int stim_index=model.marker_b[n-1][i-1].second;  
	int node_index=model.marker_b[n-1][i-1].first;
	ret[index++]=string("b_"+num2str(n)+"_"+num2str(node_index)+"_"+num2str(stim_index));
      }
    
    return ret;
  }  

  const vector<string> Subject_Model::get_names_c() const
  {
    Tracer_Plus trace("Subject_Model::get_names_c");
    
    //    vector<string> ret(model.get_num_c_params()+1);
    vector<string> ret(model.get_num_c_params());
    
    int index=0;
     for(int n=1; n<=nnodes; n++)
	for(unsigned int i=1; i<=value_c[n-1].size(); i++)
	  {
	    int stim_index=model.marker_c[n-1][i-1];

	    ret[index++]=string("c_"+num2str(n)+"_"+num2str(stim_index));	      
	  }   

//     ret[index++]=string("sigmaa");

    return ret;
  } 

  const vector<string> Subject_Model::get_names_c(int n) const
  {
    Tracer_Plus trace("Subject_Model::get_names_c");
    
    vector<string> ret(get_num_c_node_params(n));
    
    int index=0;
    for(unsigned int i=1; i<=value_c[n-1].size(); i++)
      {
	int stim_index=model.marker_c[n-1][i-1];
	
	ret[index++]=string("c_"+num2str(n)+"_"+num2str(stim_index));	      
      }   

    return ret;
  } 

  const vector<string> Subject_Model::get_names_c_amp_mod(int n, int i) const
  {
    Tracer_Plus trace("Subject_Model::get_names_c_amp_mod");
    
    vector<string> ret(get_num_c_amp_mod_params(n,i));
    
    if(ret.size()>0)
      {
	int index=0;
	
	int stim_index=model.marker_c[n-1][i-1];

	if(stimuli[stim_index-1].size()>1)
	  for(unsigned int e=1; e<=stimuli[stim_index-1].size(); e++)  
	    { 
	      ret[index++]=string("c_"+num2str(n)+"_"+num2str(stim_index)+"_"+num2str(e));
	    }
      }
    return ret;
  }   

  const vector<string> Subject_Model::get_names_d() const
  {
    Tracer_Plus trace("Subject_Model::get_names_d");
    
    vector<string> ret(model.get_num_d_params());
    
    int index=0;
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=model.marker_d[n-1].size(); i++)
	{
	  int node_index2=model.marker_d[n-1][i-1].second;  
	  int node_index=model.marker_d[n-1][i-1].first;
	  ret[index++]=string("d_"+num2str(n)+"_"+num2str(node_index)+"_"+num2str(node_index2));
	  //	  string("d_"+num2str(n)+"_"+num2str(model.marker_d[n-1][i-1]));
	}

    return ret;
  }

  const vector<string> Subject_Model::get_names_d(int n) const
  {
    Tracer_Plus trace("Subject_Model::get_names_d");
    
    vector<string> ret(get_num_d_node_params(n));
    
    int index=0;
    for(unsigned int i=1; i<=model.marker_d[n-1].size(); i++)
      {
	int node_index2=model.marker_d[n-1][i-1].second;  
	  int node_index=model.marker_d[n-1][i-1].first;
	  ret[index++]=string("d_"+num2str(n)+"_"+num2str(node_index)+"_"+num2str(node_index2));
	  //	  string("d_"+num2str(n)+"_"+num2str(model.marker_d[n-1][i-1]));
      }
    
    return ret;
  }

  const vector<string> Subject_Model::get_names_logsigmaa() const
  {
    Tracer_Plus trace("Subject_Model::get_names_logsigmaa");
    
    //    vector<string> ret(model.get_num_sigmaa_params()+1);
    vector<string> ret(1);
    
    ret[0]=string("logsigmaa");

    return ret;
  } 

  void Subject_Model::evaluate_neuronal_activity() 
  {
    Tracer_Plus trace("Subject_Model::evaluate_neuronal_activity");

    // Passed in a,b,c,d values are dcm style, i.e. dz/dt = Az+Bzu+Cu+Dzz, (apart from diag(a) and sigmaa) 
    // - need to convert them to be useful with the first order approximation to the diff eqn, i.e. dz/dt \approx (z_t-z_{t-1})/res:
    // diag(A) = (-1)*res*sigmaa + 1; 
    // otherwise A=A*res*sigmaa;
    // B=B*res*sigmaa;
    // C=C*res;
    // D=D*res*sigmaa;

    z.resize(nnodes);
    for(int n=1; n<=nnodes; n++)
      z[n-1].ReSize(ntpts);

    vector<float> S_tmo(nnodes,0);
    vector<float> S_tmt(nnodes,0);
 
    vector<float> z_t(nnodes,0);
    vector<float> z_tmo(nnodes,0);
    vector<float> z_tmt(nnodes,0);
    vector<float> z_error(nnodes,0);

    ////////////
    // adjust matrices for sigmaa
    vector<vector<float> > sigmaa_value_a=value_a;
    vector<vector<vector<float> > > sigmaa_value_b_amp_mod=value_b_amp_mod;
    vector<vector<vector<float> > > value_c_amp_mod_new=value_c_amp_mod;
    vector<vector<float> > sigmaa_value_d=value_d;

//     if(debuglevel==5)
//       {
// 	LOGOUT("Subject_Model::evaluate_neuronal_activity");
// 	for(int n=1; n<=nnodes; n++)
// 	  {
// 	    for(unsigned int i=1; i<=value_a[n-1].size(); i++)
// 	      LOGOUT(value_sigmaa);
	
// 	    for(unsigned int i=1; i<=value_b_amp_mod[n-1].size(); i++)		
// 	      for(unsigned int e=1; e<=value_b_amp_mod[n-1][i-1].size(); e++)
// 		LOGOUT(value_b[n-1][i-1]);
	
// 	    for(unsigned int i=1; i<=value_c_amp_mod[n-1].size(); i++)
// 	      for(unsigned int e=1; e<=value_c_amp_mod[n-1][i-1].size(); e++)
// 		LOGOUT(value_c[n-1][i-1]);
	    
//  	    for(unsigned int i=1; i<=value_d[n-1].size(); i++)	  
//  	      LOGOUT(value_sigmaa); 
// 	  }
//  	LOGOUT("Subject_Model::evaluate_neuronal_activity END");
//      }

    float value_sigmaa=std::exp(value_logsigmaa);

    for(int n=1; n<=nnodes; n++)
      {
	for(unsigned int i=1; i<=value_a[n-1].size(); i++)
	  sigmaa_value_a[n-1][i-1]*=value_sigmaa;
	
	for(unsigned int i=1; i<=value_b_amp_mod[n-1].size(); i++)		
	  for(unsigned int e=1; e<=value_b_amp_mod[n-1][i-1].size(); e++)
	    sigmaa_value_b_amp_mod[n-1][i-1][e-1]=(value_b[n-1][i-1]+value_b_amp_mod[n-1][i-1][e-1])*value_sigmaa;
	
	for(unsigned int i=1; i<=value_c_amp_mod[n-1].size(); i++)
	  for(unsigned int e=1; e<=value_c_amp_mod[n-1][i-1].size(); e++)
	    value_c_amp_mod_new[n-1][i-1][e-1]=(value_c[n-1][i-1]+value_c_amp_mod[n-1][i-1][e-1]);
	
	for(unsigned int i=1; i<=value_d[n-1].size(); i++)	  
	  sigmaa_value_d[n-1][i-1]*=value_sigmaa;	      
      }
    ////////////

    bool valid=true;  
    int jacob_count=0;
    Matrix J(nnodes,nnodes);
    Matrix Jinv_etc(nnodes,nnodes);
    float S=0;

    {
    Tracer_Plus trace("Subject_Model::evaluate_neuronal_activity_tloop");

    for(int t=1; t<=ntpts; t++)
      {	
	for(int n=0; n<nnodes; n++)
	  {
	    S=0;

	    // do diagonal of a
	    S+=-value_sigmaa*z_tmo[n];

	    for(unsigned int i=0; i<value_a[n].size(); i++)
	      {	
		S+=sigmaa_value_a[n][i]*z_tmo[model.marker_a[n][i]-1];
	      }
      
	    for(unsigned int i=0; i<value_b_amp_mod[n].size(); i++)
	      {	
		int stim_index=model.marker_b[n][i].second;
		for(unsigned int e=0; e<stimuli[stim_index-1].size(); e++)
		  {
		    S+=sigmaa_value_b_amp_mod[n][i][e]*z_tmo[model.marker_b[n][i].first-1]*stimuli[stim_index-1][e](t);
		  }
	      }
	    
	    for(unsigned int i=0; i<value_c_amp_mod_new[n].size(); i++)
	      {
		int stim_index=model.marker_c[n][i];
		for(unsigned int e=0; e<stimuli[stim_index-1].size(); e++)
		  {
		    S+=value_c_amp_mod_new[n][i][e]*stimuli[stim_index-1][e](t);
		  }
	      }
	    
	    for(unsigned int i=0; i<value_d[n].size(); i++)
	      {
		S+=sigmaa_value_d[n][i]*z_tmo[model.marker_d[n][i].first-1]*z_tmo[model.marker_d[n][i].second-1];
	      }
	    	   
	    z_t[n]=z_tmo[n]+model.res*S;

	    // calc error
	    float z_t2;
	    z_t2=z_tmt[n]+model.res*2*S_tmt[n];
	    z_error[n]=abs(z_t[n]-z_t2);	

	    S_tmo[n]=S;
	  }

	bool err=false;
	for(int n=0; n<nnodes; n++)
	  {
	    if(z_error[n]>0.01) err=true;
	  }

	if(err && t>2)	
	  {
	    jacob_count++;
	    
	    J=0;

	    //use jacobian
	    for(int n=0; n<nnodes; n++)  
	      {
		// n=n2
		J(n+1,n+1)+=-value_sigmaa;
		for(unsigned int i=0; i<value_a[n].size(); i++)
		  {	
		    J(n+1,model.marker_a[n][i])+=sigmaa_value_a[n][i];
		  }
      
		for(unsigned int i=0; i<value_b_amp_mod[n].size(); i++)
		  {	
		    int stim_index=model.marker_b[n][i].second;
		    for(unsigned int e=0; e<stimuli[stim_index-1].size(); e++)
		      {
			J(n+1,model.marker_b[n][i].first)+=sigmaa_value_b_amp_mod[n][i][e]*stimuli[stim_index-1][e](t);
		      }
		  }	    
	    
		for(unsigned int i=0; i<value_d[n].size(); i++)
		  {
		    J(n+1,model.marker_d[n][i].second)+=sigmaa_value_d[n][i]*z_tmo[model.marker_d[n][i].first-1];
		    J(n+1,model.marker_d[n][i].first)+=sigmaa_value_d[n][i]*z_tmo[model.marker_d[n][i].second-1];
		  }

	      }	     

	    try
	      {		
		Jinv_etc=J.i()*(expm(J*model.res)-IdentityMatrix(nnodes));
	      }
	    catch(Exception& exc)
	      {
		valid=false;
		break;
	      }	    

	    ColumnVector S_tmo_col(nnodes);
	    for(int n=0; n<nnodes; n++)
	      S_tmo_col(n+1)=S_tmo[n];
	    
	    for(int n=0; n<nnodes; n++)   
	      z_t[n] = z_tmo[n]+(Jinv_etc.Row(n+1)*S_tmo_col).AsScalar();

	  } // end of use jacobian if clause

	S_tmt=S_tmo;

	z_tmt=z_tmo;
	z_tmo=z_t;

	for(int n=0; n<nnodes; n++)
	  {
	    z[n](t)=z_t[n]; 

	    if(isnan(z[n](t)))
	      {
		valid=false;
	      }	
	  }	

      } // end of t loop
    }

//     if(jacob_count>0)
//       LOGOUT(float(jacob_count)/float(ntpts)*100);

    if(!valid || debuglevel==5)
      {
	for(int n=1; n<=nnodes; n++)
	  {
	    OUT(vector2ColumnVector(value_a[n-1]));
	    OUT(vector2ColumnVector(value_b[n-1]));
	    OUT(vector2ColumnVector(value_c[n-1]));
	    OUT(vector2ColumnVector(value_d[n-1]));

	    write_ascii_matrix(z[n-1],LogSingleton::getInstance().appendDir("node_z"+num2str(n)+num2str(debugcount1)));

 	  }
	
	OUT("************");
	float ss=0.0;    
	for(int n=0; n<nnodes; n++)
	  {
	    ColumnVector ab(nnodes);
	    ab=0;
	    for(unsigned int i=0; i<model.marker_a[n].size(); i++)
	      {
		int node_index=model.marker_a[n][i];
		ab(node_index)+=value_a[n][i];
		OUT(node_index);		
		OUT(value_a[n][i]);
	      }
	    ss+=SumSquare(ab);
	  }
	float boundary=nnodes/float(nnodes-1);    
	
	OUT(boundary);
	OUT(ss);        
	OUT("************");

// 	LOGOUT(debugcount1);
// 	debugcount1++;
      }       

    if(!valid)
      {
	OUT("Warning: Invalid Neuronal Model");
	for(int n=1; n<=nnodes; n++)
	  {
	    z[n-1]=0;
	  }
      }

    if(haemodynamic_model=="halfcosine")
      {
	// calc zfft
	zfft_real.resize(nnodes);
	zfft_imag.resize(nnodes);
	z_mean.ReSize(nnodes);
	
	for(int n=1; n<=nnodes; n++)
	  {
	    halfcos_hrf.fft_stimulus(z[n-1],zfft_real[n-1],zfft_imag[n-1],z_mean(n));
	  }
      }

  }  

  void Subject_Model::evaluate_hrf() 
  {
    Tracer_Plus trace("Subject_Model::evaluate_hrf");
    if(haemodynamic_model=="halfcosine")
      {
	hrf_fft_real.resize(nnodes);
	hrf_fft_imag.resize(nnodes);
	for(int n=1; n<=nnodes; n++)
	  {
	    evaluate_hrf(n);
	  }
      }
  }
 
  void Subject_Model::evaluate_hrf(int n) 
  {
    Tracer_Plus trace("Subject_Model::evaluate_hrf");
    if(haemodynamic_model=="halfcosine")
      {    
	halfcos_hrf.fft_hrf(subject_nodes[n-1]->get_hrf_value()(1),subject_nodes[n-1]->get_hrf_value()(2),subject_nodes[n-1]->get_hrf_value()(3),subject_nodes[n-1]->get_hrf_value()(4),subject_nodes[n-1]->get_hrf_value()(5),subject_nodes[n-1]->get_hrf_value()(6),hrf_fft_real[n-1],hrf_fft_imag[n-1]);
      }
  }

  void Subject_Model::get_hrf(vector<ColumnVector>& hrf, ColumnVector& t_hrf, float max_secs) 
  {
    Tracer_Plus trace("Subject_Model::get_hrf");

    // return haemodynamic response function, i.e. response if an impulse of activity is present in a node

    if(haemodynamic_model=="halfcosine")
      {
	
	hrf.resize(nnodes);
	for(int n=1; n<=nnodes; n++)
	  {	
	    halfcos_hrf.hrf_ifft(hrf_fft_real[n-1],hrf_fft_imag[n-1],hrf[n-1]);    
	  }
	
	// zoom in on time window
	int ntpts_window=int(max_secs/model.res);
	
	for(int n=1; n<=nnodes; n++)
	  hrf[n-1]=hrf[n-1].Rows(1,ntpts_window);       

      }
    else
      {
	// backup actual values
	vector<vector<ColumnVector> > stimuli_old=stimuli; // nstim*nevents*ntpts
	vector<ColumnVector> z_old=z; // nstim*ntpts
		
	    // set stimuli to be impulses at time zero, and set b's to zero
	    for(unsigned int s=1; s<=stimuli.size(); s++)
	      for(unsigned int e=1; e<=stimuli[s-1].size(); e++)
		{
		  stimuli[s-1][e-1]=0;
		  stimuli[s-1][e-1](1)=1;
		  //stimuli[s-1][e-1].Rows(1,10)=1;
		}
	    evaluate_neuronal_activity();

// 	for(unsigned int s=1; s<=z.size(); s++)
// 	  {
// 	    z[s-1]=0;
// 	    z[s-1].Rows(1,10)=1;
// 	  }
	
	evaluate_node_bold();
	
	hrf=node_bold_high_res;   
	
	// restore actual values
	stimuli=stimuli_old;
	z=z_old;
	
	evaluate_neuronal_activity();
	evaluate_node_bold();
	
	// zoom in on time window
	int ntpts_window=int(::round(max_secs/model.res));
	
	for(int n=1; n<=nnodes; n++)
	  hrf[n-1]=hrf[n-1].Rows(1,ntpts_window);
	
      }

    // zoom in on time window
    int ntpts_window=int(::round(max_secs/model.res));

    t_hrf.ReSize(ntpts_window);    
    float secs=model.res;
    for(int t=1; t<=ntpts_window; t++)
      {
	t_hrf(t)=secs;
	secs+=model.res;
      }
    
  }

  void Subject_Model::get_nrf(vector<vector<ColumnVector> >& nrf, ColumnVector& t_nrf, float max_secs) 
  {
    Tracer_Plus trace("Subject_Model::get_nrf");
    
    // return neural response function, i.e. response if an impulse of activity is fed through the network
    nrf.resize(stimuli.size());

    // backup actual values
    vector<vector<ColumnVector> > stimuli_old=stimuli; // nstim*nevents*ntpts
    vector<vector<float> > value_b_old=value_b; // nnodes*(inplay nnodes*nstim)

    int ntpts_window=int(::round(max_secs/model.res));
	
    // iterate through stimuli and for each one set the stimulus in question to be an impulse at time zero and all other stimuli to zero, and set b's to zero    
    for(unsigned int i=1; i<=value_b.size(); i++)
      for(unsigned int j=1; j<=value_b[i-1].size(); j++)
	value_b[i-1][j-1]=0;

    for(unsigned int s1=1; s1<=stimuli.size(); s1++)
      {
	
	for(unsigned int s2=1; s2<=stimuli.size(); s2++)
	  for(unsigned int e=1; e<=stimuli[s2-1].size(); e++)
	    {
	      stimuli[s2-1][e-1]=0;
	    }
      
	for(unsigned int e=1; e<=stimuli[s1-1].size(); e++)
	    {
	      stimuli[s1-1][e-1](1)=1;
	    }

	evaluate_neuronal_activity();
	nrf[s1-1]=z;    
	
	// restore actual values
	stimuli=stimuli_old;
	value_b=value_b_old;
	evaluate_neuronal_activity();
	
	// zoom in on time window	
	for(int n=1; n<=nnodes; n++)
	  nrf[s1-1][n-1]=nrf[s1-1][n-1].Rows(1,ntpts_window);
      }
    
    t_nrf.ReSize(ntpts_window);    
    float secs=model.res;
    for(int t=1; t<=ntpts_window; t++)
      {
	t_nrf(t)=secs;
	secs+=model.res;
      }

  }

  void Subject_Model::evaluate_node_bold() 
  {
    Tracer_Plus trace("Subject_Model::evaluate_node_bold");

    for(int n=1; n<=nnodes; n++)
      {
	evaluate_node_bold(n);
      }
  } 

   void Subject_Model::evaluate_node_bold(int n) 
  {
    Tracer_Plus trace("Subject_Model::evaluate_node_bold(n)");

    if(haemodynamic_model=="halfcosine")
      evaluate_node_bold_halfcosine(n);
    else
      {
	evaluate_node_cbf_balloon(n);
	evaluate_node_bold_balloon(n);	
      }
  }

  void Subject_Model::evaluate_node_bold_halfcosine(int n) 
  {
    Tracer_Plus trace("Subject_Model::evaluate_node_bold(n)");

    halfcos_hrf.convolve_hrf_fft(zfft_real[n-1],zfft_imag[n-1],z_mean(n),hrf_fft_real[n-1],hrf_fft_imag[n-1],node_bold[n-1]);

    // regress out confound EVs
    node_bold[n-1] = residual_forming_confound_evs*node_bold[n-1];
  }
  
  void Subject_Model::evaluate_node_cbf_balloon(int n) 
  {
    Tracer_Plus trace("Subject_Model::evaluate_node_cbf_balloon(n)");

    float kappa=subject_nodes[n-1]->get_kappa_value();
    float gamma=subject_nodes[n-1]->get_gamma_value();

//     if(debuglevel==5)
//       {
// 	LOGOUT(kappa);
// 	LOGOUT(gamma);
//       }

    float s_tmo=0; // s_{t-1}
    float f_tmo=1; // f_{t-1}
 
    float s_t=s_tmo;
    float log_f_tmo=log(f_tmo); 
    float log_f_t=log_f_tmo; 
    float s_tmt=s_tmo;
    float log_f_tmt=log(f_tmo);

    float fs_tmo=0;
    float fs_tmt=0;
    float ff_tmo=0;
    float ff_tmt=0;

    ColumnVector s_ts(ntpts);
    ColumnVector f_ts(ntpts);
    ColumnVector stims(ntpts);

    bool high_f=false;

    // Jacobian
    // syms z f log_f s kappa gam real    
    // simplify(diff(z-kappa*s-gam*(exp(log_f)-1),s))
    // simplify(diff(z-kappa*s-gam*(exp(log_f)-1),log_f))
    // simplify(diff(s/exp(log_f),log_f))
    // simplify(diff(s/exp(log_f),s))
    Matrix J(2,2); 
    J=0;
    //ds_ds
    J(1,1) = -kappa;

    Matrix Jinv_etc(2,2);
    ColumnVector fs_ff(2);
    DiagonalMatrix eye2=IdentityMatrix(2);
    int jacob_count=0;
    bool valid=true;

    for(int t=1; t<=ntpts; t++)
      {	
	float stim=z[n-1](t);

	fs_tmo=stim-kappa*s_tmo-gamma*(f_tmo-1);
	ff_tmo=s_tmo/f_tmo;

	s_t=s_tmo+model.res*fs_tmo;
	log_f_t=log_f_tmo+model.res*ff_tmo;      

	// test error
	float s_t2=s_tmt+model.res*2*fs_tmt;
	float s_error=abs(s_t-s_t2);
	float log_f_t2=log_f_tmt+model.res*2*ff_tmt;
	float log_f_error=abs(log_f_t-log_f_t2);

	if((s_error>0.01 || log_f_error>0.01) && t>1 && !isnan(stim))	
	  {
	    jacob_count++;

	    //use jacobian

	    //ds_ds
	    //J(1,1) = -kappa;
	    
	    //ds_df 
	    J(1,2) = -gamma*f_tmo;
	    
	    //df_df
	    J(2,2) = -s_tmo/f_tmo;
	    
	    //df_ds 
	    J(2,1) = 1/f_tmo;	    

// 	    try
// 	      {
// 		Jinv_etc=J.i();
// 		//Jinv_etc=J.i()*(expm(J*model.res)-eye2);
// 	      }
// 	    catch(Exception& exc)
// 	      {
// 		OUT("CBF inv(J) failed");
// 		OUT(J);
// 		OUT(stim);
// 		OUT(kappa);
// 		OUT(gamma);
// 		OUT(f_tmo);
// 		OUT(s_tmo);

// 		exit(1);
// 		valid=false;
// 		break;
// 	      }
	    
// 	    try
// 	      {
// 		Jinv_etc=expm(J*model.res);
// 		//Jinv_etc=J.i()*(expm(J*model.res)-eye2);
// 	      }
// 	    catch(Exception& exc)
// 	      {
// 		OUT("CBF expm(J) failed");
// 		OUT(J*model.res);
// 		OUT(stim);
// 		OUT(kappa);
// 		OUT(gamma);
// 		OUT(f_tmo);
// 		OUT(s_tmo);

// 		exit(1);
// 		valid=false;
// 		break;
// 	      }
	    
	    try
	      {
		Jinv_etc=J.i()*(expm(J*model.res)-eye2);
	      }
	    catch(Exception& exp)
	      {
// 		OUT(J);
		valid=false;
		break;
	      }

	    //fs
	    fs_ff(1)=fs_tmo;

	    //ff
	    fs_ff(2)=ff_tmo;
	    
 	    s_t = s_tmo+(Jinv_etc.Row(1)*fs_ff).AsScalar();
	    log_f_t = log_f_tmo+(Jinv_etc.Row(2)*fs_ff).AsScalar();
	  }

	log_f_tmt=log_f_tmo;
	s_tmt=s_tmo;

	s_tmo=s_t;
	f_tmo=exp(log_f_t);

	if(f_tmo>8) 
	  {
	    f_tmo=8; 
	    high_f=true;
	  }
	if(f_tmo<1.0/8.0) 
	  {
	    f_tmo=1.0/8.0; 
	    high_f=true;
	  }

	log_f_tmo=log_f_t;

	fs_tmt=fs_tmo;
	ff_tmt=ff_tmo;

	stims(t)=stim;
	s_ts(t)=s_tmo;
	node_cbf[n-1](t)=f_tmo;
      }

//     if(jacob_count>0)
//       LOGOUT(float(jacob_count)/float(ntpts)*100);
  
    for(int t = 1; t <= node_cbf[n-1].Nrows(); t++)
      {
	if(isnan(node_cbf[n-1](t)))
	  {
	    valid=false;	    
	  }
	 
// 	OUT(t);
      }

//     if(!valid || debuglevel==5)
//       {
//  	OUT(n);
//  	OUT(gamma);OUT(kappa);
// 	write_ascii_matrix(stims,LogSingleton::getInstance().appendDir("stims"+num2str(n)));   
// 	write_ascii_matrix(s_ts,LogSingleton::getInstance().appendDir("s_t"+num2str(n)));
// 	write_ascii_matrix(node_cbf[n-1],LogSingleton::getInstance().appendDir("node_cbf"+num2str(n)));
//       }

    if(!valid)
      {
	OUT("Invalid CBF Balloon Model");      
	node_cbf[n-1]=0;
      }

    // load s_t1; load stims1; load node_cbf1; figure;plot(node_cbf1);
    // load s_t2; load stims2; load node_cbf2; figure;plot(node_cbf2);
    // load s_t3; load stims3; load node_cbf3; figure;plot(node_cbf3);

    // kappa = 0.65;gamma = 0.41;tau=0.98;alpha=0.32;E0=0.34
    // gamma=0.346773;kappa=0.502268;E0=0.400839;alpha=0.286174;tau=0.592265
    // st=stims2; N=length(st);res=0.2; [s,f,v,q,b]=my_solve2(st,alpha,tau,E0,kappa,gamma,res);t=(res:res:res*N);plot(t,b);

     if(high_f)
       LOGOUT("WARNING - implied % CBF is greater than 800%");
  }

  void Subject_Model::evaluate_node_bold_balloon(int n) 
  {
    Tracer_Plus trace("Subject_Model::evaluate_node_bold_balloon(n)");

    float tau=subject_nodes[n-1]->get_tau_value();
    float invalpha=1.0/subject_nodes[n-1]->get_alpha_value();
    float E0=subject_nodes[n-1]->get_E0_value();
    float epsilon=std::exp(subject_nodes[n-1]->get_logepsilon_value());
 
    float V0;
    float k1;
    float k2;
    float k3;
 
    if(haemodynamic_model=="balloon")
      { 	      
	//from DCM paper:
	V0=100*0.04;
	k1=7*E0;
	k2=2;
	k3=2*E0-0.2;
      }
    else
      {
	// see stephan NI 2007 and spm_gx_dcm.m
	V0=100*0.04;
	float TE=0.040;
	float r0=25;
	float nu0=40.3;
	
	k1=4.3*nu0*E0*TE;
	k2=epsilon*r0*E0*TE;
	k3=1-epsilon;
      }

    float v_tmo=1;
    float q_tmo=1;

    float log_v_t=log(v_tmo);
    float log_q_t=log(v_tmo);
    float log_v_tmt=log(v_tmo);
    float log_q_tmt=log(q_tmo);

    float log_v_tmo=log_v_t;
    float log_q_tmo=log_q_t;

    float fv_tmo=0;
    float fv_tmt=0;
    float fq_tmo=0;
    float fq_tmt=0;

    ColumnVector v_ts(ntpts);
    ColumnVector q_ts(ntpts);

    // Jacobian
    // syms E0fterm log_v f  log_q tau invalpha E0 real
    // simplify(diff((f-exp(log_v)^(invalpha))/(tau*exp(log_v)),log_v))
    // simplify(diff((f-exp(log_v)^(invalpha))/(tau*exp(log_v)),log_q))
    // simplify(diff((E0fterm-exp(log_v)^(invalpha)*exp(log_q)/exp(log_v))/(tau*exp(log_q)),log_q))
    // simplify(diff((E0fterm-exp(log_v)^(invalpha)*exp(log_q)/exp(log_v))/(tau*exp(log_q)),log_v))
    Matrix J(2,2);
    J=0;
    Matrix Jinv_etc(2,2);
    ColumnVector fv_fq(2);
    DiagonalMatrix eye2=IdentityMatrix(2);
    int jacob_count=0;
    bool valid=true;  

    for(int t=1; t<=ntpts; t++)
      {	
	float f=node_cbf[n-1](t);
	float E0fterm=f*(1-pow(1-E0,1/f))/E0; 
	float v_to_invalpha=pow(v_tmo,invalpha);

	fv_tmo=(f-v_to_invalpha)/(tau*v_tmo);
	fq_tmo=(E0fterm - v_to_invalpha*q_tmo/v_tmo)/(tau*q_tmo);
	log_v_t=log_v_tmo+model.res*fv_tmo;
	log_q_t=log_q_tmo+model.res*fq_tmo;

	// test error
	float log_v_t2=log_v_tmt+model.res*2*fv_tmt;
	float log_v_error=abs(log_v_t-log_v_t2);
	float log_q_t2=log_q_tmt+model.res*2*fq_tmt;
	float log_q_error=abs(log_q_t-log_q_t2);

//  	LOGOUT(log_v_error);
//  	LOGOUT(log_q_error);

	if((log_v_error>0.01 || log_q_error>0.01) && t>1 && !isnan(f))
	  {	
	    jacob_count++;

	    //use jacobian

	    //dv_dv 
	    J(1,1) = -( exp(log_v_tmo*(invalpha - 1)) * (invalpha - 1) + f/v_tmo )/tau;
	    
	    //dv_dq 
	    //J(1,2) = 0;
	    
	    //dq_dq 
	    J(2,2) = -E0fterm/(q_tmo*tau);
	    
	    //dq_dv 
	    J(2,1) = -exp(log_v_tmo*(invalpha-1))*(invalpha-1)/tau;
       
// 	    try
// 	      {
// 		Jinv_etc=J.i();
// 	      }
// 	    catch(Exception& exc)
// 	      {		
// 		OUT("BOLD J.i() failed");
// 		OUT(J);
// // 		OUT(tau);
// // 		OUT(invalpha);
// // 		OUT(f);
// // 		OUT(E0);
// 		OUT(q_tmo);
// 		OUT(v_tmo);
// 		valid=false;
// 		break;
// 	      }
	    
// 	    try
// 	      {
// 		Jinv_etc=expm(J*model.res);
// 	      }
// 	    catch(Exception& exc)
// 	      {		
// 		OUT("BOLD expm(J) failed");
// 		OUT(J*model.res);
// // 		OUT(tau);
// // 		OUT(invalpha);
// // 		OUT(f);
// // 		OUT(E0);
// 		OUT(q_tmo);
// 		OUT(v_tmo);
// 		valid=false;
// 		break;
// 	      }	    

	     try
	      {
		Jinv_etc=J.i()*(expm(J*model.res)-eye2);
	      }
	    catch(Exception& exp)
	      {
	 	OUT(J);
		OUT(f);
		OUT(q_tmo);
		OUT(v_tmo);
		OUT(E0);
		valid=false;
		break;
	      }
	    
	    //fv
	    fv_fq(1)=fv_tmo;
	    //fq
	    fv_fq(2)=fq_tmo;
  
 	    log_v_t= log_v_tmo+(Jinv_etc.Row(1)*fv_fq).AsScalar();
 	    log_q_t = log_q_tmo+(Jinv_etc.Row(2)*fv_fq).AsScalar();
	  }	  

	log_v_tmt=log_v_tmo;
	log_q_tmt=log_q_tmo;

	log_v_tmo=log_v_t;
	log_q_tmo=log_q_t;

	v_tmo=exp(log_v_tmo);
	q_tmo=exp(log_q_tmo);

	fv_tmt=fv_tmo;
	fq_tmt=fq_tmo;

	v_ts(t)=v_tmo;
	q_ts(t)=q_tmo;

	// multiply by 100 to get as % signal change
	node_bold_high_res[n-1](t)=V0*(k1*(1-q_tmo)+k2*(1-q_tmo/v_tmo)+k3*(1-v_tmo));

      }

//     if(jacob_count>0)
//       LOGOUT(float(jacob_count)/float(ntpts)*100);

    // sub sample  
    for(int scan = 1; scan <= node_bold[n-1].Nrows(); scan++)
      {

	node_bold[n-1](scan) = node_bold_high_res[n-1](int(::round(scan*model.tr/model.res)));
	if(isnan(node_bold[n-1](scan)))
	  {
	    valid=false;
	  }
	 
// 	OUT(scan);
// 	OUT(int(::round(scan*model.tr/model.res)));
// 	OUT(node_bold_high_res[n-1](int(::round(scan*model.tr/model.res))));
      }

    if(!valid || debuglevel==10)
      {
	OUT(n);
	float alpha=1.0/invalpha;
	OUT(E0);OUT(alpha);OUT(tau);OUT(epsilon);
// 	write_ascii_matrix(node_bold_high_res[n-1],LogSingleton::getInstance().appendDir("node_bold_hr"+num2str(n)));       
// 	write_ascii_matrix(v_ts,LogSingleton::getInstance().appendDir("v_t"+num2str(n)));
// 	write_ascii_matrix(q_ts,LogSingleton::getInstance().appendDir("q_t"+num2str(n)));	
      }

    if(!valid)
      {
	OUT("Invalid BOLD Balloon Model");
	node_bold[n-1]=0;
      }

    // load v_t1; load q_t1; load stims1; load node_bold_hr1; figure;plot(node_bold_hr1);
    ///load v_t2; load q_t2; load stims2; load node_bold_hr2; figure;plot(node_bold_hr2);
    // load v_t3; load q_t3; load stims3; load node_bold_hr3; figure;plot(node_bold_hr3);
    // kappa = 0.65;gamma = 0.41;tau=0.98;alpha=0.32;E0=0.34
    // gamma=0.346773;kappa=0.502268;E0=0.400839;alpha=0.286174;tau=0.592265
    // st=stims2; N=length(st);res=0.2; [s,f,v,q,b]=my_solve2(st,alpha,tau,E0,kappa,gamma,res);t=(res:res:res*N);plot(t,b);

    // regress out confound EVs
    node_bold[n-1] = residual_forming_confound_evs*node_bold[n-1];

    //     if(debuglevel==5)
    //       {
    // 	write_ascii_matrix(node_bold[n-1],LogSingleton::getInstance().appendDir("node_bold"+num2str(n)+num2str(debugcount2)));
    
    // 	LOGOUT(debugcount2);
    // 	debugcount2++;
    
    //       }

  }
  
  void Subject_Model::evaluate_decoded_node_data() 
  {
    Tracer_Plus trace("Subject_Model::evaluate_decoded_node_data");

    // in this case actually evaluate node_data=voxelwise_data*pvf  (nscans*nvox)*(nvox*1) for each node
    
    decoded_node_data.resize(nnodes);
    
    //const Matrix& voxelwise_data=subject.get_voxelwise_data();
    for(int n=1; n<=nnodes; n++)
      {
	evaluate_decoded_node_data(n);
      }
  }

  void Subject_Model::evaluate_decoded_node_data(int n) 
  {
    Tracer_Plus trace("Subject_Model::evaluate_decoded_node_data(n)");
    
    voxelwise_pvf[n-1]=0;
    decoded_node_data[n-1].ReSize(voxelwise_data.Ncols());
    decoded_node_data[n-1]=0;

    float sum_pvf=0;
    float max_pvf=0;
    int max_i=1;

    for(unsigned int i =1; i<=get_voxel_coordinates().size(); i++) // loop through voxels
      {
	if(is_voxel_in_roi(i,n)) // need to include contributions to this voxel from roi
	  {
	    RowVector vox_data=voxelwise_data.Row(i);
	    float pvf=get_partial_vol_fraction(i,n);
	    voxelwise_pvf[n-1](i)=pvf;
	    decoded_node_data[n-1]+=pvf*vox_data.t();
	    sum_pvf+=pvf;

	    if(pvf>max_pvf) 
	      {
		max_pvf=pvf;
		max_i=i;
	      }
	  }
      }
    
    decoded_node_data[n-1]/=sum_pvf;

    float tmp=stdev(node_data[n-1]).AsScalar();
    float tmp2=stdev(decoded_node_data[n-1]).AsScalar();
    decoded_node_data[n-1]=tmp*decoded_node_data[n-1]/tmp2;

    //     RowVector vox_data=voxelwise_data.Row(max_i);
    //     float tmp=stdev(vox_data.t()).AsScalar();
    //     float tmp2=stdev(decoded_node_data[n-1]).AsScalar();
    //     decoded_node_data[n-1]=tmp*decoded_node_data[n-1]/tmp2;
  }
  
  void Subject_Model::evaluate_voxelwise_bold() 
  {
    Tracer_Plus trace("Subject_Model::evaluate_voxelwise_bold()");
    
	// Matrix& voxelwise_bold is num_voxels_in_allnodes*nscans
	vector<RowVector> node_bold_tmp(nnodes);
	for(int n=1; n<=nnodes; n++)
	  {
	    voxelwise_pvf[n-1]=0;
	    node_bold_tmp[n-1]=node_bold[n-1].t();
	  }
	
	//    voxelwise_bold.ReSize(get_voxel_coordinates().size(),halfcos_hrf.get_nscans());	  
	voxelwise_bold=0;
	
	for(unsigned int i =1; i<=get_voxel_coordinates().size(); i++)
	  {
	    for(int n2=1; n2<=nnodes; n2++)
	      {	
		if(is_voxel_in_roi(i,n2)) // need to include contributions to this voxel from all rois
		  {
		    float pvf=get_partial_vol_fraction(i,n2);
		    voxelwise_pvf[n2-1](i)=pvf;
		    voxelwise_bold.Row(i)+=pvf*node_bold_tmp[n2-1];
		  }
	      }	  
	  }
      
  }
 
  void Subject_Model::evaluate_voxelwise_bold(int n) 
  {
    Tracer_Plus trace("Subject_Model::evaluate_voxelwise_bold(n)");

	// Matrix& voxelwise_bold is num_voxels_in_allnodes*nscans
	
	vector<RowVector> node_bold_tmp(nnodes);
	for(int n2=1; n2<=nnodes; n2++)
	  {
	    node_bold_tmp[n2-1]=node_bold[n2-1].t();
	  }
	
	for(unsigned int i=1; i<=get_voxel_coordinates().size(); i++)
	  {
	    if(is_voxel_in_roi(i,n)) 	   // only need to bother with voxel if it is in this ROI
	      {
		voxelwise_bold.Row(i)=0;
		for(int n2=1; n2<=nnodes; n2++)
		  {	
		    if(is_voxel_in_roi(i,n2)) // need to include contributions to this voxel from all rois
		      {
			float pvf=get_partial_vol_fraction(i,n2);
			voxelwise_pvf[n2-1](i)=pvf;
			voxelwise_bold.Row(i)+=pvf*node_bold_tmp[n2-1];  
		      }
		  }
	      }
	  }
      
  }

  void Subject_Model::setup_voxelwise_data(const volume4D<float>& data, Matrix& pvoxelwise_data) 
  {
    Tracer_Plus trace("Subject_Model::setup_voxelwise_data");

    // note that voxelwise data is in functional space  
    voxelwise_data.ReSize(voxel_coordinates.size(), data.tsize());
    voxelwise_data=0;   

    for(unsigned int i=1; i<=voxel_coordinates.size(); i++)
      {
	ColumnVector coords=voxel_coordinates[i-1];
	ColumnVector tmp=data.voxelts(int(coords(1)),int(coords(2)),int(coords(3)));
	float mn=mean(tmp).AsScalar();
	float sd=sqrt(var(tmp).AsScalar());

	voxelwise_data.Row(i)=tmp.t();
	if(mn>sd)	    
	  {
	    voxelwise_data.Row(i)=100*voxelwise_data.Row(i)/mn;
	  }
	else if(haemodynamic_model!="halfcosine")
	  {
	    LOGOUT(mn);
	    LOGOUT(sd);
	    LOGOUT("Mean(data)<std(data). Data has a mean and standard deviation that is not consistent with BOLD FMRI data. NMA requires data to be inputted un-demeaned so that it can convert to approximate % signal change data.");
	    throw Exception("Mean(data)<std(data). Data has a mean and standard deviation that is not consistent with BOLD FMRI data. NMA requires data to be inputted un-demeaned so that it can convert to approximate % signal change data.");
	  }
      }
	
    pvoxelwise_data=voxelwise_data;
  }

  void Subject_Model::establish_voxel_coordinates(const volume<float>& mask) 
  {
    Tracer_Plus trace("Subject_Model::establish_voxel_coordinates");
        
    // mask in func space
    // coordinates are in functional space
    // need to establish coordinates of voxels within rois and which voxels belong to which rois
    ColumnVector coords(3);
    
    vector<ColumnVector> min_voxel_coords(nnodes);
    vector<ColumnVector> max_voxel_coords(nnodes);
    for(int r=0; r<nnodes; r++)
      {
	min_voxel_coords[r].ReSize(3);min_voxel_coords[r]=1e10;
	max_voxel_coords[r].ReSize(3);max_voxel_coords[r]=0;
      }
 
    //    save_volume(subject_nodes[1]->get_func_mask(), LogSingleton::getInstance().appendDir("VTA_my_func_mask"));	  
 
    for(int x=0; x<mask.xsize(); x++)
      for(int y=0; y<mask.ysize(); y++)
	for(int z=0; z<mask.zsize(); z++)
	  {
	    ColumnVector in_roi(nnodes);in_roi=0;
	    for(int r=0; r<nnodes; r++)
	      {
		if(subject_nodes[r]->get_func_mask()(x,y,z)>0 && mask(x,y,z)>0)
		  {
		    if(r==1)
		      {
			OUT(x);
			OUT(y);
			OUT(z);
		      }

		    in_roi(r+1)=1;
		  }
	      }

	    if(sum(in_roi).AsScalar()>0)
	      {
		coords(1)=x;coords(2)=y;coords(3)=z;		
		voxel_coordinates.push_back(coords);
		voxel_in_roi.push_back(in_roi);

		for(int r=0; r<nnodes; r++)
		  if(in_roi(r+1)>0)
		    {
		      if(coords(1)<min_voxel_coords[r](1)) min_voxel_coords[r](1)=coords(1);
		      if(coords(1)>max_voxel_coords[r](1)) max_voxel_coords[r](1)=coords(1);
		      if(coords(2)<min_voxel_coords[r](2)) min_voxel_coords[r](2)=coords(2);
		      if(coords(2)>max_voxel_coords[r](2)) max_voxel_coords[r](2)=coords(2);
		      if(coords(3)<min_voxel_coords[r](3)) min_voxel_coords[r](3)=coords(3);
		      if(coords(3)>max_voxel_coords[r](3)) max_voxel_coords[r](3)=coords(3);
		    }
	      }
	  }

    for(int r=0; r<nnodes; r++)
      {
// 	subject_nodes[r]->set_min_func_space_voxel_coordinate(min_voxel_coords[r]);
// 	subject_nodes[r]->set_max_func_space_voxel_coordinate(max_voxel_coords[r]);
	subject_nodes[r]->set_range_func_space_voxel_coordinate(min_voxel_coords[r],max_voxel_coords[r], model.single_timeseries);


	OUT("====================== Subject_Model::establish_voxel_coordinates");
	OUT(subject_nodes[r]->node.get_name());
	OUT(min_voxel_coords[r]);
	OUT(max_voxel_coords[r]);
	OUT("======================");
      }

  }
}
