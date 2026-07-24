/*  subject.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT  */

#include "subject.h"
#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"
#include "mcmc_mh.h"
#include "model.h"
#include "node.h"
#include "nma_manager.h"
#include <algorithm>

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace std;
using namespace MISCPLOT;
using namespace MISCPIC;

namespace Nma {  

  Subject::Subject(const string& pname, const string&  pdata_dir, const Model& pmodel, bool psingle_timeseries, int pdata_mode, string phaemodynamic_model, bool prandom_initialise,int pdebuglevel, bool pdecode, bool pphi_every_voxel) :
    single_timeseries(psingle_timeseries),
    haemodynamic_model(phaemodynamic_model),
    decode(pdecode),
    phi_every_voxel(pphi_every_voxel),
    data_mode(pdata_mode),
    name(pname),
    data_dir(pdata_dir),
    subject_model(pname,pmodel,phaemodynamic_model,pdebuglevel),
    already_setup_vb_data(false),
    random_initialise(prandom_initialise),
    debuglevel(pdebuglevel)
  {
    setup();
  }

  void Subject::setup()
  {
    Tracer_Plus trace("Subject::setup");

    nnodes=subject_model.get_subject_nodes().size();
    node_data.resize(nnodes);

    // load data and set nsecs
    if(data_mode==0)
      {
	for(int r=0; r<nnodes; r++)
	  {	    
	    node_data[r]=read_ascii_matrix(string(data_dir+"/"+name+"/"+subject_model.get_model().get_node_names()[r]+"_data.txt"));

	    subject_model.set_nsecs(node_data[r].Nrows()*subject_model.model.tr);
	    LOGOUT(node_data[0].Nrows());	    
	  }
      }
    else  
      {
	// load data
	read_volume4D(data,data_dir+string("/")+name+"/data");

	// set func dim
	funcvoxdim.ReSize(3);
	funcvoxdim(1)=data.xdim();
	funcvoxdim(2)=data.ydim();
	funcvoxdim(3)=data.zdim();

 	LOGOUT(funcvoxdim);

	subject_model.set_nsecs(data.tsize()*subject_model.model.tr);
	LOGOUT(data.tsize());
      }


    LOGOUT(subject_model.nsecs);
    LOGOUT(subject_model.model.tr);

    subject_model.setup(); // now that we have set nsecs we can call setup on the model

    // scale data and regress out confound EVs
    if(data_mode==0)
      {	

	for(int r=0; r<nnodes; r++)
	  {
	    // approximate data as % signal change by dividing by the mean (but only if data has not been demeaned
	    // otherwise do nothing and assume that values are percent signal change scaled)
	    float mn=mean(node_data[r]).AsScalar();
	    float sd=sqrt(var(node_data[r]).AsScalar());

//  	    OUT(mn);
//  	    OUT(sd);
	    if(mn>sd)
	      node_data[r]=100*node_data[r]/mn;
 	    else
 	      node_data[r]=node_data[r];

	    node_data[r] = subject_model.get_residual_forming_confound_evs()*node_data[r];
	  }
      }
    else
      {
	// load subject mask
	read_volume(func_mask,data_dir+string("/")+name+"/mask");  
	if((func_mask.xdim() != funcvoxdim(1)) ||(func_mask.ydim() != funcvoxdim(2)) || (func_mask.zdim() != funcvoxdim(3)))
	  {
	    LOGOUT("func voxel dimensions != func data voxel dimensions");
	    throw Exception(string("func voxel dimensions != func data voxel dimensions").data());		    
	  }

	// load xform stuff for this subject
	std2func_xform=read_ascii_matrix(string(data_dir+"/"+name+"/standard2example_func.mat"));
	if(std2func_xform.Nrows()!=4 || std2func_xform.Ncols()!=4)
	  {
	    LOGOUT(size(std2func_xform));
	    LogSingleton::getInstance().str() << "std2func_xform needs to be a 4x4 matrix" << endl;
	    throw Exception(string("std2func_xform needs to be a 4x4 matrix").data());
	  }

	// load standard brain
	read_volume(standard_brain,data_dir+string("/standard_brain"));  

	// set std dim
	stdvoxdim.ReSize(3);
	stdvoxdim(1)=standard_brain.xdim();
	stdvoxdim(2)=standard_brain.ydim();
	stdvoxdim(3)=standard_brain.zdim();

 	LOGOUT(stdvoxdim);

	for(int r=0; r<nnodes; r++)
	  {	    
	    // setup xforms between std space and func space
	    subject_model.get_subject_nodes()[r]->set_std2func_info(std2func_xform, stdvoxdim, funcvoxdim, func_mask);	
	  }

	// establish voxel coordinates and which voxels belong to which roi
	subject_model.establish_voxel_coordinates(func_mask);

	// setup voxelwise data
	subject_model.setup_voxelwise_data(data,voxelwise_data);
	
	// regress out confound EVs and compute average time series across node
	for(int r=0; r<nnodes; r++)
	  {
	    node_data[r].ReSize(data.tsize());
	    node_data[r] = 0; // will set node_data as average voxelwise_data
	  }

	ColumnVector count(nnodes);count=0;	
	vector<Matrix> voxdata(nnodes);

	for(int r=0; r<nnodes; r++)
	  {
	    voxdata[r].ReSize(voxelwise_data.Nrows(),voxelwise_data.Ncols());
	  }
	
	for(int v=1; v<=voxelwise_data.Nrows(); v++)
	  {
	    ColumnVector tmp = subject_model.get_residual_forming_confound_evs()*voxelwise_data.Row(v).t();
	    
// 	    float stdev=std::sqrt(var(tmp).AsScalar());
// 	    tmp/=stdev;

	    voxelwise_data.Row(v) = tmp.t();
	    for(int r=0; r<nnodes; r++)
	      if(subject_model.is_voxel_in_roi(v,r+1))
		{
		  node_data[r]+=tmp;
		  count(r+1)++;
		  voxdata[r].Row(int(count(r+1)))=tmp.t();
		}	    
	  }
	    
	for(int r=0; r<nnodes; r++)
	  {
	    voxdata[r]=voxdata[r].Rows(1,int(count(r+1)));
	  }
	
	// set node_data as average voxelwise_data
	for(int r=0; r<nnodes; r++)
	  node_data[r]/=count(r+1);      	

// 	////////////////////////////////////
// 	// set node data as principal component       
// 	DiagonalMatrix eigenvals;
// 	Matrix eigenvecs;  
// 	for(int r=0; r<nnodes; r++)
// 	  {
// // 	    write_ascii_matrix(voxdata[r],LogSingleton::getInstance().appendDir("voxdata_"+num2str(r)));

// 	    try
// 	      {
// 		SymmetricMatrix corr;
// 		corr << voxdata[r].t()*voxdata[r];
		
// 		EigenValues(corr, eigenvals, eigenvecs);
// 	      }
// 	    catch(Exception& e) 
// 	      {
// 		cerr << endl << e.what() << endl;
		
// 		throw e;
// 	      }	
	    
// // 	    write_ascii_matrix(eigenvecs, LogSingleton::getInstance().appendDir("eigenvecs.txt"));
// // 	    write_ascii_matrix(diag(eigenvals), LogSingleton::getInstance().appendDir("eigenvals.txt"));
// // 	    OUT(size(eigenvals));
// // 	    OUT(size(eigenvecs));

// 	    // need to flip cos of the way EigenValues does its output
// 	    //eigenvecs = eigenvecs.Reverse();
// 	    //eigenvals = eigenvals.Reverse();
// 	    //node_data[r] = eigenvecs.Column(1).Reverse();
// 	    node_data[r] = eigenvecs.Column(eigenvecs.Ncols());
	 
// 	  }

//       }
//     ////////////////////////////////////

	subject_model.set_voxelwise_data(voxelwise_data);
	subject_model.set_node_data(node_data);
      }

    // initialise
    subject_model.initialise(random_initialise);
    
    //    if(!use_balloon_model)

    if(random_initialise)
      fit_model_using_c_matrix_only(*this,0);
  }
   
  void calculate_noise_precision(const Subject& subject, ColumnVector& noise_precision)
  {
    Tracer_Plus trace("calculate_noise_precision");
    
    const Subject_Model& subject_model=subject.get_subject_model();
   
    if(subject.is_single_timeseries())    
      {
	const vector<ColumnVector>& node_bold = subject.get_subject_model().get_node_bold();
	const vector<ColumnVector>& node_data = subject.get_node_data();

	for(unsigned int r=0; r<node_bold.size(); r++)
	  {
	    double ss = SumSquare(node_bold[r]-node_data[r]);
	    double est_node_phi=(node_data[r].Nrows())/ss;

	    ColumnVector ones=node_bold[r];ones=1;

	    if(r==0)
	      noise_precision=est_node_phi*ones;
	    else
	      noise_precision&=est_node_phi*ones;	    
	  }
      }   
    else if(subject.is_decode())
      {
	const vector<ColumnVector>& node_bold = subject.get_subject_model().get_node_bold();
	const vector<ColumnVector>& node_data = subject.get_subject_model().get_decoded_node_data();

	for(unsigned int r=0; r<node_bold.size(); r++)
	  {
	    double ss = SumSquare(node_bold[r]-node_data[r]);
	    double est_node_phi=(node_data[r].Nrows())/ss;

	    ColumnVector ones=node_bold[r];ones=1;

	    if(r==0)
	      noise_precision=est_node_phi*ones;
	    else
	      noise_precision&=est_node_phi*ones;	    
	  }  
      }
    else
      {
	const Matrix& voxelwise_bold = subject.get_subject_model().get_voxelwise_bold();
	const vector<ColumnVector>& voxelwise_pvf = subject_model.get_voxelwise_pvf();
	Matrix noiseprec_tmp=voxelwise_bold;

	for(int n=0; n<voxelwise_bold.Nrows(); n++)
	  {		
	    double variance;
	    if(subject.is_phi_every_voxel())
	      variance = 1.0/std::exp(subject_model.get_log_phi_voxel()(1)+subject_model.get_log_phi_every_voxel()(n+1));
	    else
	      variance = 1.0/std::exp(subject_model.get_log_phi_voxel()(1));
		

	    for(unsigned int r=0; r<voxelwise_pvf.size(); r++)
	      {	       	     
		variance += 1.0/subject_model.get_phi_node()[r]*Sqr(voxelwise_pvf[r](n+1));
	      }
	    
	    float est_node_phi=1.0/variance;

	    noiseprec_tmp.Row(n+1)=est_node_phi;

// 	    RowVector tmp=voxelwise_bold.Row(n+1);
// 	    ColumnVector ones=tmp.t();ones=1;
	    
// 	    if(n==0)
// 	      noise_precision=est_node_phi*ones;
// 	    else
// 	      noise_precision&=est_node_phi*ones;
	  }

	noise_precision=reshape(noiseprec_tmp, noiseprec_tmp.Nrows()*noiseprec_tmp.Ncols(),1);	    
      }
    
//     OUT(subject.get_subject_model().get_voxelwise_bold().Nrows());
//     OUT(subject.get_subject_model().get_voxelwise_bold().Ncols());
//     OUT(size(noise_precision));
//     write_ascii_matrix(noise_precision, LogSingleton::getInstance().appendDir("noise_precision"));
    
//     noise_precision_diag.ReSize(noise_precision.Nrows());
//     for( int r=1; r<=noise_precision.Nrows(); r++)
//       {
// 	noise_precision_diag(r,r)=noise_precision(r);
//       }
// //     noise_precision_diag << diag(noise_precision);
//     OUT(size(noise_precision_diag));

  }

  void calculate_forward_model(const Subject& subject, ColumnVector& forward_model)
  {
    Tracer_Plus trace("calculate_forward_model");    

    if(subject.is_single_timeseries() || (!subject.is_single_timeseries() && subject.is_decode()))    
      {
	const vector<ColumnVector>& node_bold = subject.get_subject_model().get_node_bold();

	for(unsigned int r=0; r<node_bold.size(); r++)
	  {
	    if(r==0)
	      forward_model=node_bold[r];
	    else
	      forward_model&=node_bold[r];
	  }
      }
    else
      {
	const Matrix& voxelwise_bold = subject.get_subject_model().get_voxelwise_bold();

	//	OUT(size(voxelwise_bold));
	forward_model=reshape(voxelwise_bold, voxelwise_bold.Nrows()*voxelwise_bold.Ncols(),1);

	//	OUT(size(forward_model));

// 	for(int n=0; n<voxelwise_bold.Nrows(); n++)
// 	  {
// 	    RowVector tmp=voxelwise_bold.Row(n+1);
// 	    if(n==0)
// 	      forward_model=tmp.t();
// 	    else
// 	      forward_model&=tmp.t();
// 	  }

// 	OUT(size(forward_model));
      }

  }

  const ColumnVector& Subject::get_vb_data()
  {
    Tracer_Plus trace("get_vb_data");
   
    if(!already_setup_vb_data)
      {
	already_setup_vb_data=true;

	if(is_single_timeseries())    
	  {
	    const vector<ColumnVector>& node_data = get_node_data();

	    for(unsigned int r=0; r<node_data.size(); r++)
	      {
		if(r==0)
		  vb_data=node_data[r];
		else
		  vb_data&=node_data[r];
	      }
	  }
	else if(is_decode())
	  {
	    const vector<ColumnVector>& node_data = get_subject_model().get_decoded_node_data();
	    
	    for(unsigned int r=0; r<node_data.size(); r++)
	      {
		if(r==0)
		  vb_data=node_data[r];
		else
		  vb_data&=node_data[r];
	      }	
	  }
	else
	  {
	    const Matrix& voxelwise_data = get_voxelwise_data();
	    
	    vb_data=reshape(voxelwise_data, voxelwise_data.Nrows()*voxelwise_data.Ncols(),1);
	    
	    // 	for(int n=0; n<voxelwise_data.Nrows(); n++)
	    // 	  {
	    // 	    RowVector tmp=voxelwise_data.Row(n+1);
	    // 	    if(n==0)
	    // 	      vb_data=tmp.t();
	    // 	    else
	    // 	      vb_data&=tmp.t();
	    // 	  }
	  }
      }

    return vb_data;

  }
  
  void evaluate_model_voxelwise_energy(const Subject& subject, int vox, int debuglevel)
  {
    Tracer_Plus trace("evaluate_model_energy");

    const Subject_Model& subject_model=subject.get_subject_model();

    const Matrix& voxelwise_bold = subject_model.get_voxelwise_bold();
    const Matrix& voxelwise_data = subject.get_voxelwise_data();
    const vector<ColumnVector>& voxelwise_pvf = subject_model.get_voxelwise_pvf();
	
    RowVector tmp=voxelwise_bold.Row(vox);
    double ss = SumSquare(tmp-voxelwise_data.Row(vox));
	    
    double variance;
    if(subject.is_phi_every_voxel())
      variance = 1.0/std::exp(subject_model.get_log_phi_voxel()(1)+subject_model.get_log_phi_every_voxel()(vox));
    else
      variance = 1.0/std::exp(subject_model.get_log_phi_voxel()(1));		

    for(unsigned int r=0; r<voxelwise_pvf.size(); r++)
      {	       	     
	variance += 1.0/subject_model.get_phi_node()[r]*Sqr(voxelwise_pvf[r](vox));
      }

    double voxenergy=0.5*tmp.Ncols()*std::log(variance)+0.5*ss/variance;

    subject_model.set_voxelwise_energy(vox, voxenergy);	  	
  }

  double addup_model_voxelwise_energy(const Subject& subject, int debuglevel)
  {
    Tracer_Plus trace("evaluate_model_energy");

    const Subject_Model& subject_model=subject.get_subject_model();
    const Matrix& voxelwise_bold = subject_model.get_voxelwise_bold();

    double energy=0;

    for(int i=1; i<=voxelwise_bold.Nrows(); i++) // loop through voxels	  
      {
	energy += subject_model.get_voxelwise_energy(i);
      }

    //    energy += calculate_model_prior_energy(subject, debuglevel);

    return energy;
  }

  double calculate_model_energy(const Subject& subject, int debuglevel)
  {
    Tracer_Plus trace("calculate_model_energy");

    double energy=0;

    const Subject_Model& subject_model=subject.get_subject_model();
    
    if(subject.is_single_timeseries())    
      {
	const vector<ColumnVector>& node_bold = subject_model.get_node_bold();
	const vector<ColumnVector>& node_data = subject.get_node_data();

// 	for(unsigned int n=1; n<=node_data.size(); n++)      
// 	  {
// 	    write_ascii_matrix(node_bold[n-1],LogSingleton::getInstance().appendDir("node_bold"+num2str(n)));  
// 	    write_ascii_matrix(node_data[n-1],LogSingleton::getInstance().appendDir("node_data"+num2str(n)));   
// 	    exit(0);
// 	  }

	for(unsigned int r=0; r<node_bold.size(); r++)
	  {

	    double ss = SumSquare(node_bold[r]-node_data[r]);
	    energy += node_bold[r].Nrows()/2.0*std::log(ss);
	  }
      }
    else if(subject.is_decode())
      {
	const vector<ColumnVector>& node_bold = subject_model.get_node_bold();
	const vector<ColumnVector>& decoded_node_data = subject.get_subject_model().get_decoded_node_data();

	for(unsigned int r=0; r<node_bold.size(); r++)
	  {
	    ColumnVector tmp_bold=node_bold[r];
	    //	    tmp_bold/=(stdev(tmp_bold).AsScalar());
	    ColumnVector tmp_data=decoded_node_data[r];
	    //	    tmp_data/=(stdev(tmp_data).AsScalar());

	    //	    double ss = SumSquare(node_bold[r]-decoded_node_data[r]);
	    double ss = SumSquare(tmp_bold-tmp_data);
	    energy += node_bold[r].Nrows()/2.0*std::log(ss);
	  }

// 	const vector<ColumnVector>& node_bold = subject_model.get_node_bold();
// 	const vector<ColumnVector>& node_data = subject.get_subject_model().get_decoded_node_data();
// 	const vector<ColumnVector>& voxelwise_pvf = subject_model.get_voxelwise_pvf();
// 	const Matrix& voxelwise_data = subject.get_voxelwise_data();

// 	for(unsigned int r=0; r<node_bold.size(); r++)
// 	  {
// 	    float sumvar=0;
// 	    float sumpvf=0;

// 	    for(int i=1; i<=voxelwise_data.Nrows(); i++)
// 	      {
// 		if(subject_model.is_voxel_in_roi(i,r+1))
// 		  {

// 		    //	    sumvar+=Sqr(voxelwise_pvf[r](i))/subject_model.get_phi_every_voxel()(i);
// 		    sumvar+=Sqr(voxelwise_pvf[r](i))/subject_model.get_phi_node()[r];
// 		    sumpvf+=voxelwise_pvf[r](i);
// 		  }		
// 	      }

// 	    //double variance = sumvar/Sqr(sumpvf);
// 	    double variance = sumvar;
// 	    double ss = SumSquare(node_bold[r]-node_data[r]);
// 	    energy +=  0.5*(node_bold[r].Nrows())*std::log(variance)+0.5*ss/variance;
// 	  }

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

	    for(unsigned int r=0; r<voxelwise_pvf.size(); r++)
	      {	       	     
		variance += 1.0/subject_model.get_phi_node()[r]*Sqr(voxelwise_pvf[r](i));
	      }

	    double voxenergy=0.5*tmp.Ncols()*std::log(variance)+0.5*ss/variance;
	    energy += voxenergy;
	    subject_model.set_voxelwise_energy(i, voxenergy);
	  }	

// 	for(int n=0; n<voxelwise_bold.Nrows(); n++)
// 	  {
// 	    RowVector tmp=voxelwise_bold.Row(n+1);
// 	    double ss = SumSquare(tmp-voxelwise_data.Row(n+1));
	    
// 	    energy += tmp.Ncols()/2.0*std::log(ss);
// 	  }

      }

//     if(isnan(energy))
//       {
// 	LOGOUT("Invalid model");
		   
// 	const vector<ColumnVector>& node_bold = subject.get_subject_model().get_node_bold();
// 	const vector<ColumnVector>& node_data = subject.get_node_data();
// 	for(unsigned int n=1; n<=node_data.size(); n++)      
// 	  {
// 	    write_ascii_matrix(node_bold[n-1],LogSingleton::getInstance().appendDir("node_bold"+num2str(n)));  
// 	    write_ascii_matrix(node_data[n-1],LogSingleton::getInstance().appendDir("node_data"+num2str(n)));   
// 	  }

// 	exit(0);
	
//       }
    
//    energy += calculate_model_prior_energy(subject, debuglevel);

    return energy;
  }

  double calculate_model_prior_energy(const Subject& subject, int debuglevel)
  {
    Tracer_Plus trace("calculate_model_prior_energy");

    double energy=0;

    const Model& model=subject.get_subject_model().get_model();
    const Subject_Model& subject_model=subject.get_subject_model();
    
    if(subject.is_phi_every_voxel() && !(subject.is_single_timeseries() || subject.is_decode()))
      {
	
	// voxelwise ARD prior
	const Matrix& voxelwise_bold = subject_model.get_voxelwise_bold();
	double ss=0;

	for(int i=1; i<=voxelwise_bold.Nrows(); i++) // loop through voxels	  
	  {
	    ss += Sqr(subject_model.get_log_phi_every_voxel()(i));	   
	  }
	energy += voxelwise_bold.Nrows()/2.0*std::log(ss);

      }

    //////////
    // now do amplitude modulation ARD prior

    // b amp mod components
    for(unsigned int n=1; n<=model.get_nodes().size(); n++)
      for(unsigned int i=1; i<=model.get_marker_b()[n-1].size(); i++)
	{
	  if(model.is_b_amp_mod(n,i))
	    {
	      //LOGOUT(subject_model.get_value_b()[n-1][i-1]);
	      double ss=0;
	      int stim_index=model.get_marker_b()[n-1][i-1].second;
	      for(unsigned int e=1; e<=subject_model.get_stimuli()[stim_index-1].size(); e++)
		{
		  //LOGOUT(subject_model.get_value_b_amp_mod()[n-1][i-1][e-1]);
		  ss += Sqr(subject_model.get_value_b_amp_mod()[n-1][i-1][e-1]);	   
		}  

	      energy += subject_model.get_stimuli()[stim_index-1].size()/2.0*std::log(ss);	    
	    }
	}

    // c amp mod components
    for(unsigned int n=1; n<=model.get_nodes().size(); n++)
      for(unsigned int i=1; i<=model.get_marker_c()[n-1].size(); i++)
	{
	  if(model.is_c_amp_mod(n,i))
	    {
	      double ss=0;
	      int stim_index=model.get_marker_c()[n-1][i-1];
	      for(unsigned int e=1; e<=subject_model.get_stimuli()[stim_index-1].size(); e++)
		{
		  ss += Sqr(subject_model.get_value_c_amp_mod()[n-1][i-1][e-1]);	   
		}  
	      energy += subject_model.get_stimuli()[stim_index-1].size()/2.0*std::log(ss);	    
	    }
	}

    ////////////
    // connection parameters ARD priors
    
    //     vector<float> values = subject.get_subject_model().get_value_a_vec();
    //     for(unsigned int p=0; p<values.size(); p++)
    //       energy += std::log(std::abs(values[p]));
    
    //     values = subject.get_subject_model().get_value_b_vec();
    //     for(unsigned int p=0; p<values.size(); p++)
    //       energy += std::log(std::abs(values[p]));
    
    
    //     values = subject.get_subject_model().get_value_c_vec();
    //     for(unsigned int p=0; p<values.size(); p++)     
    //       energy += std::log(std::abs(values[p]));  	  
    
    int nnodes=subject_model.get_subject_nodes().size();  
    for(int n=1; n<=nnodes; n++)
      for(unsigned int i=1; i<=subject.get_subject_model().get_value_c()[n-1].size(); i++)
	{	     
	  if(model.is_c_ard(n,i))
	    {
	      energy += std::log(std::abs(subject.get_subject_model().get_value_c()[n-1][i-1]));	      
	    }
	}

    //     values = subject.get_subject_model().get_value_d_vec();
    //     for(unsigned int p=0; p<values.size(); p++)
    //       energy += std::log(std::abs(values[p])); 
    
    ////////////////////// 	
    
    //     const vector<vector<float> >& value_a=subject_model.get_value_a();
    //     const vector<vector<float> >& value_b=subject_model.get_value_b();
    //     const vector<vector<float> >& value_d=subject_model.get_value_d();
    //     const vector<vector<int> >& marker_a=subject_model.get_model().get_marker_a();
    //     const vector<vector<pair<int,int> > >& marker_b=subject_model.get_model().get_marker_b();
    //     const vector<vector<pair<int,int> > >& marker_d=subject_model.get_model().get_marker_d();
    
    //     int nnodes=subject_model.get_subject_nodes().size();    
    
    /////////////////////////////
    //     // sum squares of all a_{ij} pairs needs to be less than nnodes/(nnodes-1) (see appendix a.3 in dcm paper)
    //     float ss=0.0;    
    //     for(int n=0; n<nnodes; n++)
    //       {
    // 	ColumnVector ab(nnodes);
    // 	ab=0;
    // 	for(unsigned int i=0; i<marker_a[n].size(); i++)
    // 	  {
    // 	    int node_index=marker_a[n][i];
    // 	    ab(node_index)+=value_a[n][i];
    // 	  }
    // 	for(unsigned int i=0; i<marker_b[n].size(); i++)
    // 	  {    
    // 	    int node_index=marker_b[n][i].first;
    // 	    ab(node_index)+=value_b[n][i];
    // 	  }
    // 	for(unsigned int i=0; i<marker_d[n].size(); i++)
    // 	  {    
    // 	    int node_index=marker_d[n][i].first;
    // 	    ab(node_index)+=value_d[n][i];
    // 	  }
    
    // 	ss+=SumSquare(ab);
    //       }
    
    //     float boundary=nnodes/float(nnodes-1);
    
    //if(!subject.is_single_timeseries())    
    //{
    // 	OUT(boundary);
    // 	OUT(ss);
    //	OUT(-std::log(1-1.0/(1+std::exp(-25*(ss-(boundary))))));	
      //}
    
    //    x=0:0.01:2;figure;plot(x,1-1./(1+exp(-25*(x-(1-0.1))));
    
    //    energy+=-std::log(1-1.0/(1+std::exp(-25*(ss-(boundary-0.2)))));
    //energy+=-std::log(1-1.0/(1+std::exp(-10*(ss-(boundary-0.1)))));
    
    //     if(ss>boundary) 
    //       {
    // //     	LOGOUT(ss);
    // //     	LOGOUT(boundary);
    
    //     	energy=1e32;
    //     	return energy;
    //       }
    //////////////////////////
    
//////////////////////////
//     // sum squares of all a_{ij} pairs needs to be less than nnodes/(nnodes-1) (see appendix a.3 in dcm paper)
//     float ss=0.0;    
//     for(int n=0; n<nnodes; n++)
//       {
// 	ColumnVector ab(nnodes);
// 	ab=0;
// 	for(unsigned int i=0; i<marker_a[n].size(); i++)
// 	  {
// 	    int node_index=marker_a[n][i];
// 	    ab(node_index)+=value_a[n][i];
// 	  }
// 	ss+=SumSquare(ab);
//       }
//     float boundary=nnodes/float(nnodes-1);
    
//     //    x=0:0.01:2;figure;plot(x,1-1./(1+exp(-25*(x-(1-0.1)))));
//     //    energy+=-std::log(1-1.0/(1+std::exp(-25*(ss-(boundary-0.2)))));

//     if(ss>boundary) 
//       {
// // 	LOGOUT(ss);
// // 	LOGOUT(boundary);

// 	energy=1e32;
// 	return energy;
//       }

//     // for b
//     ss=0.0;    
//     for(int n=0; n<nnodes; n++)
//       {
// 	ColumnVector ab(nnodes);
// 	ab=0;
// 	for(unsigned int i=0; i<marker_b[n].size(); i++)
// 	  {    
// 	    int node_index=marker_b[n][i].first;
// 	    ab(node_index)+=value_b[n][i];
// 	  }
// 	ss+=SumSquare(ab);
//       }
//     //    boundary=nnodes/float(nnodes-1);

//     //    x=0:0.01:2;figure;plot(x,1-1./(1+exp(-25*(x-(1-0.1)))));
//     //    energy+=-std::log(1-1.0/(1+std::exp(-25*(ss-(boundary-0.2)))));

//     if(ss>boundary) 
//       {
// // 	LOGOUT(ss);
// // 	LOGOUT(boundary);
// 	energy=1e32;
// 	return energy;
//       }
    
//     // for d
//     ss=0.0;    
//     for(int n=0; n<nnodes; n++)
//       {
// 	ColumnVector ab(nnodes);
// 	ab=0;
// 	for(unsigned int i=0; i<marker_d[n].size(); i++)
// 	  {    
// 	    int node_index=marker_d[n][i].first;
// 	    ab(node_index)+=value_d[n][i];
// 	  }
// 	ss+=SumSquare(ab);
//       }
//     boundary=nnodes/float(nnodes-1);

//     //    x=0:0.01:2;figure;plot(x,1-1./(1+exp(-25*(x-(1-0.1)))));
//     //    energy+=-std::log(1-1.0/(1+std::exp(-25*(ss-(boundary-0.2)))));

//     if(ss>boundary) 
//       {
// // 	LOGOUT(ss);
// // 	LOGOUT(boundary);
// 	energy=1e32;
// 	return energy;
//       }
    /////////////////////

//     // sigmaa
//     energy += Sqr(subject_model.get_value_sigmaa()-1)/(2.0*Sqr(0.2));
    
//     ////////////
//     // do priors on non-neuronal connectivity params

//     for(int n=1; n<=nnodes; n++)
//       {
// 	Subject_Node& subject_node=*(subject_model.get_subject_nodes()[n-1]);

// 	// HRFs
// 	if(subject_model.get_haemodynamic_model()=="balloon" || subject_model.get_haemodynamic_model()=="balloon_epsilon")      
// 	  {
// 	    energy+=subject_node.get_balloon_cbf_prior_energy(vector2ColumnVector(subject_node.get_value_balloon_cbf()));

// 	    energy+=subject_node.get_balloon_prior_energy(vector2ColumnVector(subject_node.get_value_balloon()));

// 	    if(subject_model.get_haemodynamic_model()=="balloon_epsilon")
// 	      energy+=subject_node.get_balloon_prior_energy2(vector2ColumnVector(subject_node.get_value_balloon2()));
// 	  }
// 	else
// 	  {
// 	    vector<float> values = subject_node.get_value_hrf();
// 	    vector<float> mins = subject_node.get_mins_hrf();
// 	    vector<float> maxs = subject_node.get_maxs_hrf();
// 	    for(unsigned int p=0; p<values.size(); p++)
// 	      {
// 		float mn=(mins[p]+maxs[p])/2.0;
// 		float st=(maxs[p]-mins[p])/2.0;
// 		energy += Sqr(values[p]-mn)/(2.0*Sqr(st));
// 	      }
// 	  }

// // 	// MVN off-diag covariance ARD priors
// // 	if(!subject.is_single_timeseries())    
// // 	  {
// // 	    values = subject_node.get_value_mvn_sqrt_cov_offdiag();
// // 	    if(values.size()>0)
// // 	      {
// // 		float ss=0;
// // 		for(unsigned int p=0; p<values.size(); p++)
// // 		  {
// // 		    ss += Sqr(values[p]);	   
// // 		  }  
// // 		energy += values.size()/2.0*std::log(ss);
// // 	      }
// // 	  }
//       }

//     // put in penalty for roi's overlapping too much
//     if(!subject.is_single_timeseries())    
//       for(int n=1; n<=nnodes; n++)
// 	{
// 	  Subject_Node& subject_node=*(subject_model.get_subject_nodes()[n-1]);
// 	  ColumnVector mean1=subject_node.get_mvn_mean_func_space_value();
// 	  for(int n2=n+1; n2<=nnodes; n2++)
// 	    {
// 	      Subject_Node& subject_node2=*(subject_model.get_subject_nodes()[n2-1]);
	      
// 	      ColumnVector mean2=subject_node2.get_mvn_mean_func_space_value();
	      
// 	      // P(MVN_MEAN)
// 	      ColumnVector tmp=mean1-mean2;
	      
// 	      float distance=sqrt(sum(SP(tmp,tmp)).AsScalar());
	      
// 	      //x=1:20;figure;plot(x,1./(1+exp(-10*(x-10))));
// 	      energy+=std::log(1+std::exp(-10*(distance-15.0)));

// 	      // P(MVN_COV|MVN_MEAN)
// // 	      ColumnVector mid_coord=(mean1+mean2)/2;
	      
// // 	      float measure=std::sqrt(subject_node2.get_partial_vol_fraction(mid_coord)*subject_node.get_partial_vol_fraction(mid_coord));
	      
// // 	      // x=0:0.1:1;figure;plot(x,1-1./(1+exp(-25*(x-0.2))));
// // 	      energy+=-std::log(1-1.0/(1+std::exp(-25*(measure-0.2))));
	      
// 	    }
// 	}
    
    //     if(debuglevel==3)
    //       {
    // 	LOGOUT(energy);
    //       }

    return energy;
  }
  
  
  void create_report_fit(Subject& subject, string base_name)
  {
    Tracer_Plus trace("create_report_fit");
    
    Subject_Model& subject_model = subject.get_subject_model();

    //int number;

//     LogSingleton::getInstance().str()<< "Creating report fit"<< endl;

    ///////////////////
    // setup html report file
    string logfilename=LogSingleton::getInstance().getLogFileName();
    LogSingleton::getInstance().setLogFile(string("report_")+base_name+string(".html"));
    LogSingleton::getInstance().set_stream_to_cout(false);
    Log& htmllog = LogSingleton::getInstance();
	
    htmllog << "<HTML> " << endl
	    << "<TITLE>NMA " << base_name << " fits for " << subject.get_name() << "</TITLE>" << endl
	    << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	    << "/doc/images/fsl-bg.jpg\">" << endl 
	    << "<hr><CENTER><H1>NMA " << base_name << " fits for <br>" << subject.get_name() << " </H1>"<< endl
	    << "<hr><p>" << endl;
    
    // output model energy
    double energy = calculate_model_energy(subject, 0);
    htmllog << "Energy=" << energy << "<p>" << endl;

    //  LOGOUT(energy);
 
    Nma_Mcmc_Log_Likelihood nma_mcmc_log_likelihood(subject);
    
    htmllog << "Likelihood=" << nma_mcmc_log_likelihood.evaluate() << "<p>" << endl;
   
    int nnodes=subject.get_subject_model().get_subject_nodes().size();
    for(int n=1; n<=nnodes; n++)      
      {
	htmllog << "Likelihood for " << subject.get_subject_model().get_subject_nodes()[n-1]->get_node().get_name() << "=" << nma_mcmc_log_likelihood.evaluate(n) << "; " << endl;
      }
    htmllog << "<p>" << endl;

//     /////////////////////
//     // write HRFs
//     vector<ColumnVector> hrf;
//     ColumnVector t_hrf;
//     subject_model.get_hrf(hrf,t_hrf,25);
//     write_ascii_matrix(t_hrf,LogSingleton::getInstance().appendDir(base_name+"_t_hrf.txt"));

//     for(int n=1; n<=nnodes; n++)      
//       {	
// 	miscplot tsplot;
// 	tsplot.add_xlabel("seconds");    
// 	tsplot.set_xysize(300,150);
// 	tsplot.set_minmaxscale(1);
// 	tsplot.timeseries(hrf[n-1].t(),LogSingleton::getInstance().appendDir(base_name+"_hrf"+num2str(n)), string("HRF ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), t_hrf(1),400,3,0,false);	    
// 	htmllog << "<img BORDER=0 SRC=\"" <<  base_name+"_hrf"+num2str(n) << ".png\">"<< endl;
// 	write_ascii_matrix(hrf[n-1].t(), LogSingleton::getInstance().appendDir(base_name+"_hrf"+num2str(n)+".txt"));
				 
//       } 
//     htmllog << "<p>" << endl;

    /////////////////////
    // write neural response functions
    vector<vector<ColumnVector> > nrf;
    ColumnVector t_nrf;
    subject_model.get_nrf(nrf,t_nrf,25);
    write_ascii_matrix(t_nrf,LogSingleton::getInstance().appendDir(base_name+"_t_nrf.txt"));
    for(unsigned int s=1; s<=nrf.size(); s++)      
      {	
	htmllog << "Neural response functions for stimulus " << subject_model.get_model().get_stimuli_names()[s-1] << endl;
	htmllog << "<p>" << endl;
	for(int n=1; n<=nnodes; n++)      
	  {
	    string tmp_name=base_name+"_nrf_stim"+num2str(s)+"_node"+num2str(n);

	    miscplot tsplot;
	    tsplot.add_xlabel("seconds");    
	    tsplot.set_xysize(300,150);
	    tsplot.set_minmaxscale(1);
	    tsplot.timeseries(nrf[s-1][n-1].t(),LogSingleton::getInstance().appendDir(tmp_name), string("NRF ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), t_nrf(1),400,3,0,false);	    
	    htmllog << "<img BORDER=0 SRC=\"" <<  LogSingleton::getInstance().appendDir(tmp_name + ".png\">") << endl;
	    write_ascii_matrix(nrf[s-1][n-1].t(), LogSingleton::getInstance().appendDir(tmp_name+".txt"));
	    
	  } 
	htmllog << "<p>" << endl;
      }

    if(!subject.is_single_timeseries())// && (opts.data_mode.value()==1 || opts.data_mode.value()==3))
      {

	/////////////////////
	// write ROI shape and position and voxelwise bold
	const volume4D<float>& data=subject.get_data();
	volume<float> pvf;
	pvf=subject.get_func_mask();
	pvf=0;

	//	volume4D<float> voxelwise_bold_4D(data.xsize(),data.ysize(),data.zsize(),data.tsize());
	//	voxelwise_bold_4D=0;
	const Matrix& voxelwise_bold = subject_model.get_voxelwise_bold();
	const Matrix& voxelwise_data = subject.get_voxelwise_data();

	vector<float> max_pvf(nnodes,0.0);
	vector<float> sum_pvf(nnodes,0.0);
	vector<ColumnVector> max_bold(nnodes);
	vector<ColumnVector> max_data(nnodes);
	vector<ColumnVector> average_bold(nnodes);
	vector<ColumnVector> average_data(nnodes);
	vector<ColumnVector> max_pvf_coords(nnodes);
	vector<volume<float> > pvf_separate(nnodes);

	for( int r=0; r<nnodes; r++) // indexes nodes
	  {
	    pvf_separate[r]=subject.get_func_mask();
	    pvf_separate[r]=0;

	    average_bold[r].ReSize(data.tsize());
	    average_bold[r]=0;
	    average_data[r].ReSize(data.tsize());
	    average_data[r]=0;
	    max_pvf[r]=0.0;
	    sum_pvf[r]=0.0;
	  }

	const vector<ColumnVector>&  coords=subject_model.get_voxel_coordinates();
	    	    
	for(unsigned int n=0; n<coords.size(); n++) // indexes voxels
	  {
	    RowVector tmp=voxelwise_bold.Row(n+1);
	    //	    voxelwise_bold_4D.setvoxelts(tmp.t(),int(coords[n](1)),int(coords[n](2)),int(coords[n](3)));

	    for(int r=0; r<nnodes; r++) // indexes nodes
	      {		
		if(subject_model.is_voxel_in_roi(n+1,r+1))
		  {
		    pvf(int(coords[n](1)),int(coords[n](2)),int(coords[n](3)))+=subject_model.get_partial_vol_fraction(n+1,r+1);

	// 	    if(r==1)
// 		      subject_model.print_partial_vol_fraction(n+1,r+1);

		    pvf_separate[r](int(coords[n](1)),int(coords[n](2)),int(coords[n](3)))=subject_model.get_partial_vol_fraction(n+1,r+1);

		    float pvf_val=pvf_separate[r](int(coords[n](1)),int(coords[n](2)),int(coords[n](3)));

		    // addup weighted average fit:
		    average_bold[r]+=pvf_val*tmp.t();
		    average_data[r]+=pvf_val*voxelwise_data.Row(n+1).t();
		    sum_pvf[r]+=pvf_val;

		    if(max_pvf[r]<pvf_val)
		      {
			max_pvf[r]=pvf_val;
			max_pvf_coords[r]=coords[n];
			max_bold[r]=tmp.t();
			max_data[r]=voxelwise_data.Row(n+1).t();
		      }
		  }
	      }
	    
	  }

	for( int r=0; r<nnodes; r++) // indexes nodes
	  {
	    average_bold[r]/=sum_pvf[r];
	    average_data[r]/=sum_pvf[r];
	  }

	vector<volume<float> > pvf_std(nnodes);

	for( int r=0; r<nnodes; r++) // indexes nodes
	  {
	   
	    save_volume(pvf_separate[r], LogSingleton::getInstance().appendDir("pvf_"+num2str(r+1)));

	    // convert into std space and save those out too
	    
	    const Subject_Node& subject_node = *(subject.get_subject_model().get_subject_nodes()[r]);

	    pvf_std[r]=subject_node.get_std_mask();

	    affine_transform(pvf_separate[r],pvf_std[r],subject_node.std2func_xform.i());
	 
   	    save_volume(pvf_std[r], LogSingleton::getInstance().appendDir("pvf_std_"+num2str(r+1)));

	  }      

	if(subject.is_decode())
	  {
	    /////////////////////
	    // write decoded node bold
	    const vector<ColumnVector>& node_bold=subject_model.get_node_bold();
	    const vector<ColumnVector>& node_data=subject_model.get_decoded_node_data();
	    
	    for(int n=1; n<=nnodes; n++)      
	      {
		string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
		string name=base_name+"_"+node_name+"_bold";
		write_ascii_matrix(node_bold[n-1],LogSingleton::getInstance().appendDir(name));    
		
		ColumnVector tmp_bold=node_bold[n-1];
		//tmp_bold/=(stdev(tmp_bold).AsScalar());
		ColumnVector tmp_data=node_data[n-1];
		//		tmp_data/=(stdev(tmp_data).AsScalar());

		miscplot tsplot;
		tsplot.add_xlabel("TRs");    
		tsplot.set_xysize(300,200);
		tsplot.set_minmaxscale(1);
		Matrix tmp=tmp_data.t() & tmp_bold.t();
		tsplot.timeseries(tmp,LogSingleton::getInstance().appendDir(name), string("node bold ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), 0,400,3,0,false);	    
		htmllog << "<img BORDER=0 SRC=\"" <<  name << ".png\">" << endl;
		write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(name+".txt"));
		
	      }  
	    
	    htmllog << "<p>" << endl;    
	  }
	else
	  {
	    /////////////////////
	    // write voxelwise bold fit and data at max voxel	
	    for(int n=1; n<=nnodes; n++)      
	      {
		const Subject_Node& subject_node=*(subject_model.get_subject_nodes()[n-1]);   
		
		miscplot tsplot;
		tsplot.add_xlabel("TRs");    
		tsplot.set_xysize(300,200);
		tsplot.set_minmaxscale(1);
		Matrix tmp=max_data[n-1].t() & max_bold[n-1].t();
		write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(base_name+"_bold"+num2str(n)+".txt"));
		tsplot.timeseries(tmp,LogSingleton::getInstance().appendDir(base_name+"_bold"+num2str(n)), string("Max voxel fit ")+subject_node.get_node().get_name(), 0,400,3,0,false);	    
		htmllog << "<img BORDER=0 SRC=\"" <<  base_name+"_bold"+num2str(n) << ".png\">" << endl;
		
		
	      }   
	    htmllog << "<p>" << endl;    	    

	    /////////////////////
	    // write weighted average bold fit and data	
	    for(int n=1; n<=nnodes; n++)      
	      {
		const Subject_Node& subject_node=*(subject_model.get_subject_nodes()[n-1]);   
		
		miscplot tsplot;
		tsplot.add_xlabel("TRs");    
		tsplot.set_xysize(300,200);
		tsplot.set_minmaxscale(1);
		Matrix tmp=average_data[n-1].t() & average_bold[n-1].t();
		write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(base_name+"_average_bold"+num2str(n)+".txt"));
		tsplot.timeseries(tmp,LogSingleton::getInstance().appendDir(base_name+"_average_bold"+num2str(n)), string("Average voxel fit ")+subject_node.get_node().get_name(), 0,400,3,0,false);	    
		htmllog << "<img BORDER=0 SRC=\"" <<  base_name+"_average_bold"+num2str(n) << ".png\">" << endl;
		
		
	      }   
	    htmllog << "<p>" << endl;    
	  }	

	/////////////////////
	// write images of ROI shapes and positions in func space and std space
	const volume<float>& epivol=subject.get_data()[int(data.tsize()/2.0)];
	//	volumeinfo& data_volinfo=subject.get_data_volinfo();

	for(int n=1; n<=nnodes; n++)      
	  {
	    {	 
	      // func space
	      string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
	      string name=base_name+"_"+node_name+"_func_pvf";
	      create_pvf_overlays(name, pvf_separate[n-1], epivol, subject.get_func_mask(), max_pvf_coords[n-1]);//, data_volinfo);
	      htmllog << node_name << "<br>" << endl;	    
	      htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_x.png") << "\"> " << endl;
	      htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_y.png") << "\">" << endl;
	      htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_z.png") << "\">" << endl;
	      htmllog << "<p>" << endl; 
	    }
	    
// 	    if(!(subject_model.get_subject_nodes()[n-1]->std2func_xform(1,1)==1 && subject_model.get_subject_nodes()[n-1]->std2func_xform(2,2)==1 && subject_model.get_subject_nodes()[n-1]->std2func_xform(3,3)==1))
// 	    {
// 	      // std space
// 	      string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
// 	      string name=base_name+"_"+node_name+"_std_pvf";
// 	      ColumnVector max_pvf_std_coords = xform_coords(max_pvf_coords[n-1],subject_model.get_subject_nodes()[n-1]->std2func_xform.i(),subject.get_funcvoxdim(),subject.get_stdvoxdim());
// 	      create_pvf_overlays(name, pvf_std[n-1], subject.get_standard_brain(), subject.get_standard_brain(), max_pvf_std_coords);
// 	      htmllog << node_name << "<br>" << endl;	    
// 	      htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_x.png") << "\"> " << endl;
// 	      htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_y.png") << "\">" << endl;
// 	      htmllog << "<img BORDER=0 SRC=\"" <<  name+string("_z.png") << "\">" << endl;
// 	      htmllog << "<p>" << endl; 
// 	    }

 	  }

 	htmllog << "<p>" << endl; 

	//	save_volume4D(voxelwise_bold_4D, LogSingleton::getInstance().appendDir(base_name+"_voxelwise_bold"));	  
	
	save_volume(pvf, LogSingleton::getInstance().appendDir(base_name+"_pvf"));	  

      }
    else
      {
	/////////////////////
	// write node bold
	const vector<ColumnVector>& node_bold=subject_model.get_node_bold();
	const vector<ColumnVector>& node_data=subject.get_node_data();

	for(int n=1; n<=nnodes; n++)      
	  {
	    string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
	    string name=base_name+"_"+node_name+"_bold";
	    write_ascii_matrix(node_bold[n-1],LogSingleton::getInstance().appendDir(name));    
	    
	    miscplot tsplot;
	    tsplot.add_xlabel("TRs");    
	    tsplot.set_xysize(300,200);
	    tsplot.set_minmaxscale(1);
	    Matrix tmp=node_data[n-1].t() & node_bold[n-1].t();
	    tsplot.timeseries(tmp,LogSingleton::getInstance().appendDir(name), string("node bold ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), 0,400,3,0,false);	    
	    htmllog << "<img BORDER=0 SRC=\"" <<  name << ".png\">" << endl;
	    write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(name+".txt"));

	  }  
 
	htmllog << "<p>" << endl;    
      }

    /////////////////////
    // write node z
    const vector<ColumnVector>& node_z=subject_model.get_z();

    for(int n=1; n<=nnodes; n++)      
      {	    
	string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
	string name=base_name+"_"+node_name+"_z";
	    
	miscplot tsplot;
	tsplot.add_xlabel("");    
	tsplot.set_xysize(300,200);
	tsplot.set_minmaxscale(1.1);
	Matrix tmp=node_z[n-1].t();
	//OUT(tmp);
	tsplot.timeseries(tmp,LogSingleton::getInstance().appendDir(name), string("neuronal activity ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), 0,400,3,0,false);	    
	htmllog << "<img BORDER=0 SRC=\"" <<  name << ".png\">" << endl;
	write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(name+".txt"));

      }  

    htmllog << "<p>" << endl;  

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


  void create_pvf_overlays(const string& name, const volume<float>& pvf, const volume<float>& bgvol, const volume<float>& bgvolmask, const ColumnVector& coords)//, volumeinfo& data_volinfo)
  {
    Tracer_Plus trace("create_pvf_overlays");
    
    volume<float> newvol; 
    miscpic newpic;
    volume<float> bgvoltmp = bgvol;
    volume<float> pvftmp=pvf;

    pvftmp.threshold(0.001);

//     volumeinfo volinfo=data_volinfo;
//     FslSetDim(&volinfo,bgvol.xsize(),bgvol.ysize(),bgvol.zsize(),1);

    save_volume(pvftmp, LogSingleton::getInstance().appendDir("pvftmp"));
    save_volume(bgvoltmp, LogSingleton::getInstance().appendDir("bgvoltmp"));

    save_volume(bgvol, LogSingleton::getInstance().appendDir("bgvol"));
    save_volume(bgvolmask, LogSingleton::getInstance().appendDir("bgvolmask"));

    // overlay pvf onto bg volume
    newpic.overlay(newvol, bgvoltmp, pvftmp, pvftmp, 
		   bgvol.percentile(0.01,bgvolmask),bgvol.percentile(0.99,bgvolmask),
		   //float(0.1), float(1.0), float(0.001), float(1.0),
		   float(0.001), float(1.0), float(0.001), float(1.0),
		   0, 0);//, &volinfo);

//     save_volume(newvol, LogSingleton::getInstance().appendDir("newvol"),volinfo,true);
 
    // now do x
    {
      string instr(" -l render1 -s 3 -x ");
      instr+=num2str(coords(1)/float(bgvol.xsize())) + " " + LogSingleton::getInstance().appendDir(name+string("_x.png"));     
      newpic.set_cbar(string("y"));

      char instrtmp[10000];
      sprintf(instrtmp,instr.data());
      newpic.slicer(newvol, instrtmp);//, &data_volinfo); 	
    }

    // now do y
    {
      string instr(" -l render1 -s 3 -y ");
      instr+=num2str(coords(2)/float(bgvol.ysize())) + " " + LogSingleton::getInstance().appendDir(name+string("_y.png"));     
      
      // 	  string tit = "Activation prob map";
      // 	  if(thresh>0)
      // 	    tit += " thresholded at p>" + num2str(thresh);
      // 	  newpic.set_title(tit);	    
      
      newpic.set_cbar(string("y"));
      char instrtmp[10000];
      sprintf(instrtmp,instr.data());
      newpic.slicer(newvol, instrtmp);//, &data_volinfo); 	
    }
    
    // now do z
    {
      string instr(" -l render1 -s 3 -z ");
      instr+=num2str(coords(3)/float(bgvol.zsize())) + " " + LogSingleton::getInstance().appendDir(name+string("_z.png"));     
      
      // 	  string tit = "Activation prob map";
      // 	  if(thresh>0)
      // 	    tit += " thresholded at p>" + num2str(thresh);
      // 	  newpic.set_title(tit);	    
      
      newpic.set_cbar(string("y"));
      char instrtmp[10000];
      sprintf(instrtmp,instr.data());
      newpic.slicer(newvol, instrtmp);//, &data_volinfo); 	
    }
  }

void create_report_signal(Subject_Model& subject_model, string base_name)
{
  ///////////////////
  // setup html report file
  string logfilename=LogSingleton::getInstance().getLogFileName();
  LogSingleton::getInstance().setLogFile(string("report_")+base_name+string(".html"));
  LogSingleton::getInstance().set_stream_to_cout(false);
  Log& htmllog = LogSingleton::getInstance();
  
  htmllog << "<HTML> " << endl
	  << "<TITLE>NMA " << base_name << " generated data </TITLE>" << endl
	  << "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	  << "/doc/images/fsl-bg.jpg\">" << endl 
	  << "<hr><CENTER><H1>NMA " << base_name << " generated data </H1>" << endl
	  << "<hr><p>" << endl;
  
   int nnodes=subject_model.get_subject_nodes().size();

//   /////////////////////
//   // write HRFs
//   vector<ColumnVector> hrf;
//   ColumnVector t_hrf;
//   subject_model.get_hrf(hrf,t_hrf,25);
//   write_ascii_matrix(t_hrf,LogSingleton::getInstance().appendDir(base_name+"_t_hrf"));
//   for(int n=1; n<=nnodes; n++)      
//     {	
//       miscplot tsplot;
//       tsplot.add_xlabel("seconds");    
// 	tsplot.set_xysize(300,150);
// 	tsplot.set_minmaxscale(1);
// 	tsplot.timeseries(hrf[n-1].t(),LogSingleton::getInstance().appendDir(base_name+"_hrf"+num2str(n)), string("HRF ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), t_hrf(1),400,3,0,false);	    
// 	htmllog << "<img BORDER=0 SRC=\"" <<  base_name+"_hrf"+num2str(n) << ".png\">"<< endl;
// 	write_ascii_matrix(hrf[n-1].t(), LogSingleton::getInstance().appendDir(base_name+"_hrf"+num2str(n)+".txt"));
	
//     } 
//   htmllog << "<p>" << endl;

    /////////////////////
    // write neural response functions
    vector<vector<ColumnVector> > nrf;
    ColumnVector t_nrf;
    subject_model.get_nrf(nrf,t_nrf,25);
    write_ascii_matrix(t_nrf,LogSingleton::getInstance().appendDir(base_name+"_t_nrf.txt"));
    for(unsigned int s=1; s<=nrf.size(); s++)      
      {	
	htmllog << "Neural response functions for stimulus " << subject_model.get_model().get_stimuli_names()[s-1] << endl;
	htmllog << "<p>" << endl;
	for(int n=1; n<=nnodes; n++)      
	  {
	    string tmp_name=base_name+"_nrf_stim"+num2str(s)+"_node"+num2str(n);
	    
	    miscplot tsplot;
	    tsplot.add_xlabel("seconds");    
	    tsplot.set_xysize(300,150);
	    tsplot.set_minmaxscale(1);
	    tsplot.timeseries(nrf[s-1][n-1].t(),LogSingleton::getInstance().appendDir(tmp_name), string("NRF ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), t_nrf(1),400,3,0,false);	    
	    htmllog << "<img BORDER=0 SRC=\"" <<  LogSingleton::getInstance().appendDir(tmp_name + ".png\">") << endl;
	    write_ascii_matrix(nrf[s-1][n-1].t(), LogSingleton::getInstance().appendDir(tmp_name+".txt"));
	    
	  } 
	htmllog << "<p>" << endl;
      }
  
  const vector<ColumnVector>& node_bold=subject_model.get_node_bold();
  
  for(int n=1; n<=nnodes; n++)      
    {	    
      string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
      string name=base_name+"_"+node_name+"_bold";
      
      miscplot tsplot;
      tsplot.add_xlabel("TRs");    
      tsplot.set_xysize(300,200);
      tsplot.set_minmaxscale(1);
      Matrix tmp=node_bold[n-1].t();
      tsplot.timeseries(tmp,LogSingleton::getInstance().appendDir(name), string("node bold ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), 0,400,3,0,false);	    
      htmllog << "<img BORDER=0 SRC=\"" <<  name << ".png\">" << endl;
      write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(name+".txt"));
      
    }  
  
  htmllog << "<p>" << endl;  
  
  /////////////////////
  // write node z
  const vector<ColumnVector>& node_z=subject_model.get_z();
  
  for(int n=1; n<=nnodes; n++)      
    {	    
      string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
      string name=base_name+"_"+node_name+"_z";   
      
      miscplot tsplot;
      tsplot.add_xlabel("");    
      tsplot.set_xysize(300,200);
      tsplot.set_minmaxscale(1.1);
      Matrix tmp=node_z[n-1].t();
      tsplot.timeseries(tmp,LogSingleton::getInstance().appendDir(name), string("neuronal activity ")+subject_model.get_subject_nodes()[n-1]->get_node().get_name(), 0,400,3,0,false);	    
      htmllog << "<img BORDER=0 SRC=\"" <<  name << ".png\">" << endl;
      write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(name+".txt"));
      
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


}
