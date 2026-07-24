/*  node.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT  */

#include "utils/log.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "warpfns/warpfns.h"
#include "newimage/newimageall.h"
#include "utils/tracer_plus.h"
#include "node.h"
#include "mcmc_mh.h"
#include "model.h"
#include <algorithm>

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace std;

namespace Nma {

//    ReturnMatrix xform_coords(const ColumnVector& in, const Matrix& xform, const ColumnVector& in_voxdim, const ColumnVector& out_voxdim)
//   {
//     ColumnVector tmp(4);
//     tmp(1)=in(1)*in_voxdim(1);
//     tmp(2)=in(2)*in_voxdim(2);
//     tmp(3)=in(3)*in_voxdim(3);
//     tmp(4)=1;

//     tmp=xform*tmp;

//     tmp=tmp.Rows(1,3);
//     tmp(1)=tmp(1)/out_voxdim(1);
//     tmp(2)=tmp(2)/out_voxdim(2);
//     tmp(3)=tmp(3)/out_voxdim(3);
    
//     tmp.Release();
//     return tmp;
//   }

//    ReturnMatrix xform_coords(const vector<float>& in, const Matrix& xform, const ColumnVector& in_voxdim, const ColumnVector& out_voxdim)
//   {
//     ColumnVector tmp(4);
//     tmp(1)=in[0]*in_voxdim(1);
//     tmp(2)=in[1]*in_voxdim(2);
//     tmp(3)=in[2]*in_voxdim(3);
//     tmp(4)=1;

//     tmp=xform*tmp;

//     tmp=tmp.Rows(1,3);
//     tmp(1)=tmp(1)/out_voxdim(1);
//     tmp(2)=tmp(2)/out_voxdim(2);
//     tmp(3)=tmp(3)/out_voxdim(3);
    
//     tmp.Release();
//     return tmp;
//   }
  
  Hrf::Hrf(const ColumnVector& pprior_mean, const ColumnVector& pprior_min, const ColumnVector& pprior_max) :    
    prior_mean(pprior_mean),
      prior_min(pprior_min),
      prior_max(pprior_max)
    {
      for(int i=1; i <= prior_mean.Nrows(); i++)
	if(prior_max(i)>prior_min(i))
	  {
	    marker.push_back(i);      
	  }
    }

  void Hrf::set(const ColumnVector& pprior_mean, const ColumnVector& pprior_min, const ColumnVector& pprior_max)    
   
    {
      prior_mean=pprior_mean;
      prior_min=pprior_min;
      prior_max=pprior_max;
      
      marker.clear();
      for(int i=1; i <= prior_mean.Nrows(); i++)
	if(prior_max(i)>prior_min(i))
	  {
	    marker.push_back(i);      
	  }
    }

  Node::Node(const string& pname, const string&  pdata_dir, bool psingle_timeseries, bool pdo_load_mask) :
    single_timeseries(psingle_timeseries),
    do_load_mask(pdo_load_mask),
    name(pname),
    data_dir(pdata_dir)
  {
    setup();
  }

  void Node::setup()
  {
    Tracer_Plus trace("Node::setup");

 

  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////

  Subject_Node::Subject_Node(string& psubject_name, const Node& pnode) :
    subject_name(psubject_name),
    node(pnode),
    min_func_space_voxel_coord(3),
    max_func_space_voxel_coord(3),
    decode(false)
  {
    setup();
  }

  void Subject_Node::load_mask()
  {

   if(get_node().do_load_mask)
      {
	// load node mask 
	read_volume(std_mask,get_node().data_dir+string("/")+subject_name+string("/")+get_node().name);//+string("_mask"));

      }

    if(!get_node().single_timeseries)
      {
	// set mean, min and max for mean MVN from mask (mask is in std space)  

	ColumnVector min_standard_space_coord(3);	ColumnVector max_standard_space_coord(3);
	min_standard_space_coord=1e10;
	max_standard_space_coord=0;
	for(int x=0; x<std_mask.xsize(); x++)
	  for(int y=0; y<std_mask.ysize(); y++)
	    for(int z=0; z<std_mask.zsize(); z++)
	      {
		if(std_mask(x,y,z))
		  {
		    if(x<min_standard_space_coord(1)) min_standard_space_coord(1)=x;
		    if(x>max_standard_space_coord(1)) max_standard_space_coord(1)=x;
		    if(y<min_standard_space_coord(2)) min_standard_space_coord(2)=y;
		    if(y>max_standard_space_coord(2)) max_standard_space_coord(2)=y;
		    if(z<min_standard_space_coord(3)) min_standard_space_coord(3)=z;
		    if(z>max_standard_space_coord(3)) max_standard_space_coord(3)=z;
		  }
	      }

	mvn.update.ReSize(3);

// 	OUT("====================== Node::setup");
// 	OUT(name);

        for(int i=1; i<=3; i++)
	  {
// 	    OUT(i);
// 	    OUT(max_standard_space_coord(i));
// 	    OUT(min_standard_space_coord(i));

	    //	    if(min_standard_space_coord(i)==max_standard_space_coord(i)) mvn.update(i)=0;
	    if((max_standard_space_coord(i)-min_standard_space_coord(i))<6) mvn.update(i)=0;
	    else mvn.update(i)=1;
	  }

	//mvn.mean_standard_space_prior_mean=min_standard_space_coord+1; mvn.mean_standard_space_prior_mean(3)=0;
	mvn.mean_standard_space_prior_mean=(max_standard_space_coord+min_standard_space_coord)/2.0;

// 	OUT(mvn.mean_standard_space_prior_mean);

	mvn.mean_standard_space_prior_min=min_standard_space_coord;
	mvn.mean_standard_space_prior_max=max_standard_space_coord;

	mvn.sqrt_cov_prior_mean.ReSize(3);
	mvn.sqrt_cov_prior_min.ReSize(3);
	mvn.sqrt_cov_prior_max.ReSize(3);

	mvn.sqrt_cov_prior_mean << IdentityMatrix(3);

        for(int i=1; i<=3; i++)
	  {
	    if(mvn.update(i))
	      {
		//	      mvn.sqrt_cov_prior_mean(i,i)=2;
		mvn.sqrt_cov_prior_mean(i,i)=(max_standard_space_coord(i)-min_standard_space_coord(i));
		for(int j=i+1; j<=3; j++)
		  if(mvn.update(j))
		    mvn.sqrt_cov_prior_mean(i,j)=0.01;
	      }
	  }

	mvn.sqrt_cov_prior_min=1;
	mvn.sqrt_cov_prior_min(1,2)=-100;mvn.sqrt_cov_prior_min(1,3)=-100;mvn.sqrt_cov_prior_min(2,3)=-100;
	mvn.sqrt_cov_prior_max=100;

// 	OUT("======================");

      }
  }
  
  void Subject_Node::setup()
  {
    Tracer_Plus trace("Subject_Node::setup");

    load_mask();

    // initialise HRF params to prior means
    hrf_value=node.gethrf().prior_mean;

    // initialise balloon model values and priors

    //////////////////////////////////
    // priors taken from spm_hdm_priors.m with num modes, h=3.
//       {
//     balloon_cbf_prior_mean.ReSize(2);
//     balloon_cbf_prior_mean(1)=0.65;
//     balloon_cbf_prior_mean(2)=0.41;

//     balloon_prior_mean.ReSize(3);
//     balloon_prior_mean(1)=0.98;
//     balloon_prior_mean(2)=0.32;
//     balloon_prior_mean(3)=0.34;
    
//     balloon2_prior_mean.ReSize(1);
//     balloon2_prior_mean(1)=0;

//     balloon_cbf_prior_cov.ReSize(2);
//     balloon_cbf_prior_cov.Row(1) << 0.0148;
//     balloon_cbf_prior_cov.Row(2) << -0.0001  << 0.0018;

//     balloon_prior_cov.ReSize(3);
//     balloon_prior_cov.Row(1)  << 0.0514;
//     balloon_prior_cov.Row(2)  << 0.0021 << 0.0002;
//     balloon_prior_cov.Row(3)  << 0.0018 <<  0.0001 << 0.0001;
 
// //     balloon_prior_cov.Row(1) << 0.0150;
// //     balloon_prior_cov.Row(2) << 0.0052  << 0.0020;
// //     balloon_prior_cov.Row(3) << 0.0283  << 0.0104  << 0.0568;
// //     balloon_prior_cov.Row(4) << 0.0002  << 0.0004  << 0.0010  << 0.0013;
// //     balloon_prior_cov.Row(5) << -0.0027 << -0.0013 << -0.0069 << -0.0010 << 0.0024;

//     balloon2_prior_cov.ReSize(1);
//     balloon2_prior_cov.Row(1) << 0.0312;
//   }

    ///////////////////////////////////////////

  }  

  void Subject_Node::set_std2func_info(const Matrix& pstd2func_xform, const ColumnVector& pstdvoxdim, const ColumnVector& pfuncvoxdim, const volume<float>& pfunc_subject_mask)
  {
    Tracer_Plus trace("Subject_Node::set_std2func_info");

    std2func_xform=pstd2func_xform;
    stdvoxdim=pstdvoxdim;
    funcvoxdim=pfuncvoxdim;

    // mask is in std space
    // coordinates need to be in functional space
    // so affine xform mask into func space and then extract coords

    // setup container for node mask in functional space
    func_mask=pfunc_subject_mask;
    
    if((get_std_mask().xdim() != stdvoxdim(1)) ||(get_std_mask().ydim() != stdvoxdim(2)) || (get_std_mask().zdim() != stdvoxdim(3)))
      {
	LOGOUT(get_std_mask().xdim());
	LOGOUT(get_std_mask().xsize());

	LOGOUT(string("Voxel dimensions for std mask of ")+node.name+string(" != voxel dimensions of passed in standard brain"));
	throw Exception(string(string("Voxel dimensions for std mask of ")+node.name+string(" != voxel dimensions of passed in standard brain")).data());	
      }

    affine_transform(get_std_mask(),func_mask,std2func_xform);

    save_volume(func_mask, LogSingleton::getInstance().appendDir(node.get_name()+"_func_mask"));	  
    save_volume(get_std_mask(), LogSingleton::getInstance().appendDir(node.get_name()+"_std_node_mask"));

 //    volume<float> std_test=get_std_mask();

//     affine_transform(func_mask,std_test,std2func_xform.i());
 
//     save_volume(get_std_mask(), LogSingleton::getInstance().appendDir(node.get_name()+"_test_std_node_mask"));
    
    if(!node.single_timeseries)
      {
	// initialise MVN params to prior means
	mvn_cov_value << IdentityMatrix(3);    
	
	mvn_sqrt_cov_value=mvn.sqrt_cov_prior_mean;
	
	mvn_mean_standard_space_value=ColumnVector2vector(mvn.mean_standard_space_prior_mean);


	// take square
	mvn_cov_value << mvn_sqrt_cov_value*mvn_sqrt_cov_value;
	
	// mean is passed in as a vector of values in standard space
	// convert mean from std space to functional space for subject:
	//	mvn_mean_func_space_value=xform_coords(mvn_mean_standard_space_value,std2func_xform,stdvoxdim,funcvoxdim);
	ColumnVector tmp=vector2ColumnVector(mvn_mean_standard_space_value);
	mvn_mean_func_space_value=NewimageCoord2NewimageCoord(std2func_xform,get_std_mask(),func_mask,tmp);

// 	OUT("====================== Subject_Node::set_std2func_info");

// 	OUT(node.name);
// 	OUT(mvn_mean_standard_space_value);	
// 	OUT(mvn_mean_func_space_value);
// 	OUT(funcvoxdim);
// 	OUT(stdvoxdim);

// 	// do sanity check
// 	//ColumnVector sanity_mean_std_space=xform_coords(mvn_mean_func_space_value,std2func_xform.i(),funcvoxdim,stdvoxdim);
// 	ColumnVector sanity_mean_std_space=NewimageCoord2NewimageCoord(std2func_xform.i(),func_mask,get_std_mask(),mvn_mean_func_space_value);
// 	OUT(sanity_mean_std_space);
       
// 	OUT("======================");
	
	// calculate normalisation factor (in func space)
	if(decode)
	  mvn_norm=1.0/mvnpdf(mvn_mean_func_space_value.t(), mvn_mean_func_space_value.t(), mvn_cov_value);
	else
	  mvn_norm=1;

	try
	  {		
	    covar_inv=mvn_cov_value.i();
	  }
	catch(Exception& exc)
	  {

	    LOGOUT("Warning: invalid mvn_cov_value in Subject_Node::set_std2func_info");
	    LOGOUT(mvn_cov_value);
	    mvn_norm=1e20;
	  }
      }

   
  }


  /////

  void Subject_Node::set_value_hrf(const vector<float>& hrf_vals) // this sets the hrf values for those being estimated
  {
    for(unsigned int i=1; i <= node.gethrf().marker.size(); i++)
      hrf_value(node.gethrf().marker[i-1])=hrf_vals[i-1];
  }

  const vector<float> Subject_Node::get_value_hrf() const // this returns the hrf values for those being estimated
  {
    vector<float> ret(node.gethrf().marker.size());

    for(unsigned int i=1; i <= ret.size(); i++)
      ret[i-1]=hrf_value(node.gethrf().marker[i-1]);

    return ret;
  }

  const vector<string> Subject_Node::get_names_hrf() const // this returns the hrf names for those being estimated
  {
    vector<string> ret(node.gethrf().marker.size());

    for(unsigned int i=1; i <= ret.size(); i++)
      ret[i-1]=string(node.get_name()+"_hrf_"+num2str(node.gethrf().marker[i-1]));

    return ret;
  }
  
  const vector<float> Subject_Node::get_mins_hrf() const // this returns the hrf mins for those being estimated
  {
    vector<float> ret(node.gethrf().marker.size());

    for(unsigned int i=1; i <= ret.size(); i++)
      ret[i-1]=node.gethrf().prior_min(node.gethrf().marker[i-1]);

    return ret;
  }
 
  const vector<float> Subject_Node::get_maxs_hrf() const // this returns the hrf maxs for those being estimated
  {
    vector<float> ret(node.gethrf().marker.size());

    for(unsigned int i=1; i <= ret.size(); i++)
      ret[i-1]=node.gethrf().prior_max(node.gethrf().marker[i-1]);

    return ret;
  }

  /////

  void Subject_Node::set_value_balloon_cbf(const vector<float>& balloon_vals) // this sets the balloon values for those being estimated
  {
     kappa_value=balloon_vals[0];
     gamma_value=balloon_vals[1];
  }

  void Subject_Node::set_value_balloon(const vector<float>& balloon_vals) // this sets the balloon values for those being estimated
  {
     tau_value=balloon_vals[0];
     alpha_value=balloon_vals[1];
     E0_value=balloon_vals[2];
  }

  void Subject_Node::set_value_balloon2(const vector<float>& balloon_vals) // this sets the balloon values for those being estimated
  {
     logepsilon_value=balloon_vals[0];
  }

  const vector<float> Subject_Node::get_value_balloon_cbf() const // this returns the balloon values for those being estimated
  {
    vector<float> ret;

    ret.push_back(kappa_value);
    ret.push_back(gamma_value);

    return ret;
  }

  const vector<float> Subject_Node::get_value_balloon() const // this returns the balloon values for those being estimated
  {
    vector<float> ret;

    ret.push_back(tau_value);
    ret.push_back(alpha_value);
    ret.push_back(E0_value);

    return ret;
  }

  const vector<float> Subject_Node::get_value_balloon2() const // this returns the balloon values for those being estimated
  {
    vector<float> ret;

    ret.push_back(logepsilon_value);

    return ret;
  }


  const vector<string> Subject_Node::get_names_balloon_cbf() const // this returns the balloon names for those being estimated
  {
    vector<string> ret;

    ret.push_back(string(node.get_name())+"_kappa");
    ret.push_back(string(node.get_name())+"_gamma");
    return ret;
  }

  const vector<string> Subject_Node::get_names_balloon() const // this returns the balloon names for those being estimated
  {
    vector<string> ret;

    ret.push_back(string(node.get_name())+"_tau");
    ret.push_back(string(node.get_name())+"_alpha");
    ret.push_back(string(node.get_name())+"_E0");
    return ret;
  }
  
  const vector<string> Subject_Node::get_names_balloon2() const // this returns the balloon names for those being estimated
  {
    vector<string> ret;

    ret.push_back(string(node.get_name())+"_logepsilon");

    return ret;
  }

  const vector<float> Subject_Node::get_mins_balloon_cbf() const // this returns the balloon mins for those being estimated
  {
    vector<float> ret(balloon_cbf_prior_mean.Nrows());

    for(unsigned int i=1; i<=ret.size(); i++)
      ret[i-1]=balloon_cbf_prior_mean(i)-10*sqrt(balloon_cbf_prior_cov(i,i));

    return ret;
  }

  const vector<float> Subject_Node::get_mins_balloon() const // this returns the balloon mins for those being estimated
  {
    vector<float> ret(balloon_prior_mean.Nrows());

    for(unsigned int i=1; i<=ret.size(); i++)
      ret[i-1]=balloon_prior_mean(i)-10*sqrt(balloon_prior_cov(i,i));

    return ret;
  }
 
  const vector<float> Subject_Node::get_mins_balloon2() const // this returns the balloon mins for those being estimated
  {
    vector<float> ret(balloon2_prior_mean.Nrows());

    for(unsigned int i=1; i<=ret.size(); i++)
      ret[i-1]=balloon2_prior_mean(i)-10*sqrt(balloon2_prior_cov(i,i));

    return ret;
  }
 
  const vector<float> Subject_Node::get_maxs_balloon_cbf() const // this returns the balloon maxs for those being estimated
  {
    vector<float> ret(balloon_cbf_prior_mean.Nrows());

    for(unsigned int i=1; i<=ret.size(); i++)
      ret[i-1]=balloon_cbf_prior_mean(i)+10*sqrt(balloon_cbf_prior_cov(i,i));

    return ret;
  }

  const vector<float> Subject_Node::get_maxs_balloon() const // this returns the balloon maxs for those being estimated
  {
    vector<float> ret(balloon_prior_mean.Nrows());

    for(unsigned int i=1; i<=ret.size(); i++)
      ret[i-1]=balloon_prior_mean(i)+10*sqrt(balloon_prior_cov(i,i));

    return ret;
  }

  const vector<float> Subject_Node::get_maxs_balloon2() const // this returns the balloon maxs for those being estimated
  {
    vector<float> ret(balloon2_prior_mean.Nrows());

    for(unsigned int i=1; i<=ret.size(); i++)
      ret[i-1]=balloon2_prior_mean(i)+10*sqrt(balloon2_prior_cov(i,i));

    return ret;
  }

  /////

  void Subject_Node::set_value_mvn_mean_standard_space(const vector<float>& pmean_standard_space_value) // this sets the mvn_mean_standard_space values for those being estimated
  {
    //    LOGOUT(mvn_mean_standard_space_value);

    int index=0;
    for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))
	mvn_mean_standard_space_value[i-1]=pmean_standard_space_value[index++];

//     LOGOUT("++++++++++");
//     LOGOUT("In Subject_Node::set_value_mvn_mean_standard_space");
//     LOGOUT(node.name);
//     LOGOUT(pmean_standard_space_value);
//     LOGOUT(std2func_xform);
//     LOGOUT(funcvoxdim);
//     LOGOUT(stdvoxdim);

    // mean is passed in as a vector of values in standard space
    // convert mean from std space to functional space for subject:
    //    mvn_mean_func_space_value=xform_coords(mvn_mean_standard_space_value,std2func_xform,stdvoxdim,funcvoxdim);
    mvn_mean_func_space_value=NewimageCoord2NewimageCoord(std2func_xform,get_std_mask(),func_mask,vector2ColumnVector(mvn_mean_standard_space_value));

//     LOGOUT(mvn_mean_func_space_value);
//     LOGOUT("++++++++++");

    // calculate normalisation factor (in func space)
    if(decode)
      mvn_norm=1.0/mvnpdf(mvn_mean_func_space_value.t(), mvn_mean_func_space_value.t(), mvn_cov_value);
    else
      mvn_norm=1;

    try
      {	
	covar_inv=mvn_cov_value.i();
      }
    catch(Exception& exc)
      {
	
	LOGOUT("Warning: invalid mvn_cov_value in Subject_Node::set_value_mvn_mean_standard_space");
	LOGOUT(mvn_cov_value);
	mvn_norm=1e20;
      }
  }

  void Subject_Node::set_full_value_mvn_mean_func_space(const ColumnVector& pmean_func_space_value) // this sets the mvn_mean_func_space values 
  {
    mvn_mean_func_space_value=pmean_func_space_value;

    LOGOUT("************");
    LOGOUT("In Subject_Node::set_full_value_mvn_mean_func_space");
//     LOGOUT(node.name);
    LOGOUT(mvn_mean_func_space_value);
//     LOGOUT(std2func_xform);
//     LOGOUT(funcvoxdim);
//     LOGOUT(stdvoxdim);

    const volume<float>& std_mask=get_std_mask();

    //  mvn_mean_standard_space_value=ColumnVector2vector(xform_coords(ColumnVector2vector(mvn_mean_func_space_value),std2func_xform.i(),funcvoxdim,stdvoxdim));
    ColumnVector mvn_mean_standard_space_value_col=NewimageCoord2NewimageCoord(std2func_xform.i(),func_mask,std_mask,mvn_mean_func_space_value);
    
    LOGOUT(mvn_mean_standard_space_value);
    
    // check if centre is in std space mask:
    bool inbounds=(std_mask(::round(mvn_mean_standard_space_value_col(1)),::round(mvn_mean_standard_space_value_col(2)),::round(mvn_mean_standard_space_value_col(3)))>0);
    
    if(!inbounds)
      {
	// find the nearest voxel to mvn_mean_standard_space_value that is in the std space mask
	float min_dist=1e32;
	ColumnVector nearest_coord(3);
	ColumnVector coord(3);

	for(int x=mvn.mean_standard_space_prior_min(1); x<=mvn.mean_standard_space_prior_max(1); x++)
	  for(int y=mvn.mean_standard_space_prior_min(2); y<=mvn.mean_standard_space_prior_max(2); y++)
	    for(int z=mvn.mean_standard_space_prior_min(3); z<=mvn.mean_standard_space_prior_max(3); z++)
	      {
		if(std_mask(x,y,z))
		  {
		    coord(1)=x; coord(2)=y; coord(3)=z;
		    float dist=(coord-mvn_mean_standard_space_value_col).SumSquare();
		    if(dist < min_dist)
		      {
			min_dist=dist;
			nearest_coord=coord;
		      }
		  }
	      }

	mvn_mean_standard_space_value_col=nearest_coord;

	// sanity check to see if centre is now in std space mask:
	inbounds=(std_mask(::round(mvn_mean_standard_space_value_col(1)),::round(mvn_mean_standard_space_value_col(2)),::round(mvn_mean_standard_space_value_col(3)))>0);

	if(!inbounds)
	  {
	    LOGOUT("Still not inbounds in Subject_Node::set_full_value_mvn_mean_func_space");
	    throw Exception("Still not inbounds in Subject_Node::set_full_value_mvn_mean_func_space");
	  }
    
      }
    
    mvn_mean_standard_space_value=ColumnVector2vector(mvn_mean_standard_space_value_col);
    LOGOUT(mvn_mean_standard_space_value);
    
    //  mvn_mean_func_space_value=xform_coords(mvn_mean_standard_space_value,std2func_xform,stdvoxdim,funcvoxdim);
    mvn_mean_func_space_value=NewimageCoord2NewimageCoord(std2func_xform,get_std_mask(),func_mask,vector2ColumnVector(mvn_mean_standard_space_value));
    
    LOGOUT(mvn_mean_func_space_value);

    LOGOUT("************");    

    // calculate normalisation factor (in func space)
     if(decode)
      mvn_norm=1.0/mvnpdf(mvn_mean_func_space_value.t(), mvn_mean_func_space_value.t(), mvn_cov_value);
    else
      mvn_norm=1;

     try
       {	
	 covar_inv=mvn_cov_value.i();
       }
     catch(Exception& exc)
       {

	 LOGOUT("Warning: invalid mvn_cov_value in Subject_Node::set_full_value_mvn_mean_func_space");
	 LOGOUT(mvn_cov_value);
	 mvn_norm=1e20;
       }	
  }

  const vector<float> Subject_Node::get_value_mvn_mean_standard_space() const // this returns the mvn_mean_standard_space values for those being estimated
  {
    vector<float> ret;
    for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))
	ret.push_back(mvn_mean_standard_space_value[i-1]);

    //    LOGOUT(ret);

    return ret;
  }

  const vector<string> Subject_Node::get_names_mvn_mean_standard_space() const // this returns the mvn_mean_standard_space names for those being estimated
  {
    vector<string> ret;

    if(!node.single_timeseries)
      for(unsigned int i=1; i <= 3; i++)
	if(mvn.update(i))
	  ret.push_back(string(node.get_name()+"_mvn_mean_standard_space_"+num2str(i)));
    
    return ret;
  }
  
  const vector<float> Subject_Node::get_mins_mvn_mean_standard_space() const // this returns the mvn_mean_standard_space mins for those being estimated
  {
    vector<float> ret;

    for(unsigned int i=1; i <= 3; i++)
      if(mvn.update(i))
	ret.push_back(mvn.mean_standard_space_prior_min(i));

    return ret;
  }
 
  const vector<float> Subject_Node::get_maxs_mvn_mean_standard_space() const // this returns the mvn_mean_standard_space maxs for those being estimated
  {
    vector<float> ret;

    for(unsigned int i=1; i <= 3; i++)
      if(mvn.update(i))
	ret.push_back(mvn.mean_standard_space_prior_max(i));
 
    return ret;
  }

  //////////

  void Subject_Node::set_full_value_mvn_sqrt_cov(const SymmetricMatrix& psqrt_cov_value) // this sets the mvn_sqrt_cov values
  {
 
    // cov is passed in as a vector of values, these values represent the square root of the 3D spatial covariance matrix in functional space
    
//     OUT("Subject_Node::set_full_value_mvn_sqrt_cov");
//     OUT(psqrt_cov_value);

    // convert to matrix
    mvn_sqrt_cov_value << psqrt_cov_value;

    // take square
    mvn_cov_value << mvn_sqrt_cov_value*mvn_sqrt_cov_value;
 
    // calculate normalisation factor (in func space)
    if(decode)
      mvn_norm=1.0/mvnpdf(mvn_mean_func_space_value.t(), mvn_mean_func_space_value.t(), mvn_cov_value);
    else
      mvn_norm=1;

    try
      {	
	covar_inv=mvn_cov_value.i();
      }
    catch(Exception& exc)
      {

	LOGOUT("Warning: invalid mvn_cov_value in Subject_Node::set_full_value_mvn_sqrt_cov");
	LOGOUT(mvn_cov_value);
	mvn_norm=1e20;
      }
  }

  //////////

  void Subject_Node::set_value_mvn_sqrt_cov(const vector<float>& psqrt_cov_value) // this sets the mvn_sqrt_cov values for those being estimated
  {
 
    // cov is passed in as a vector of values, these values represent the square root of the 3D spatial covariance matrix in functional space
    
    // convert to matrix
    int index=0;
    for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))
	mvn_sqrt_cov_value(i,i)=psqrt_cov_value[index++];

    // take square
    mvn_cov_value << mvn_sqrt_cov_value*mvn_sqrt_cov_value;

    // calculate normalisation factor (in func space)
    if(decode)
      mvn_norm=1.0/mvnpdf(mvn_mean_func_space_value.t(), mvn_mean_func_space_value.t(), mvn_cov_value);
    else
      mvn_norm=1;

    try
      {	
	covar_inv=mvn_cov_value.i();
      }
    catch(Exception& exc)
      {


	LOGOUT("Warning: invalid mvn_cov_value in Subject_Node::set_value_mvn_sqrt_cov");
	LOGOUT(mvn_cov_value);
	mvn_norm=1e20;
      }
  }  

  const vector<float> Subject_Node::get_value_mvn_sqrt_cov() const // this returns the mvn_sqrt_cov values for those being estimated
  {
    vector<float> ret;
    
    for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))	
	ret.push_back(mvn_sqrt_cov_value(i,i));
    
    return ret;
  }

  const vector<string> Subject_Node::get_names_mvn_sqrt_cov() const // this returns the mvn_sqrt_cov names for those being estimated
  {
    vector<string> ret;

    if(!node.single_timeseries)
      for(unsigned int i=1; i<=3; i++)
	if(mvn.update(i))
	  ret.push_back(string(node.get_name()+"_mvn_sqrt_cov_"+num2str(i)));
    
    return ret;
  }
  
  const vector<float> Subject_Node::get_mins_mvn_sqrt_cov() const // this returns the mvn_sqrt_cov mins for those being estimated
  {
    vector<float> ret;

    for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))
	ret.push_back(mvn.sqrt_cov_prior_min(i,i));

    return ret;
  }
 

  const vector<float> Subject_Node::get_maxs_mvn_sqrt_cov() const // this returns the mvn_sqrt_cov maxs for those being estimated
  {
    vector<float> ret;

     for(unsigned int i=1; i<=3; i++)
       if(mvn.update(i))
	 ret.push_back(mvn.sqrt_cov_prior_max(i,i));
 
    return ret;
  }

  //////////

  void Subject_Node::set_value_mvn_sqrt_cov_offdiag(const vector<float>& psqrt_cov_offdiag_value) // this sets the mvn_sqrt_cov_offdiag values for those being estimated
  {
 
    // cov is passed in as a vector of values, these values represent the square root of the 3D spatial covariance matrix in functional space
    
    // convert to matrix
    int index=0;
    for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))
	for(unsigned int j=i+1; j<=3; j++)
	  //	  for(unsigned int j=i; j<=i; j++)
	    if(mvn.update(j))
	      mvn_sqrt_cov_value(i,j)=psqrt_cov_offdiag_value[index++];

    // take square
    mvn_cov_value << mvn_sqrt_cov_value*mvn_sqrt_cov_value;

    // calculate normalisation factor (in func space)
    if(decode)
      mvn_norm=1.0/mvnpdf(mvn_mean_func_space_value.t(), mvn_mean_func_space_value.t(), mvn_cov_value);
    else
      mvn_norm=1;

    try
      {
	covar_inv=mvn_cov_value.i();
      }
    catch(Exception& exc)
      {

	LOGOUT("Warning: invalid mvn_cov_value in Subject_Node::set_value_mvn_sqrt_cov_offdiag");
	LOGOUT(mvn_cov_value);
	mvn_norm=1e10;
	
	mvn_cov_value=0;
	for(int p=1; p<=mvn_cov_value.Nrows(); p++)
	  mvn_cov_value(p,p)=1;
	covar_inv=mvn_cov_value;
      }
  }  

  const vector<float> Subject_Node::get_value_mvn_sqrt_cov_offdiag() const // this returns the mvn_sqrt_cov_offdiag values for those being estimated
  {
    vector<float> ret;
    
    for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))
	for(unsigned int j=i+1; j<=3; j++)
	  if(mvn.update(j))
	    ret.push_back(mvn_sqrt_cov_value(i,j));
    
    return ret;
  }

  const vector<string> Subject_Node::get_names_mvn_sqrt_cov_offdiag() const // this returns the mvn_sqrt_cov_offdiag names for those being estimated
  {
    vector<string> ret;

    if(!node.single_timeseries)
      for(unsigned int i=1; i<=3; i++)
	if(mvn.update(i))
	  for(unsigned int j=i+1; j<=3; j++)
	    if(mvn.update(j))
	      ret.push_back(string(node.get_name()+"_mvn_sqrt_cov_"+num2str(i)+"_"+num2str(j)));
    
    return ret;
  }
  
  const vector<float> Subject_Node::get_mins_mvn_sqrt_cov_offdiag() const // this returns the mvn_sqrt_cov_offdiag mins for those being estimated
  {
    vector<float> ret;

    for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))
	for(unsigned int j=i+1; j<=3; j++)
	  if(mvn.update(j))
	    ret.push_back(mvn.sqrt_cov_prior_min(i,j));

    return ret;
  }
 

  const vector<float> Subject_Node::get_maxs_mvn_sqrt_cov_offdiag() const // this returns the mvn_sqrt_cov_offdiag maxs for those being estimated
  {
    vector<float> ret;

     for(unsigned int i=1; i<=3; i++)
      if(mvn.update(i))
	for(unsigned int j=i+1; j<=3; j++)
	  if(mvn.update(j))
	    ret.push_back(mvn.sqrt_cov_prior_max(i,j));
 
    return ret;
  }

  ////////////////

  void Subject_Node::copy_hrf_values(const Subject_Node& sub_node)
  {
    hrf_value=sub_node.hrf_value;
    kappa_value=sub_node.kappa_value;
    gamma_value=sub_node.gamma_value;
    tau_value=sub_node.tau_value;
    alpha_value=sub_node.alpha_value;
    E0_value=sub_node.E0_value;
    logepsilon_value=sub_node.logepsilon_value;
  }

  void Subject_Node::copy_roi_values(const Subject_Node& sub_node)
  {
    Tracer_Plus trace("Subject_Node::copy_roi_values");
    
    mvn_norm=sub_node.mvn_norm;
    covar_inv=sub_node.covar_inv;
    mvn_mean_func_space_value=sub_node.mvn_mean_func_space_value;
    mvn_mean_standard_space_value=sub_node.mvn_mean_standard_space_value;
    mvn_cov_value=sub_node.mvn_cov_value;
    mvn_sqrt_cov_value=sub_node.mvn_sqrt_cov_value;

    //    OUT("==============");
//     OUT(mvn_norm);
//     OUT(mvn_sqrt_cov_value);
//    OUT(mvn_mean_func_space_value);
//     OUT(mvn_mean_standard_space_value);
//    OUT("==============");

  }

  ///////////////////
  
  void Subject_Node::setup_mvn_basis_set()
  {
    Tracer_Plus trace("Subject_Node::setup_mvn_basis_set");
    
    ColumnVector diff=max_func_space_voxel_coord-min_func_space_voxel_coord;
    
    OUT("====================== Subject_Node::setup_mvn_basis_set");
    OUT(node.name);
    OUT(diff);

    basis_nums.clear();
    basis_res.clear();
    basis_centres.clear();
    
    int total_num=9;   

    // calc num in each dimension
    for(int n=1; n<=3; n++)
      {
	OUT(mvn.update(n));

	if(mvn.update(n))// && diff(n)>4)
	  {
	    // calculate number of centres needed in this dimension based on wanting approx total_num
	    // int numtmp=int(Max(total_num*std::sqrt(diff(n)/Sum(diff)),1));
	    // x/y=X/Y; x/z=X/Z; xyz=N; => x^3=NX^2/YZ
	    float YZ=1;float X2=1;
	    int count=1;
	    for(int n2=1; n2<=3; n2++)
	      {
		if(mvn.update(n2) && n!=n2)
		  {
		    count++;
		    YZ*=diff(n2);
		    X2*=diff(n);
		  }
	      }

	    int numtmp=int(std::floor(0.5+Max(std::pow(total_num*X2/YZ,float(1.0/count)),1.0)));
	    basis_nums.push_back(numtmp);
	    basis_res.push_back(((max_func_space_voxel_coord(n)-min_func_space_voxel_coord(n))/float(numtmp)));

// 	    OUT(YZ);
// 	    OUT(X2);
// 	    OUT(total_num);
//  	    OUT(count);
// 	    OUT(std::pow(total_num*X2/YZ,float(1.0/count)));
// 	    OUT(Max(std::pow(total_num*X2/YZ,float(1.0/count)),1.0));
// 	    OUT(std::floor(0.5+Max(std::pow(total_num*X2/YZ,float(1.0/count)),1.0)));
// 	    OUT(numtmp);

// 	    OUT(max_func_space_voxel_coord(n));
// 	    OUT(min_func_space_voxel_coord(n));

	    OUT(basis_nums[n-1]);
	    OUT(basis_res[n-1]);

	    ColumnVector centrestmp(numtmp); 
	    for(int m=1; m<=numtmp; m++)
	      {
		centrestmp(m)=min_func_space_voxel_coord(n)+basis_res[n-1]/2.0+(m-1)*basis_res[n-1];
	      }

	    basis_centres.push_back(centrestmp);	

	    OUT(centrestmp);	  
	  }
	else// if(mvn.update(n))
	  {
	    int numtmp=1;
	    basis_nums.push_back(numtmp);
	    basis_res.push_back(3);

	    OUT(max_func_space_voxel_coord(n));
	    OUT(min_func_space_voxel_coord(n));

	    ColumnVector centrestmp(numtmp); 
	    centrestmp(1)=0.5*(max_func_space_voxel_coord(n)+min_func_space_voxel_coord(n));	      

	    basis_centres.push_back(centrestmp);

	    OUT(centrestmp);
	  }
// 	else
// 	  {
// 	    int numtmp=1;
// 	    basis_nums.push_back(numtmp);
// 	    basis_res.push_back(3);

// 	    ColumnVector centrestmp(numtmp); 
// 	    centrestmp(1)=min_func_space_voxel_coord(n);	      

// 	    basis_centres.push_back(centrestmp);		  

// 	OUT(centrestmp);
// 	  }

      }
    
    basis_centres_list.clear();

    // setup list containing each basis fn centre (in func space)
    for(int i=1; i<=basis_nums[0]; i++)
      for(int j=1; j<=basis_nums[1]; j++)
	for(int k=1; k<=basis_nums[2]; k++)		
	  {
	    ColumnVector coords(3);
	    coords(1)=basis_centres[0](i);
	    coords(2)=basis_centres[1](j);
	    coords(3)=basis_centres[2](k);

	    basis_centres_list.push_back(coords);
	  }
    
    num_basis=basis_centres_list.size();

    OUT(num_basis);

    OUT("======================");

  }
  
  void Subject_Node::set_mvn_basis(int b)
  {    
    Tracer_Plus trace("Subject_Node::set_mvn_basis");

//     ColumnVector centre_std_coords=xform_coords(basis_centres_list[b-1],std2func_xform.i(),funcvoxdim,stdvoxdim);

    OUT(basis_centres_list[b-1]);

    set_full_value_mvn_mean_func_space(basis_centres_list[b-1]);

    SymmetricMatrix sqrt_cov_value;
    sqrt_cov_value << IdentityMatrix(3);
    sqrt_cov_value(1,1)=basis_res[0]/2.0;
    sqrt_cov_value(2,2)=basis_res[1]/2.0;
    sqrt_cov_value(3,3)=basis_res[2]/2.0;

    OUT(sqrt_cov_value);

    set_full_value_mvn_sqrt_cov(sqrt_cov_value);
  }

}
