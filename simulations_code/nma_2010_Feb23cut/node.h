/*  node.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(node_h)
#define node_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "nmaoptions.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

namespace Nma {
  
//   ReturnMatrix xform_coords(const ColumnVector& in, const Matrix& xform, const ColumnVector& in_voxdim, const ColumnVector& out_voxdim);

//   ReturnMatrix xform_coords(const vector<float>& in, const Matrix& xform, const ColumnVector& in_voxdim, const ColumnVector& out_voxdim);

  class Hrf
  {
  public:
    
    Hrf() {}
    Hrf(const ColumnVector& pprior_mean, const ColumnVector& pprior_min, const ColumnVector& pprior_max);

    void set(const ColumnVector& pprior_mean, const ColumnVector& pprior_min, const ColumnVector& pprior_max);   
      
    ColumnVector prior_mean; // m1,m2,m3,m4,c1,c2 for halfcos hrf
    ColumnVector prior_min; 
    ColumnVector prior_max;

    vector<int> marker;
  };
  
  class Mvn
  {
  public:
    Mvn(){}

    Mvn(const ColumnVector& pmean_standard_space_prior_mean, const ColumnVector& pmean_standard_space_prior_min, const ColumnVector& pmean_standard_space_prior_max, const SymmetricMatrix& psqrt_cov_prior_mean, const SymmetricMatrix& psqrt_cov_prior_min, const SymmetricMatrix& psqrt_cov_prior_max) :
      mean_standard_space_prior_mean(pmean_standard_space_prior_mean),
      mean_standard_space_prior_min(pmean_standard_space_prior_min),
      mean_standard_space_prior_max(pmean_standard_space_prior_max),
      sqrt_cov_prior_mean(psqrt_cov_prior_mean),
      sqrt_cov_prior_min(psqrt_cov_prior_min),
      sqrt_cov_prior_max(psqrt_cov_prior_max)
    {
    }

    ColumnVector update; // update(i)=0 if no estimation needed in dimension i

    // these are passed-in (estimated) in std space
    ColumnVector mean_standard_space_prior_mean;
    ColumnVector mean_standard_space_prior_min;
    ColumnVector mean_standard_space_prior_max;

    // covs passed-in (estimated) as the sqrt of the cov matrix in functional space 
    SymmetricMatrix sqrt_cov_prior_mean;
    SymmetricMatrix sqrt_cov_prior_min;
    SymmetricMatrix sqrt_cov_prior_max;

  };

  class Subject_Node;

  class Node
  {
  public:

    Node(const string& pname, const string&  pdata_dir, bool psingle_timeseries, bool do_load_mask);  

    void setup();

    void set_hrf_prior(const ColumnVector& hrf_prior_mean,const ColumnVector& hrf_prior_min,const ColumnVector& hrf_prior_max) {
      hrf.set(hrf_prior_mean,hrf_prior_min,hrf_prior_max);
    }

    const Hrf& gethrf() const {return hrf;}

    const string& get_name() const {return name;}
    

    // Destructor
    virtual ~Node() {}
 
  private:

    // constructor
    Node();
    const Node& operator=(Node& par);     
    Node(Node& des);      

    // data members
    bool single_timeseries;
    bool do_load_mask;

    string name; 
    string data_dir;
 
    Hrf hrf;

    friend class Subject_Node;

  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  class Subject_Node
  {
  public:
    
    Subject_Node(string& psubject_name, const Node& pnode);
    
    void setup();
    void load_mask();

    void set_std2func_info(const Matrix& pstd2func_xform, const ColumnVector& pstdvoxdim, const ColumnVector& pfuncvoxdim, const volume<float>& pfunc_mask);

    void setup_mvn_basis_set();
    void set_mvn_basis(int b);
    int get_num_basis() const {return num_basis;}
    const Mvn& getmvn() const {return mvn;}

    const volume<float>& get_std_mask() const {return std_mask;}  // in std space

    const ColumnVector& get_min_func_space_voxel_coordinate() const {return min_func_space_voxel_coord;}
    const ColumnVector& get_max_func_space_voxel_coordinate() const {return max_func_space_voxel_coord;}

    void set_range_func_space_voxel_coordinate(const ColumnVector& min_in, const ColumnVector& max_in, bool is_single_timeseries) { 
      min_func_space_voxel_coord=min_in; 
      max_func_space_voxel_coord=max_in;

      if(!is_single_timeseries)
	{
	  ColumnVector factor=max_in-min_in;

	  for(int r=1; r<=3; r++)
	    {
	      if(factor(r)>5)
		factor(r)=(factor(r))/2.0;     
	  
	      mvn.sqrt_cov_prior_max(r,r)=factor(r);
	    }

	  float tmpfactor=Maximum(factor);
	  mvn.sqrt_cov_prior_min(1,2)=-tmpfactor;mvn.sqrt_cov_prior_min(1,3)=-tmpfactor;mvn.sqrt_cov_prior_min(2,3)=-tmpfactor;
	  mvn.sqrt_cov_prior_max(1,2)=tmpfactor;mvn.sqrt_cov_prior_max(1,3)=tmpfactor;mvn.sqrt_cov_prior_max(2,3)=tmpfactor;
	}
	
    }


    void set_balloon_prior(const ColumnVector& pballoon_cbf_means,const SymmetricMatrix& pballoon_cbf_covs,const ColumnVector& pballoon_means,const SymmetricMatrix& pballoon_covs,const ColumnVector& pballoon2_means,const SymmetricMatrix& pballoon2_covs) {
      balloon_cbf_prior_mean=pballoon_cbf_means;
      balloon_cbf_prior_cov=pballoon_cbf_covs;
      balloon_prior_mean=pballoon_means;
      balloon_prior_cov=pballoon_covs;
      balloon2_prior_mean=pballoon2_means;
      balloon2_prior_cov=pballoon2_covs; 

//       OUT(balloon_cbf_prior_mean);
//       OUT(balloon_cbf_prior_cov);
//       OUT(balloon_prior_mean);
//       OUT(balloon_prior_cov);
//       OUT(balloon2_prior_mean);
//       OUT(balloon2_prior_cov);
//       exit(1);

    balloon_cbf_prior_cov_inv=balloon_cbf_prior_cov.i();
    balloon_prior_cov_inv=balloon_prior_cov.i();
    balloon2_prior_cov_inv=balloon_prior_cov.i();

    kappa_value=balloon_cbf_prior_mean(1);
    gamma_value=balloon_cbf_prior_mean(2);

    tau_value=balloon_prior_mean(1);
    alpha_value=balloon_prior_mean(2);
    E0_value=balloon_prior_mean(3);

    logepsilon_value=balloon2_prior_mean(1);
   }

//     void set_max_func_space_voxel_coordinate(const ColumnVector& in) { max_func_space_voxel_coord=in; }

    float get_balloon_cbf_prior_energy(const ColumnVector& values)
    {
      ColumnVector tmp=values-balloon_cbf_prior_mean;

//       OUT(tmp);
//       OUT(mvn_norm);
//       OUT(covar_inv);

      return 0.5*(tmp.t()*balloon_cbf_prior_cov_inv*tmp).AsScalar();
    }

    float get_balloon_prior_energy(const ColumnVector& values)
    {
      ColumnVector tmp=values-balloon_prior_mean;

//       OUT(tmp);
//       OUT(mvn_norm);
//       OUT(covar_inv);

      return 0.5*(tmp.t()*balloon_prior_cov_inv*tmp).AsScalar();
    }

    float get_balloon_prior_energy2(const ColumnVector& values)
    {
      ColumnVector tmp=values-balloon2_prior_mean;

//       OUT(tmp);
//       OUT(mvn_norm);
//       OUT(covar_inv);

      return 0.5*(tmp.t()*balloon2_prior_cov_inv*tmp).AsScalar();
    }

    float get_partial_vol_fraction(const ColumnVector& voxel_coordinate)
    {
//       LOGOUT(voxel_index);
//       LOGOUT(mvn_norm);
//       LOGOUT(voxel_coordinates.size());
//       LOGOUT(voxel_coordinates[voxel_index-1].t());
//       LOGOUT(mvn_mean_func_space_value.t());
      //return mvnpdf(voxel_coordinates[voxel_index-1].t(), mvn_mean_func_space_value.t(), mvn_cov_value)/mvn_norm;

//       OUT(voxel_coordinates[voxel_index-1]);
//       OUT(mvn_mean_func_space_value);

      ColumnVector tmp=voxel_coordinate-mvn_mean_func_space_value;

//       OUT(tmp);
//       OUT(mvn_norm);
//       OUT(covar_inv);


//       float tmp2= std::exp(-0.5*(tmp.t()*covar_inv*tmp).AsScalar())/mvn_norm;
//       if(tmp2>0.25) tmp2=1; 
//       else tmp2=0;

//       float tmp2=1;

//       return tmp2;

      return std::exp(-0.5*(tmp.t()*covar_inv*tmp).AsScalar())/mvn_norm;
    }

    float print_partial_vol_fraction(const ColumnVector& voxel_coordinate)
    {

      OUT(voxel_coordinate);
      OUT(mvn_mean_func_space_value);

      ColumnVector tmp=voxel_coordinate-mvn_mean_func_space_value;

      OUT(tmp);
      OUT(mvn_norm);
      OUT(covar_inv);

      float ret=std::exp(-0.5*(tmp.t()*covar_inv*tmp).AsScalar())/mvn_norm;
      OUT(ret);
      return ret;
    }

    const ColumnVector& get_mvn_mean_func_space_value() const {return mvn_mean_func_space_value;}

    const ColumnVector& get_hrf_value() const {return hrf_value;}

    float get_kappa_value() const {return kappa_value;}
    float get_gamma_value() const {return gamma_value;}
    float get_tau_value() const {return tau_value;}
    float get_E0_value() const {return E0_value;}
    float get_alpha_value() const {return alpha_value;}
    float get_logepsilon_value() const {return logepsilon_value;}

    void set_value_hrf(const vector<float>& hrf_vals); // this sets the hrf values for those being estimated
    const vector<float> get_value_hrf() const; // this returns the hrf values for those being estimated
    const vector<string> get_names_hrf() const; // this returns the hrf names for those being estimated
    const vector<float> get_mins_hrf() const; // this returns the hrf mins for those being estimated
    const vector<float> get_maxs_hrf() const; // this returns the hrf maxs for those being estimated

    void set_value_balloon_cbf(const vector<float>& balloon_vals); // this sets the balloon values for those being estimated
    const vector<float> get_value_balloon_cbf() const; // this returns the balloon values for those being estimated
    const vector<string> get_names_balloon_cbf() const; // this returns the balloon names for those being estimated
    const vector<float> get_mins_balloon_cbf() const; // this returns the balloon mins for those being estimated
    const vector<float> get_maxs_balloon_cbf() const; // this returns the balloon maxs for those being estimated
 
    void set_value_balloon(const vector<float>& balloon_vals); // this sets the balloon values for those being estimated
    const vector<float> get_value_balloon() const; // this returns the balloon values for those being estimated
    const vector<string> get_names_balloon() const; // this returns the balloon names for those being estimated
    const vector<float> get_mins_balloon() const; // this returns the balloon mins for those being estimated
    const vector<float> get_maxs_balloon() const; // this returns the balloon maxs for those being estimated
 
    void set_value_balloon2(const vector<float>& balloon_vals); // this sets the balloon values for those being estimated
    const vector<float> get_value_balloon2() const; // this returns the balloon values for those being estimated
    const vector<string> get_names_balloon2() const; // this returns the balloon names for those being estimated
    const vector<float> get_mins_balloon2() const; // this returns the balloon mins for those being estimated
    const vector<float> get_maxs_balloon2() const; // this returns the balloon maxs for those being estimated
       
    void set_value_mvn_mean_standard_space(const vector<float>& mvn_mean_standard_space_vals); // this sets the mvn_mean_standard_space values for those being estimated
    void set_full_value_mvn_mean_func_space(const ColumnVector& pmean_func_space_value);

    const vector<float> get_value_mvn_mean_standard_space() const; // this returns the mvn_mean_standard_space values for those being estimated
   
    const vector<string> get_names_mvn_mean_standard_space() const; // this returns the mvn_mean_standard_space names for those being estimated
    const vector<float> get_mins_mvn_mean_standard_space() const; // this returns the mvn_mean_standard_space mins for those being estimated
    const vector<float> get_maxs_mvn_mean_standard_space() const; // this returns the mvn_mean_standard_space maxs for those being estimated
 
    const vector<float>& get_mvn_mean_standard_space_value() const { return mvn_mean_standard_space_value; } // means are passed-in (MCMC sampled/estimated) in std space


    void set_full_value_mvn_sqrt_cov(const SymmetricMatrix& psqrt_cov_value);

    void set_value_mvn_sqrt_cov(const vector<float>& mvn_sqrt_cov_vals); // this sets the mvn_sqrt_cov values for those being estimated    
    const vector<float> get_value_mvn_sqrt_cov() const; // this returns the mvn_sqrt_cov values for those being estimated
    const vector<string> get_names_mvn_sqrt_cov() const; // this returns the mvn_sqrt_cov names for those being estimated
    const vector<float> get_mins_mvn_sqrt_cov() const; // this returns the mvn_sqrt_cov mins for those being estimated
    const vector<float> get_maxs_mvn_sqrt_cov() const; // this returns the mvn_sqrt_cov maxs for those being estimated
    const SymmetricMatrix& get_mvn_sqrt_cov_value() const {return mvn_sqrt_cov_value;}

    void set_value_mvn_sqrt_cov_offdiag(const vector<float>& mvn_sqrt_cov_offdiag_vals); // this sets the mvn_sqrt_cov_offdiag values for those being estimated
    const vector<float> get_value_mvn_sqrt_cov_offdiag() const; // this returns the mvn_sqrt_cov_offdiag values for those being estimated
    const vector<string> get_names_mvn_sqrt_cov_offdiag() const; // this returns the mvn_sqrt_cov_offdiag names for those being estimated
    const vector<float> get_mins_mvn_sqrt_cov_offdiag() const; // this returns the mvn_sqrt_cov_offdiag mins for those being estimated
    const vector<float> get_maxs_mvn_sqrt_cov_offdiag() const; // this returns the mvn_sqrt_cov_offdiag maxs for those being estimated
 
    const Node& get_node() const { return node; }
    
    // Destructor
    virtual ~Subject_Node() {}
 
    void copy_roi_values(const Subject_Node& sub_nod);

    void copy_hrf_values(const Subject_Node& sub_nod);

    const volume<float>& get_func_mask() const {return func_mask;}

    const ColumnVector& get_prior_mean_balloon_cbf() const {return balloon_cbf_prior_mean;}
    const SymmetricMatrix& get_prior_prec_balloon_cbf() const {return balloon_cbf_prior_cov_inv;}

    const ColumnVector& get_prior_mean_balloon() const {return balloon_prior_mean;}
    const SymmetricMatrix& get_prior_prec_balloon() const {return balloon_prior_cov_inv;}

    const ColumnVector& get_prior_mean_balloon2() const {return balloon2_prior_mean;}
    const SymmetricMatrix& get_prior_prec_balloon2() const {return balloon2_prior_cov_inv;}

  private:

    // constructor
    Subject_Node();
    const Subject_Node& operator=(Subject_Node& par);     
    Subject_Node(Subject_Node& des);
    
  public:  
    string subject_name;

    const Node& node;

    ColumnVector min_func_space_voxel_coord; //nvoxels*3, in functional space
    ColumnVector max_func_space_voxel_coord; //nvoxels*3, in functional space

    Matrix std2func_xform; 
    ColumnVector stdvoxdim;
    ColumnVector funcvoxdim;
    volume<float> func_mask;  // mask for roi in functional space

    volume<float> std_mask; // mask for roi in std space

    // hrf values
    ColumnVector hrf_value; // m1,m2,m3,m4,c1,c2 for halfcos hrf

    // balloon cbf model values
    ColumnVector balloon_cbf_prior_mean;
    SymmetricMatrix balloon_cbf_prior_cov;
    SymmetricMatrix balloon_cbf_prior_cov_inv;

   // balloon model values
    ColumnVector balloon_prior_mean;
    SymmetricMatrix balloon_prior_cov;
    SymmetricMatrix balloon_prior_cov_inv;
  
    // balloon model values
    ColumnVector balloon2_prior_mean;
    SymmetricMatrix balloon2_prior_cov;
    SymmetricMatrix balloon2_prior_cov_inv;

    float kappa_value;
    float gamma_value;
    float tau_value;
    float alpha_value;
    float E0_value;
    float logepsilon_value;

    // mvn values
    float mvn_norm;
    SymmetricMatrix covar_inv;

    // MVN exists in functional space for each subject

    ColumnVector mvn_mean_func_space_value;  // note that means are passed-in (MCMC sampled/estimated) in std space
    vector<float> mvn_mean_standard_space_value;  // means are passed-in (MCMC sampled/estimated) in std space

    SymmetricMatrix mvn_cov_value;    // note that covs are passed-in (MCMC sampled/estimated) as sqrt in functional space 
    SymmetricMatrix mvn_sqrt_cov_value;    // covs are passed-in (MCMC sampled/estimated) as sqrt in functional space 

    bool decode;

    int num_basis; // total num of basis fns
    vector<int> basis_nums; // num basis fns in each dimension
    vector<float> basis_res; // resolution of basis in each dimension  (in func space)
    vector<ColumnVector> basis_centres; // coords of basis fns: num_dims*basis_nums(in each dim)  (in func space)
    vector<ColumnVector> basis_centres_list; // list containing each basis fn centre: num_basis*num_dims  (in func space)
  

    mutable Mvn mvn;
  };


}   
#endif







