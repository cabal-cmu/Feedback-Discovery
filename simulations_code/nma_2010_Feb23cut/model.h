/*  model.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(model_h)
#define model_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "newimage/newimageall.h"
#include "node.h"
#include "halfcos_hrf.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

namespace Nma {
  
  class Subject_Model;
  class Subject;

  class Model
    {
    public:

      Model(bool psingle_timeseries, int pdata_mode, string phaemodynamic_model, const string& pdata_dir, const string& pmodel_dir, const vector<string>& pnode_names, const vector<string>& pstimuli_names, const vector<int>& pstim_amp_mod, const vector<int>& pstim_ard, bool pstim_single_col_format, double pres, float ptr, bool pdecode);

      void setup();

      // saves results in logging directory
      void save();

      const vector<string>& get_node_names() const {return node_names;}
      const vector<Node*>& get_nodes() const {return nodes;}   
      const vector<int>& get_stim_amp_mod() const {return stim_amp_mod;}
      const vector<int>& get_stim_ard() const {return stim_ard;}
      const vector<string>& get_stimuli_names() const {return stimuli_names;}

      bool is_b_amp_mod(int n, int i) const // n indexes node, i indexes inplay (nnodes*nstim)
      {
	int stim_index=marker_b[n-1][i-1].second;
	return stim_amp_mod[stim_index-1];
      }
      bool is_c_amp_mod(int n, int i) const // n indexes node, i indexes inplay (nnodes*nstim)
      {
	int stim_index=marker_c[n-1][i-1];
	return stim_amp_mod[stim_index-1];
      }

      bool is_c_ard(int n, int i) const // n indexes node, i indexes inplay (nnodes*nstim)
      {
	int stim_index=marker_c[n-1][i-1];
	return stim_ard[stim_index-1];
      }

      int get_num_a_params() const { return num_a_params; }
      int get_num_b_params() const { return num_b_params; }
      int get_num_c_params() const { return num_c_params; }
      int get_num_d_params() const { return num_d_params; }

      const vector<vector<int> >& get_marker_a() const {return marker_a;}
      const vector<vector<pair<int,int> > >& get_marker_b() const {return marker_b;}
      const vector<vector<int> >& get_marker_c() const {return marker_c;}
      const vector<vector<pair<int,int> > >& get_marker_d() const {return marker_d;}


      const bool is_single_timeseries() const {return single_timeseries;}
      const bool is_decode() const {return decode;}

      // Destructor
      virtual ~Model() { 
	//	LOGOUT("~Model start");
	for(int r=0; r<nnodes; r++) delete nodes[r]; nodes.clear();
	//	LOGOUT("~Model end");
      }
 
    private:
      Model(){}
     const Model& operator=(const Model& mod);
      Model(const Model& mod)
      {
	*this = mod;
      }

      bool single_timeseries;
      string haemodynamic_model;

      bool decode;
      int data_mode;

      string data_dir;
      string model_dir;

      vector<string> node_names;
      vector<string> stimuli_names;
      vector<int> stim_amp_mod;
      vector<int> stim_ard;

      bool stim_single_col_format;
      string confound_evs_name;

      double res; // in secs
      float tr; // in secs
 
      vector<Node*> nodes;
      int nnodes;

      // model specification
      Matrix matA;
      vector<Matrix> matB;
      Matrix matC;
      vector<Matrix> matD;

      // markers
      vector<vector<int> > marker_a; // nnodes*nnodes
      vector<vector<pair<int,int> > > marker_b; // nnodes*(nnodes*nstim) pair(node_in_index, stim_index)
      vector<vector<int> > marker_c; // nnodes*nstim
      vector<vector<pair<int,int> > > marker_d; // nnodes*(nnodes*nodes) pair(node_in1_index, node_in2_index)
         
      int num_a_params;
      int num_b_params;
      int num_c_params;
      int num_d_params;

      friend class Subject_Model;
      friend class Subject;

      friend void fit_model_using_c_matrix_only(Subject& subject, int debuglevel);
      friend void find_best_vox(Subject& subject, int debuglevel);
    };

  class Subject_Model
  {
  public:

    Subject_Model(string psubject_name, const Model& pmodel, string haemodynamic_model, int pdebuglevel);

    void copy_neural_connectivity_values(const Subject_Model& sub_model);
    void copy_hrf_values(const Subject_Model& sub_model);
    void copy_roi_values(const Subject_Model& sub_model);
    void output_roi_values();

    void setup();
    void initialise(bool random_initialise);

    const Model& get_model() const {return model;}
    const vector<Subject_Node*>& get_subject_nodes() const {return subject_nodes;}   
    vector<Subject_Node*>& get_subject_nodes()  {return subject_nodes;}   

    // Destructor
    virtual ~Subject_Model() { //LOGOUT("~Subject_Model start");
      for(int r=0; r<nnodes; r++) delete subject_nodes[r]; subject_nodes.clear(); 
      //LOGOUT("~Subject_Model end");
    }
 
    void evaluate_neuronal_activity();

    void evaluate_hrf();
    void evaluate_hrf(int n);

    void get_hrf(vector<ColumnVector>& hrf, ColumnVector& t_hrf, float max_secs);
    void get_nrf(vector<vector<ColumnVector> >& nrf, ColumnVector& t_nrf, float max_secs);

    void evaluate_node_bold();
    void evaluate_node_bold(int n);
    void evaluate_node_bold_halfcosine(int n);
    void evaluate_node_cbf_balloon(int n);
    void evaluate_node_bold_balloon(int n);

    void evaluate_voxelwise_bold();
    void evaluate_voxelwise_bold(int n);

    void evaluate_decoded_node_data();
    void evaluate_decoded_node_data(int n);

    const vector<float> get_value_a_vec() const;    
    void set_value_a(const vector<float>& pvalues);

    const vector<float> get_value_a_vec(int n) const;    
    void set_value_a(int n, const vector<float>& pvalues);

    const vector<float> get_value_b_vec() const;    
    void set_value_b(const vector<float>& pvalues);

    const vector<float> get_value_b_amp_mod_vec(int n, int i) const;// n indexes node, i indexes inplay (nnodes*nstim)
    void set_value_b_amp_mod(int n, int i, const vector<float>& pvalues);

    const vector<float> get_value_b_vec(int n) const;    
    void set_value_b(int n, const vector<float>& pvalues);

    const vector<bool> get_isard_c_vec() const;    
    const vector<float> get_value_c_vec() const;    
    void set_value_c(const vector<float>& pvalues);

    const vector<float> get_value_c_amp_mod_vec(int n, int i) const;// n indexes node, i indexes inplay nstim)
    void set_value_c_amp_mod(int n, int i, const vector<float>& pvalues); 

    const vector<bool> get_isard_c_vec(int n) const;    
    const vector<float> get_value_c_vec(int n) const;    
    void set_value_c(int n, const vector<float>& pvalues);

    const vector<float> get_value_d_vec() const;    
    void set_value_d(const vector<float>& pvalues);

    const vector<float> get_value_d_vec(int n) const;    
    void set_value_d(int n, const vector<float>& pvalues);

    const vector<float> get_value_logsigmaa_vec() const;    
    void set_value_logsigmaa(const vector<float>& pvalues);
    float get_prior_mean_logsigmaa() const {return prior_mean_logsigmaa;}
    float get_prior_precision_logsigmaa() const {return prior_precision_logsigmaa;}

    const vector<string> get_names_a() const; 
    const vector<string> get_names_a(int n) const;
    const vector<string> get_names_b() const; 
    const vector<string> get_names_b_amp_mod(int n, int i) const; // n indexes node, i indexes inplay (nnodes*nstim)
    const vector<string> get_names_b(int n) const;
    const vector<string> get_names_c() const; 
    const vector<string> get_names_c_amp_mod(int n, int i) const; // n indexes node, i indexes inplay nstim)
    const vector<string> get_names_c(int n) const;
    const vector<string> get_names_d() const; 
    const vector<string> get_names_d(int n) const;
    const vector<string> get_names_logsigmaa() const; 

    const vector<ColumnVector>& get_node_cbf() const {return node_cbf;}    
    const vector<ColumnVector>& get_node_bold() const {return node_bold;}
    const Matrix& get_voxelwise_bold() const {return voxelwise_bold;}
    const vector<ColumnVector>& get_voxelwise_pvf() const {return voxelwise_pvf;}
    void set_node_bold(const vector<ColumnVector>& in)  { node_bold=in;}
    void set_node_cbf(const vector<ColumnVector>& in)  { node_cbf=in;}
    void set_voxelwise_bold(const Matrix& in)  { voxelwise_bold=in;}
    void set_voxelwise_pvf(const vector<ColumnVector>& in)  { voxelwise_pvf=in;}
 
    const Matrix& get_residual_forming_confound_evs() const { return residual_forming_confound_evs; }

    int get_num_b_amp_mod_params(int n, int i) const { return num_b_amp_mod_params[n-1][i-1]; } // n indexes node, i indexes inplay (nnodes*nstim)
    int get_num_c_amp_mod_params(int n, int i) const { return num_c_amp_mod_params[n-1][i-1]; } // n indexes node, i indexes inplay nstim)

    int get_num_a_node_params(int n) const { return num_a_node_params[n-1]; }
    int get_num_b_node_params(int n) const { return num_b_node_params[n-1]; }
    int get_num_c_node_params(int n) const { return num_c_node_params[n-1]; }
    int get_num_d_node_params(int n) const { return num_d_node_params[n-1]; }

    const vector<vector<ColumnVector> >& get_stimuli() const {return stimuli;}

    const vector<vector<float> >& get_value_a() const {return value_a;}
    const vector<vector<float> >& get_value_b() const {return value_b;}
    const vector<vector<float> >& get_value_c() const {return value_c;}
    const vector<vector<float> >& get_value_d() const {return value_d;}
    const vector<vector<vector<float> > >& get_value_b_amp_mod() const {return value_b_amp_mod;}
    const vector<vector<vector<float> > >& get_value_c_amp_mod() const {return value_c_amp_mod;}
    const float get_value_logsigmaa() const {return value_logsigmaa;}

    const vector<ColumnVector>& get_z() const {return z;}
    const vector<ColumnVector>& get_zfft_real() const {return zfft_real;} 
    const vector<ColumnVector>& get_zfft_imag() const {return zfft_imag;} 
    const ColumnVector& get_z_mean() const {return z_mean;} 
    const vector<ColumnVector>& get_hrf_fft_real() const {return hrf_fft_real;}
    const vector<ColumnVector>& get_hrf_fft_imag() const {return hrf_fft_imag;}
  
    void set_z(const vector<ColumnVector>& in)  { z=in;}
    void set_zfft_real(const vector<ColumnVector>& in)  { zfft_real=in;} 
    void set_zfft_imag(const vector<ColumnVector>& in)  { zfft_imag=in;} 
    void set_z_mean(const ColumnVector& in)  { z_mean=in;} 
    void set_hrf_fft_real(const vector<ColumnVector>& in)  { hrf_fft_real=in;}
    void set_hrf_fft_imag(const vector<ColumnVector>& in)  { hrf_fft_imag=in;}

    void setup_voxelwise_data(const volume4D<float>& data, Matrix& voxelwise_data);
    void set_voxelwise_data(const Matrix& pvoxelwise_data){voxelwise_data=pvoxelwise_data;}
    void set_node_data(const vector<ColumnVector>& pnode_data){node_data=pnode_data;}
    void establish_voxel_coordinates(const volume<float>& mask);
    const vector<ColumnVector>& get_voxel_coordinates() const {return voxel_coordinates;}

    float get_logsigmaa_prior_energy(float psig)
    {
      return Sqr(value_logsigmaa-prior_mean_logsigmaa)/(2.0*prior_precision_logsigmaa);
    }

    float get_partial_vol_fraction(int vox, int nod){
      
      int maxr=0;
      float maxpvf=0;
      for(unsigned int r=1; r<=subject_nodes.size(); r++)
	{
	  if(is_voxel_in_roi(vox,r))
	    {
	      
	      float pvftmp=subject_nodes[r-1]->get_partial_vol_fraction(voxel_coordinates[vox-1]);
	      
// 	      OUT(pvftmp);
	      
	      if(pvftmp>maxpvf)
		{
		  maxpvf=pvftmp;
		  maxr=r;
		}
	    }
	}


      float ret;
      if(maxr==nod)
	ret=maxpvf;
      else
	ret=0.0001;

//       OUT(maxpvf);
//       OUT(maxr);
//       OUT(nod);
//       OUT(ret);

      return ret;

    }

    void print_partial_vol_fraction(int vox, int nod){
      
      int maxr=0;
      float maxpvf=0;
      for(unsigned int r=1; r<=subject_nodes.size(); r++)
	{
	  if(is_voxel_in_roi(vox,r))
	    {	      
	      float pvftmp=subject_nodes[r-1]->print_partial_vol_fraction(voxel_coordinates[vox-1]);
	      	      
	      if(pvftmp>maxpvf)
		{
		  maxpvf=pvftmp;
		  maxr=r;
		}
	    }
	}


      float ret;
      if(maxr==nod)
	ret=maxpvf;
      else
	ret=0.0001;

      OUT(maxpvf);
      OUT(maxr);
      OUT(nod);
      OUT(ret);


    }
    
    bool is_voxel_in_roi(int vox, int nod) const {return bool(voxel_in_roi[vox-1](nod));}    

    void set_nsecs(double pnsecs) {
      nsecs=pnsecs;    
      ntpts=int(::round(nsecs/model.res));
    }

    string get_haemodynamic_model() const {return haemodynamic_model;}

    void set_debuglevel(int pdebuglevel) {debuglevel=pdebuglevel;}

    const ColumnVector& get_log_phi_voxel() const {return log_phi_voxel;}
    const ColumnVector& get_log_phi_every_voxel() const {return log_phi_every_voxel;}
    const vector<float>& get_phi_node() const {return phi_node;}

    const vector<float> get_value_log_phi_every_voxel_vec(int vox) const {vector<float> ret; ret.push_back(log_phi_every_voxel(vox)); return ret;}    
    void set_value_log_phi_every_voxel(const vector<float>& pvalues, int vox) {log_phi_every_voxel(vox)=pvalues[0]; }
    const vector<string> get_names_log_phi_every_voxel(int vox) const {vector<string> ret; ret.push_back(string("log_phi_every_voxel_"+num2str(vox))); return ret;}

    const vector<float> get_value_log_phi_voxel_vec() const {vector<float> ret; ret.push_back(log_phi_voxel(1)); return ret;}    
    void set_value_log_phi_voxel(const vector<float>& pvalues) {log_phi_voxel(1)=pvalues[0]; }
    const vector<string> get_names_log_phi_voxel() const {vector<string> ret; ret.push_back(string("log_phi_voxel")); return ret;}

    const vector<float> get_value_phi_node_vec(int node_index) const {vector<float> ret; ret.push_back(std::log(phi_node[node_index-1])); return ret;}    
    void set_value_phi_node(const vector<float>& pvalues, int node_index) {phi_node[node_index-1]=std::exp(pvalues[0]);}
    const vector<string> get_names_phi_node(int node_index) const {vector<string> ret; ret.push_back(string("log_phi_node_"+num2str(node_index))); return ret;}

    const vector<ColumnVector>& get_decoded_node_data() const {return decoded_node_data;}

    const ColumnVector& get_voxelwise_energy() const { return voxelwise_energy; }
    void set_voxelwise_energy(const ColumnVector& in) const {  voxelwise_energy=in; }
 
    double get_voxelwise_energy(int vox) const { return voxelwise_energy(vox); }
    void set_voxelwise_energy(int vox, double in) const {  voxelwise_energy(vox)=in; }

    
    
  private:
    
    const Model& model;

    string subject_name;

    double nsecs; // length of experiment, in secs 
    int ntpts; // high res
    int nnodes;

    vector<vector<float> > value_a; // nnodes*(inplay nnodes)
    vector<vector<float> > value_b; // nnodes*(inplay nnodes*nstim)
    vector<vector<float> > value_c; // nnodes*(inplay nstim)
    vector<vector<float> > value_d; // nnodes*(inplay nnodes*nnodes)   
    
    vector<vector<vector<float> > > value_b_amp_mod; // nnodes*(inplay nnodes*nstim)*nevents    
    vector<vector<vector<float> > > value_c_amp_mod; // nnodes*(inplay nstim)*nevents

    vector<vector<int> > num_b_amp_mod_params; // nnodes*(inplay nnodes*nstim)
    vector<vector<int> > num_c_amp_mod_params; // nnodes*(inplay nstim)
 
    vector<int> num_a_node_params;// nnodes
    vector<int> num_b_node_params;// nnodes
    vector<int> num_c_node_params;// nnodes
    vector<int> num_d_node_params;// nnodes

    float value_logsigmaa;
    float prior_mean_logsigmaa;
    float prior_precision_logsigmaa;

    Matrix confound_evs;
    Matrix residual_forming_confound_evs;
    
    vector<vector<ColumnVector> > stimuli; // nstim*nevents*ntpts

    vector<Subject_Node*> subject_nodes;

    Halfcos_Hrf halfcos_hrf;

    vector<ColumnVector> node_cbf;
    vector<ColumnVector> node_bold;
    vector<ColumnVector> node_bold_high_res;
    Matrix voxelwise_bold; // is num_voxels_in_allnodes*nscans
    vector<ColumnVector> voxelwise_pvf; // nnodes*num_voxels_in_allnodes

    vector<ColumnVector> z;
    vector<ColumnVector> zfft_real; 
    vector<ColumnVector> zfft_imag; 
    ColumnVector z_mean; 
    vector<ColumnVector> hrf_fft_real;
    vector<ColumnVector> hrf_fft_imag;

    string haemodynamic_model;

    vector<ColumnVector> voxel_coordinates; //nvoxels*3, in functional space
    vector<ColumnVector> voxel_in_roi; //nvoxels*nnodes, in functional space

    ColumnVector log_phi_voxel;
    ColumnVector log_phi_every_voxel;
    vector<float> phi_node;

    friend class Subject;

    int debuglevel;
    int debugcount1;
    int debugcount2;

    Matrix voxelwise_data;
    vector<ColumnVector> node_data;
    vector<ColumnVector> decoded_node_data;

    mutable ColumnVector voxelwise_energy;

    friend void fit_model_using_c_matrix_only(Subject& subject, int debuglevel);
  };

}   
#endif







