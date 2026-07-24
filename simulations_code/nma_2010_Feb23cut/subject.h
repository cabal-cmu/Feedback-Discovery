/*  subject.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(subject_h)
#define subject_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "nmaoptions.h"
#include "newimage/newimageall.h"
#include "model.h"

using namespace NEWIMAGE;
using namespace MISCMATHS;

namespace Nma {
    
  class Subject
    {
    public:

      // constructor
      Subject(const string& pname, const string&  pdata_dir, const Model& pmodel, bool psingle_timeseries, int pdata_mode, string phaemodynamic_model, bool prandom_initialise, int pdebuglevel, bool pdecode, bool pphi_every_voxel);

      void setup();

      // Destructor
      virtual ~Subject() {}
 
      Subject_Model& get_subject_model() { return subject_model; }

      const Subject_Model& get_subject_model() const { return subject_model; }

      bool is_single_timeseries() const { return single_timeseries; }
      bool is_decode() const { return decode; }
      bool is_phi_every_voxel() const { return phi_every_voxel; }

      const volume4D<float>& get_data() const {return data;}
      //  volumeinfo& get_data_volinfo()  {return data_volinfo;}

      const ColumnVector& get_vb_data();

      const volume<float>& get_func_mask() const { return  func_mask; }
      const volume<float>& get_standard_brain() const {return standard_brain;}
 
      const vector<ColumnVector>& get_node_data() const {return node_data;}

      const Matrix& get_voxelwise_data() const {return voxelwise_data;}
      string get_name() const {return name;}
      const ColumnVector& get_stdvoxdim() const {return stdvoxdim;}
      const ColumnVector& get_funcvoxdim() const {return funcvoxdim;}

    private:

      Subject();
      const Subject& operator=(Subject& par);     
      Subject(Subject& des);
      
      // data members
      bool single_timeseries;
      string haemodynamic_model;

      bool decode;
      bool phi_every_voxel;

      int data_mode;

      string name; 
      string data_dir;

      Subject_Model subject_model;

      volume4D<float> data;
      //      volumeinfo data_volinfo;
      volume<float> func_mask;
      volume<float> standard_brain;
      
      Matrix voxelwise_data; //  num_voxels_in_allnodes*nscans)

      vector<ColumnVector> node_data;

      ColumnVector vb_data;
      bool already_setup_vb_data;

      int nnodes;
      int nscans;

      Matrix std2func_xform; 
      ColumnVector stdvoxdim;
      ColumnVector funcvoxdim;

      bool random_initialise;

      int debuglevel;
    };

  void create_pvf_overlays(const string& name, const volume<float>& pvf, const volume<float>& epivol, const volume<float>& epivolmask, const ColumnVector& coords);//, volumeinfo& data_volinfo);

  void create_report_fit(Subject& subject, string base_name);

  void create_report_signal(Subject_Model& subject_model, string base_name);

  double calculate_model_energy(const Subject& subject, int debuglevel);
  double calculate_model_prior_energy(const Subject& subject, int debuglevel);

  void evaluate_model_voxelwise_energy(const Subject& subject, int vox, int debuglevel);

  double addup_model_voxelwise_energy(const Subject& subject, int debuglevel);

  void calculate_forward_model(const Subject&, ColumnVector&);
  void calculate_noise_precision(const Subject&, ColumnVector&);
}   
#endif







