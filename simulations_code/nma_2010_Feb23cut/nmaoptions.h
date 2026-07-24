/*  nmaoptions.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(NmaOptions_h)
#define NmaOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace Nma {

  class NmaOptions {
  public:
    static NmaOptions& getInstance();
    ~NmaOptions() { delete gopt; }
  
    Option<int> debuglevel;
    Option<bool> timingon;
    Option<bool> help;
    Option<float> tr;
    Option<int> resfactor;
    Option<string> data_directory;
    Option<string> model_directory;
    Option<string> logdir;
    Option<int> nsamps;
    Option<int> burnin;
    Option<int> sampleevery;
    Option<vector<string> > node_names;
    Option<vector<string> > subject_names;
    Option<vector<string> > stimuli_names;
    Option<vector<int> > stim_amp_mod;
    Option<vector<int> > stim_ard;
    Option<int> data_mode;
    Option<int> seed;
    Option<bool> output_samples;
    Option<bool> no_init_roi;
    Option<bool> decode;
    Option<bool> phi_every_voxel;
    Option<string> haemodynamic_model;
    Option<bool> rand_init;
    Option<bool> spatial_crosscorr;
    Option<bool> stimuli_are_single_col_format;

    void parse_command_line(int argc, char** argv, Log& logger);
  
  private:

    NmaOptions();
    const NmaOptions& operator=(NmaOptions&);
    NmaOptions(NmaOptions&);

    OptionParser options; 
      
    static NmaOptions* gopt;
  
  };

  inline NmaOptions& NmaOptions::getInstance(){
    if(gopt == NULL)
      gopt = new NmaOptions();
   
    return *gopt;
  }
  
  inline NmaOptions::NmaOptions() :
    debuglevel(string("--debug,--debuglevel"), 0,
	       string("set debug level"), 
	       false, requires_argument),
    timingon(string("--to,--timingon"), false, 
	     string("turn timing on"), 
	     false, no_argument),
    help(string("-h,--help"), false,
	 string("display this message"),
	 false, no_argument),
    tr(string("--tr"), 3,
       string("TR(secs)"),
       true, requires_argument),
    resfactor(string("--rf,--resfactor"), 10,
	string("resolution expressed as a factor, i.e. res (in secs) = tr/resfactor, recommendation is that this is chosen so that res is approx 0.2 secs."),
	true, requires_argument),
    data_directory(string("--dd,--datadirectory"), string("data_directory"),
		    string("data directory, contains standard space image, node standard space masks and in the subject's subdirectory contains standard2func space xform, 3 column format stimulus files and the functional data"),
		    true, requires_argument), 
    model_directory(string("--md,--modeldirectory"), string("model_directory"),
		    string("model directory, which contains A,B,C,D matrices as text files"),
		    true, requires_argument), 
    logdir(string("--ld,--logdir"), string("logdir"),
	   string("log directory"),
	   false, requires_argument),  
    nsamps(string("--ns,--nsamps"), 3000,
	   string("Num of samps to be made by MCMC"),
	   false, requires_argument),
    burnin(string("--bi,--burnin"), 3000,
	   string("Num of samps at start of MCMC to be discarded"),
	   false, requires_argument),
    sampleevery(string("--se,--sampleevery"), 2,
		string("Num of jumps for each sample"),
		false, requires_argument),
    node_names(string("--nodenames,--nn"), vector<string>(), 
	       string("comma separated list of node names"), 
	       true, requires_argument),
    subject_names(string("--subjectname,--sub"), vector<string>(), 
		  //string("comma separated list of subject names"), 
		  string("subject name"), 
		  true, requires_argument),
    stimuli_names(string("--stimulinames,--stim"), vector<string>(), 
		  string("comma separated list of stimuli names"), 
		  true, requires_argument),
    stim_amp_mod(string("--sam"), vector<int>(), 
		 string("Comma separated boolean list indicating which stimuli are amplitude modulated"), 
		 true, requires_argument),
    stim_ard(string("--stimard"), vector<int>(), 
		 string("Comma separated boolean list indicating which stimuli have ARD priors for c matrix inputs"), 
		 true, requires_argument),
    data_mode(string("--dm,--datamode"), 0, 
	      string("data mode, 0=uses single time series input (default), 1=uses full node data and infers ROI, 2=uses single time series data averaged across node mask input"), 
	      false, requires_argument),
    seed(string("--seed"), 10, 
	 string("seed for pseudo random number generator"), 
	 false, requires_argument),
    output_samples(string("--os,--outputsamples"), false, 
		   string("output MCMC samples"), 
		   false, no_argument),
    no_init_roi(string("--nir,--noinitroi"), false, 
		string("turns off initialisation of ROI multivariate normal"), 
		false, no_argument),
    decode(string("--decode"), false, 
		string("turns on spatial decoding"), 
		false, no_argument),
    phi_every_voxel(string("--pev"), false, 
		string("model voxelwise variance (default is model a global variance)"), 
		false, no_argument),
    haemodynamic_model(string("--hm,--haemodynamicmodel"), string("balloon"), 
		string("String specifying haemodynamic model: [halfcosine] Linear half cosine HRF, [balloon] Use balloon model with no epsilon (default), [balloon_epsilon] Use balloon model with epsilon "), 
		false, requires_argument),
    rand_init(string("--nri,--norandinit"), true, 
		string("no random initialisation - uses values passed in in the A, B, C and D matrices"), 
		false, no_argument),
    spatial_crosscorr(string("--scc,--spatial_crosscorr"), false, 
		string("model non-zero cross-correlations in the spatial 3D Gaussian models of the ROIs"), 
		false, no_argument),
    stimuli_are_single_col_format(string("--sascf"), false, 
		  string("stimuli files are single column format (default is 3 col format)"), 
		  false, no_argument),
    options("nma","nma --ld=logdir --dd=nma_data --md=nma_data/model1 --tr=3 --resfactor=16 --sub=subject1 --stim=attention,photic,motion --nn=v1,v5,ppc --sam=0,0,0 --bi=3000 --ns=3000 --se=2 --dm=2 --ubm")
  {
    try {
      options.add(debuglevel);
      options.add(timingon);
      options.add(help);
      options.add(tr);
      options.add(resfactor);
      options.add(data_directory);
      options.add(model_directory);
      options.add(logdir);
      options.add(nsamps);
      options.add(burnin);
      options.add(sampleevery);
      options.add(node_names);
      options.add(subject_names);
      options.add(stimuli_names);
      options.add(stim_amp_mod);
      options.add(stim_ard);
      options.add(data_mode);
      options.add(seed);
      options.add(output_samples);
      options.add(no_init_roi);
      options.add(decode);
      options.add(phi_every_voxel);
      options.add(haemodynamic_model);
      options.add(rand_init);
      options.add(spatial_crosscorr);
      options.add(stimuli_are_single_col_format);
    }
    catch(X_OptionError& e) {
      options.usage();
      cerr << endl << e.what() << endl;
    } 
    catch(std::exception &e) {
      cerr << e.what() << endl;
    }    
     
  }
}

#endif



