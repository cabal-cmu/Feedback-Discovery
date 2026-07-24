/*  nma_generate_options.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(NmaGenerateOptions_h)
#define NmaGenerateOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace Nma {

  class NmaGenerateOptions {
  public:
    static NmaGenerateOptions& getInstance();
    ~NmaGenerateOptions() { delete gopt; }
  
    Option<int> debuglevel;
    Option<bool> timingon;
    Option<bool> help;
    Option<float> tr;
    Option<float> sigmaa;
    Option<int> resfactor;
    Option<string> data_directory;
    Option<string> model_directory;
    Option<string> logdir;
    Option<vector<string> > node_names;
    Option<string> subject_name;
    Option<vector<string> > stimuli_names;
    Option<int> data_mode;
    Option<int> num_scans;
    Option<bool> decode;
    Option<bool> phi_every_voxel;
    Option<string> haemodynamic_model;
    Option<bool> show_fit;
    Option<bool> stimuli_are_single_col_format;

    void parse_command_line(int argc, char** argv, Log& logger);
  
  private:

    NmaGenerateOptions();
    const NmaGenerateOptions& operator=(NmaGenerateOptions&);
    NmaGenerateOptions(NmaGenerateOptions&);

    OptionParser options; 
      
    static NmaGenerateOptions* gopt;
  
  };

  inline NmaGenerateOptions& NmaGenerateOptions::getInstance(){
    if(gopt == NULL)
      gopt = new NmaGenerateOptions();
   
    return *gopt;
  }
  
  inline NmaGenerateOptions::NmaGenerateOptions() :
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
    sigmaa(string("--sigmaa"), 1,
       string("sigma_a"),
       false, requires_argument),
    resfactor(string("--rf,--resfactor"), 6,
	      string("resolution expressed as a factor, i.e. res (in secs) = tr/resfactor"),
	      true, requires_argument),
    data_directory(string("--dd,--datadirectory"), string("data_directory"),
		   string("data directory, contains node standard space masks and stimulus files"),
		   true, requires_argument), 
    model_directory(string("--md,--modeldirectory"), string("model_directory"),
		    string("model directory, contains A,B,C,D matrices with names matA, matB, matC, matD"),
		    true, requires_argument), 
    logdir(string("--ld,--logdir"), string("logdir"),
	   string("log directory"),
	   false, requires_argument),  
    node_names(string("--nodenames,--nn"), vector<string>(), 
	       string("comma separated list of node names"), 
	       true, requires_argument),
    subject_name(string("--subjectname,--sub,--sn"), string(""), 
		  string("subject name"), 
		  true, requires_argument),
    stimuli_names(string("--stimulinames,--stim"), vector<string>(), 
		  string("comma separated list of stimuli names"), 
		  true, requires_argument),
    data_mode(string("--dm,--datamode"), 0, 
	      string("data mode, 0=uses single time series input (default), 1=uses full node data and infers ROI, 2=uses single time series data averaged across node mask input"), 
	      false, requires_argument),
    num_scans(string("--numscans,--ns"), 200,
	      string("Number of scans"),
	      true, requires_argument),
    decode(string("--decode"), false, 
		string("turns on spatial decoding"), 
		false, no_argument),
    phi_every_voxel(string("--pev"), false, 
		string("model voxelwise variance (default is model a global variance)"), 
		false, no_argument),
    haemodynamic_model(string("--hm,--haemodynamicmodel"), string("balloon"), 
		       string("String specifying haemodynamic model: [halfcosine] Linear half cosine HRF, [balloon] Use balloon model with no epsilon (default), [balloon_epsilon] Use balloon model with epsilon "), 
		       false, requires_argument),
    show_fit(string("--sf,--showfit"), false, 
	     string("Show fit to subject data"), 
	     false, no_argument),
    stimuli_are_single_col_format(string("--sascf"), false, 
		  string("stimuli files are single column format (default is 3 col format)"), 
		  false, no_argument),
    options("nma_generate","")    
  {
    try {
      options.add(debuglevel);
      options.add(timingon);
      options.add(help);
      options.add(tr);
      options.add(sigmaa);
      options.add(resfactor);
      options.add(data_directory);
      options.add(model_directory);
      options.add(logdir);
      options.add(node_names);
      options.add(subject_name);
      options.add(stimuli_names);
      options.add(data_mode);
      options.add(num_scans);
      options.add(decode);
      options.add(phi_every_voxel);
      options.add(haemodynamic_model);
      options.add(show_fit);
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



