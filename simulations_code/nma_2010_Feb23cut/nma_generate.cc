/*  NMA

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH
#include "newmatap.h"
#include "newmatio.h"
#include <string>
#include <math.h>
#include "stdlib.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"
#include "libvis/miscplot.h"
#include "libvis/miscpic.h"
#include "nma_generate_options.h"
#include "model.h"
#include "subject.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace Nma;
using namespace MISCMATHS;
using namespace MISCPLOT;
using namespace MISCPIC;

int main(int argc, char *argv[])
{
  try
  {  
    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    
    // parse command line - will output arguments to logfile
    NmaGenerateOptions& opts = NmaGenerateOptions::getInstance();
    opts.parse_command_line(argc, argv, logger);

    if(opts.debuglevel.value()==1)
      Tracer_Plus::setrunningstackon();

    if(opts.timingon.value())
      Tracer_Plus::settimingon();

    /////////////////

    bool sts=true; //single_timeseries

    vector<int> amp_mods(opts.stimuli_names.value().size(),0);
    vector<int> ard(opts.stimuli_names.value().size(),0);

    double res=opts.tr.value()/opts.resfactor.value();
    Model model(sts, opts.data_mode.value(), opts.haemodynamic_model.value(), opts.data_directory.value(), opts.model_directory.value(), opts.node_names.value(), opts.stimuli_names.value(), amp_mods, ard, opts.stimuli_are_single_col_format.value(), res, opts.tr.value(), opts.decode.value());     
  
    string report_name;

    if(opts.show_fit.value())
      {
	report_name="fit";
	Subject subject(opts.subject_name.value(), opts.data_directory.value(), model, sts, opts.data_mode.value(), opts.haemodynamic_model.value(),false,opts.debuglevel.value(), opts.decode.value(), opts.phi_every_voxel.value());	

	Subject_Model& subject_model=subject.get_subject_model();

// 	vector<float> logsigmaa;
// 	logsigmaa.push_back(std::log(opts.sigmaa.value()));
// 	subject_model.set_value_logsigmaa(logsigmaa);

	subject_model.evaluate_neuronal_activity();
	subject_model.evaluate_hrf();
	subject_model.evaluate_node_bold();
	
	if(!model.is_single_timeseries())
	  {	
	    if(model.is_decode())
	      subject_model.evaluate_decoded_node_data();
	    else
	      subject_model.evaluate_voxelwise_bold();	 
	  }
	
	create_report_fit(subject, report_name);
      }
    else
      {
	report_name="signal";
	Subject_Model subject_model(opts.subject_name.value(), model, opts.haemodynamic_model.value(),opts.debuglevel.value());

	subject_model.set_nsecs(opts.num_scans.value()*opts.tr.value());
	subject_model.setup();

	subject_model.initialise(false);
	
// 	vector<float> logsigmaa;
// 	logsigmaa.push_back(std::log(opts.sigmaa.value()));
// 	subject_model.set_value_logsigmaa(logsigmaa);

//	subject_model.evaluate_neuronal_activity();
//	subject_model.evaluate_hrf();
//	subject_model.evaluate_node_bold();
	
	if(!model.is_single_timeseries())
	  {	
	    if(model.is_decode())
	      subject_model.evaluate_decoded_node_data();
	    else
	      subject_model.evaluate_voxelwise_bold();	 
	  }
	///	create_report_signal(subject_model, report_name);   /// removed this on suggestion of Wooly

        {
	  Tracer_Plus trace("newcodefromwooly");
	  int nnodes= subject_model.get_model().get_node_names().size();
  
	  const vector<ColumnVector>& node_bold=subject_model.get_node_bold();
  
	  for(int n=1; n<=nnodes; n++)      
	    {    
	      string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
	      string name=report_name+"_"+node_name+"_bold";
      
	      Matrix tmp=node_bold[n-1].t();
	      write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(name+".txt"));
      
	    }  

	  const vector<ColumnVector>& node_z=subject_model.get_z();

	  for(int n=1; n<=nnodes; n++)      
	    {    
	      string node_name=subject_model.get_subject_nodes()[n-1]->get_node().get_name();
	      string name=report_name+"_"+node_name+"_z";
	          
	      Matrix tmp=node_z[n-1].t();
	      write_ascii_matrix(tmp, LogSingleton::getInstance().appendDir(name+".txt"));

	    }  

	}

      }
        
    /////////////////

    if(opts.timingon.value())
      Tracer_Plus::dump_times(logger.getDir());

    LogSingleton::getInstance().str() << endl << "Log directory was: " << logger.getDir() << endl;

    LogSingleton::getInstance().str() << endl << "For a web page report view: " << endl << LogSingleton::getInstance().appendDir("report_"+report_name+".html")  << endl;
  }
  catch(Exception& e) 
    {
      cerr << endl << e.what() << endl;
      return 1;
    }
  catch(X_OptionError& e) 
    {
      cerr << endl << e.what() << endl;
      return 1;
    }

  return 0;
}












