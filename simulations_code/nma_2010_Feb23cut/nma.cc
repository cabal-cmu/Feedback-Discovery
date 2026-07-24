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
#include "utils/log.h"
#include "nma_manager.h"
//#include "mcmc_hmc.h"
#include "nmaoptions.h"
#include "utils/tracer_plus.h"
#include "miscmaths/miscprob.h"
#include "stdlib.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace Nma;
using namespace MISCMATHS;

int main(int argc, char *argv[])
{
  try
  { 

    // Setup logging:
    Log& logger = LogSingleton::getInstance();
    
    // parse command line - will output arguments to logfile
    NmaOptions& opts = NmaOptions::getInstance();
    opts.parse_command_line(argc, argv, logger);

    srand(NmaOptions::getInstance().seed.value());   

//     for(int i=1; i < 2; i++)
//       {
// 	LOGOUT(rand());
//         LOGOUT(normrnd(1,1,0,1).AsScalar());
//       }
//     exit(0);

    if(opts.debuglevel.value()==1)
      Tracer_Plus::setrunningstackon();

    if(opts.timingon.value())
      Tracer_Plus::settimingon();

    ////////////////////////////////////
    //////////// bivariate normal test 

//     vector<string> parameter_name(2); parameter_name[0]=string("x1"); parameter_name[1]=string("x2");
//     vector<float> initial_values(2); initial_values[0]=0.4; initial_values[1]=1;
//     vector<float> param_min(2); param_min[0]=0;param_min[1]=-100;
//     vector<float> param_max(2); param_max[0]=100;param_max[1]=100;

//     float x1=0;
//     float x2=0;

//     OUT("HERE1");

//     Bvn_Mcmc_Component bvn(x1,x2,"bvn", parameter_name, initial_values, param_min, param_max, opts.debuglevel.value(), 1);

//     OUT("HERE2");

//     vector<Mcmc_Component*> components;
//     components.push_back(&bvn);
//     Bvn_Mcmc_Log_Likelihood bvn_mcmc_log_likelihood(x1,x2);

//     OUT("HERE3");

//     Mcmc_hmc mcmc_hmc(components,bvn_mcmc_log_likelihood,opts.nsamps.value(),opts.burnin.value(),opts.sampleevery.value(),opts.debuglevel.value(),1);

//     OUT("HERE4");

//     mcmc_hmc.run();

//     OUT("HERE5");

//     mcmc_hmc.save();

//     exit(1);

    //////////// end of bivariate normal test 
    ////////////////////////////////////
  
    Nma_manager nma_manager;

    LogSingleton::getInstance().str() << "Setting up:" << endl;
    nma_manager.setup();
    LogSingleton::getInstance().str() << "Running:" << endl;
    nma_manager.run();

    if(opts.timingon.value())
      Tracer_Plus::dump_times(logger.getDir());

    LogSingleton::getInstance().str() << endl << "Log directory was: " << logger.getDir() << endl;

    LogSingleton::getInstance().str() << endl << "For a web page report view: " << endl << LogSingleton::getInstance().appendDir("report.html")  << endl;
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












