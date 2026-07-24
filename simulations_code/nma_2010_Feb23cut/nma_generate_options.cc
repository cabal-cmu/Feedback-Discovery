/*  nma_generate_options.cc

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  CCOPYRIGHT  */

#define WANT_STREAM
#define WANT_MATH

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "nma_generate_options.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"

using namespace Utilities;

namespace Nma {

NmaGenerateOptions* NmaGenerateOptions::gopt = NULL;

void NmaGenerateOptions::parse_command_line(int argc, char** argv, Log& logger)
{
  Tracer_Plus("NmaGenerateOptions::parse_command_line");

  // do once to establish log directory name
  for(int a = options.parse_command_line(argc, argv); a < argc; a++);
    
  // setup logger directory
  logger.makeDir(logdir.value(),"logfile",true,true);
  
  logger.str() << "Log directory is: " << logger.getDir() << endl;

  // do again so that options are logged
  for(int a = 0; a < argc; a++)
    logger.str() << argv[a] << " ";
  logger.str() << endl << "---------------------------------------------" << endl << endl;

  if(help.value() || ! options.check_compulsory_arguments())
    {
      options.usage();
      throw Exception("Not all of the compulsory arguments have been provided");
    }      
}

}
