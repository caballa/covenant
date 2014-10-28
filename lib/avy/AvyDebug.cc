//#include "llvm/Support/CommandLine.h"
#include <boost/program_options.hpp>

#include <avy/AvyDebug.h>

#include <string>
#include <set>

#ifndef NAVYLOG
using namespace avy;

bool avy::AvyLogFlag = false;
std::set<std::string> avy::AvyLog;

void avy::AvyEnableLog (std::string x) 
{
  if (x.empty ()) return;
  AvyLogFlag = true;
  AvyLog.insert (x); 
}

namespace avy
{
  struct LogOpt
  { void operator= (const std::string &tag) const 
    { 
      avy::AvyEnableLog (tag); 
    } 
  };
  
  LogOpt loc;
}

// namespace{
//   namespace po = boost::program_options;
//   struct LogOpts { static po::options_description opts; };

// // po::options_description log_opts("Logging Options");
// // log_opts.add_options()
// //     ("log", boost::program_options::value<std::string>(), "Enable specified log level")
//   // ;

// po::options_description LogOpts::opts("Logging Options");
// po::options_description_easy_init const dummy = opts.add_options()
//      ("log", boost::program_options::value<std::string>(), "Enable specified log level")
// ;
// // po::options_description_easy_init const dummy = cmd_opt.add_options()
// // //     ("log", boost::program_options::value<std::string>(), "Enable specified log level")
// // //     ;
// }

// static llvm::cl::opt<avy::LogOpt, true, llvm::cl::parser<std::string> > 
// LogClOption ("log",
//              llvm::cl::desc ("Enable specified log level"),
//              llvm::cl::location (avy::loc),
//              llvm::cl::value_desc ("string"),
//              llvm::cl::ValueRequired, llvm::cl::ZeroOrMore);

#else
#endif

