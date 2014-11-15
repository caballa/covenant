#include <cstdio>
#include <cstdlib>

#include <boost/program_options.hpp>

#include <Common.h>
#include <parser/ParseUtils.h>
#include <parser/ParseProblem.h>
#include <CFG.h>
#include <Solver.h>
#include <regsolver/Product.h>

using namespace std;
using namespace covenant;

namespace po = boost::program_options;

int main (int argc, char** argv){

  typedef Sym edge_sym_t;
  typedef Solver<edge_sym_t> solver_t;
  typedef RegSolver<edge_sym_t> reg_solver_t;

  // Main flags
  int max_cegar_iter;
  GeneralizationMethod gen;
  AbstractMethod abs;
  // Heuristics flags
  // - regular solver will try to find a witness starting from this length
  int shortest_witness;
  // - increase the witness length after n iterations -1 disable this option.
  int freq_incr_witness;
  // - increment the witness length by n but only if freq_incr_witness >= 0
  unsigned incr_witness;
  // Other flags
  bool is_dot_enabled=false;
  unsigned num_solutions;

  string header("Covenant: semi-decider for intersection of context-free languages\n");
  header += string ("Authors : G.Gange, J.A.Navas, P.Schachte, H.Sondergaard, and P.J.Stuckey\n");
  header += string ("Usage   : covenant [Options] file");

  po::options_description config("General options ");

  config.add_options()
      ("input-file,f",  po::value<string>(), "input file")
      ("help,h", "print help message")
      ("dot,d" , 
       "print solutions and emptiness proofs in dot format as well as the result of abstractions and refinements")
      ("verbose,v" , "verbose mode")
      ("stats,s", "Show some statistics and exit")
      ("solutions" , po::value<unsigned>(&num_solutions)->default_value(1), 
       "enumerate up to n solutions (default n=1)")
      ("iter,i",  po::value<int>(&max_cegar_iter)->default_value(-1), 
       "maximum number of CEGAR iterations (default no limit)")
      ("abs,a",  po::value<AbstractMethod>(&abs)->default_value(CYCLE_BREAKING), 
       "choose abstraction method [sigma-star|cycle-breaking]")
      ("gen,g",  po::value<GeneralizationMethod>(&gen)->default_value(GREEDY), 
       "choose generalization method [greedy|max-gen]")
      ;
 
  po::options_description cegar_config("(unsound) cegar options ");

  cegar_config.add_options()
      ("l",  po::value<int>(&shortest_witness)->default_value(0), 
       "shortest length of the witness (default 0)")
      ("freq-incr-witness",  po::value<int>(&freq_incr_witness)->default_value(-1), 
       "how often increasing witnesses' length")
      ("delta-incr-witness",  po::value<unsigned>(&incr_witness)->default_value(1), 
       "how much witnesses' length is incremented")
      ;

  po::options_description log("Logging Options");
  log.add_options()
      ("log",  po::value<std::vector<string> >(), "Enable specified log level")
      ;

  po::options_description cmmdline_options;
  cmmdline_options.add(config).add(log);
  cmmdline_options.add(cegar_config);

  po::positional_options_description p;
  p.add("input-file", -1);
  
  po::variables_map vm;

  try {
    po::store(po::command_line_parser(argc, argv).
              options(cmmdline_options).
              positional(p).
              run(), vm);
    po::notify(vm);    
  }
  catch(error &e)
  {
    cerr << "covenant error:" << e << endl;
    return 0;
  }

  if (vm.count("help"))
  {
    cout << header << endl << config << endl; 
    return 0;
  }  

  std::string in;
  bool file_opened = false;
  if (vm.count ("input-file"))
  {
    std::string infile = vm ["input-file"].as<std::string> ();
    std::ifstream fd;
    fd.open (infile.c_str ());
    if (fd.good ())
    {
      file_opened = true;
      while (!fd.eof ())
      {
        std::string line;
        getline(fd, line);
        in += line;
        in += "\n";
      }
      fd.close ();
    }
  }

  if (!file_opened)
  {
    cout << "Input file not found " << endl;
    return 0;
  }

  if (vm.count("dot"))
  {
    is_dot_enabled = true;
  }

  if (vm.count("verbose"))
  {
    avy::AvyEnableLog ("verbose");
  }

  // enable loggers
  if (vm.count("log"))
  {
    vector<string> loggers = vm ["log"].as<vector<string> > ();
    for(unsigned int i=0; i<loggers.size (); i++)
    {
      avy::AvyEnableLog (loggers [i]);
    }
  }

  solver_t::Options opts(max_cegar_iter, 
                         gen, 
                         abs, 
                         is_dot_enabled, 
                         shortest_witness, 
                         freq_incr_witness, 
                         incr_witness,
                         num_solutions);

  CFGProblem problem;
  //StreamParse input(cin);
  StrParse input(in);
  boost::shared_ptr<TerminalFactory> tfac (new TerminalFactory ());

  try
  {
    parse_problem(problem, input, tfac);
  }
  catch(error &e)
  {
    cerr << "covenant error:" << e << endl;
    return 0;
  }

  if (problem.size() == 0)
  {
    cerr << "covenant: no CFGs found." << endl;  
    return 0;
  }

  try
  {
    for(int ii = 0; ii < problem.size(); ii++)
    {
      CFGProblem::Constraint cn(problem[ii]);
      if (cn.vars.empty())
      {
        throw error("expected at least one CFG");
      }
      cn.lang.check_well_formed();
    }
  }
  catch(error &e)
  {
    cerr << "covenant error:" << e << endl;
    return 0;
  }
    
  reg_solver_t *regsolver = new Product<edge_sym_t>();
  try
  {
    solver_t s(problem, regsolver, opts);

    s.preprocess ();

    if (vm.count("stats"))
    {
      s.stats (cout);
      delete regsolver;
      return 0;
    }

    STATUS res = s.solve();

    switch (res)
    {
      case SAT:
        cout << "======\nSAT\n======\n";
        break;
      case UNSAT:
        cout << "======\nUNSAT\n======\n";
        break;
      default:
        cout << "======\nUNKNOWN\n======\n";  
        break;
    }
      
  }
  catch(error &e)
  {
    cerr << "covenant error:" << e << endl;
  }
  catch(Exit &e)
  {
    cerr << e << endl;
  }

  delete regsolver;

  return 0;
}


