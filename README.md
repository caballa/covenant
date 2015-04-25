# Covenant #

Intersection of Context-Free Languages

#About#

Covenant is a tool for testing emptiness of a set of context free
languages (CFLs). Since this problem is undecidable, Covenant
implements a CEGAR-based schema which might not terminate. Covenant
implements several refinement techniques. One of these refinements is
complete if the CFLs are regularly separable.

Read this [technical report](http://arxiv.org/abs/1411.5131) for details.

#Prerequisites#

- Boost is required. It is recommended version 1.55 or newer
- Be sure to install the Boost program_options library and set `BOOST_ROOT`

#Installation#

- `mkdir build && cd build`
- `cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=my_install_dir ..`
- `cmake --build . --target install`

The code has been tested only for X86_64 with clang++ 3.2 and g++ 4.8

#Usage#

`my_install_dir/bin/covenant --help` 

We describe the format of the CFLs through an example:

    ; this is a comment

    (  S1 -> [ A B]; 
      A  -> [ "a" "b"]; 
      A  -> [ "b" "b"]; 
      A  -> [ "a" S1 "a", "b" S1 "b" ]; 
      B  -> [ "a" "b" B]; 
      B  -> [ "a" "b"]  
    )
    
    (  S2 -> [ A B]; 
      A  -> [ "a" "a"]; 
      A  -> [ "b" "b"]; 
      A  -> [ "a" S2 "a", "b" S2 "b" ]; 
      B  -> [ "b" "a" B];
      B  -> [ "b" "a"]  
    )  

The non-terminal symbol appearing on the left-hand side in the first
grammar production is considered the start symbol of the grammar. In
the above example the start symbols are `S1` and `S2`, respectively.

Terminal symbols must be between double quotes (i.e.,
`""`). Any symbol that is not between double quotes will be
interpreted as a non-terminal symbol.

Each grammar production is of the form `A -> [ ... ];` where `A` is a nonterminal 
and `...` is a sequence of any nonterminal or terminal symbol separated by one or 
more blanks. We also allow `A -> [ ... , ... ];` to express two grammar
productions with the same left-hand side A. That is,  `A  -> [ "a" S2 "a", "b" S2 "b" ];` 
is syntactic sugar for `A  -> [ "a" S2 "a"]; A ->[ "b" S2 "b" ];`
 
Note that all the right-hand side of the productions must ends up with the symbol `;`
except the last one.

If we wrap the above example into a file test.cfg and try
`covenant test.cfg`, you should obtain:

`Finished after 5 cegar iterations.`   

` UNSAT`

Since covenant might not terminate, we provide a script `covenant-par`
that runs in parallel several configurations (i.e., heuristics) and
stops as soon as one of the them terminates:

`my_install_dir/bin/covenant-par file`

#People#

* [Jorge A. Navas](http://ti.arc.nasa.gov/profile/jorge/)
* Graeme Gange
* Peter Schachte
* Harald Sondergaard
* Peter J. Stuckey

#Publications#

- "A Tool for Intersecting Context-Free Grammars and Its Applications". G.Gange , J.A.Navas, P.Schachte, H.Sondergaard, and P.J. Stuckey. [(PDF)](http://www.clip.dia.fi.upm.es/~jorge/docs/cfg_nfm15.pdf) . NFM'15
