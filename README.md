# Covenant #

Intersection of Context-Free Languages

#About#

Covenant is a tool for testing emptiness of a set of context free
languages (CFLs). Since this problem is undecidable, Covenant
implements a CEGAR-based schema which might not terminate. Covenant
implements several refinement techiques. One of these refinements is
complete if the CFLs are regularly separable.

Read this [technical report](http://people.eng.unimelb.edu.au/gkgange/pubs/cfg_preprint.pdf) for details.

#Prerequisites#

- Boost is required. It is recommended version 1.55 or newer
- Be sure to install the Boost program_options library and set `BOOST_ROOT`

#Installation#

- `mkdir build && cd build`
- `cmake ..`
- `make`

If ninja is installed then try instead:

- `mkdir build && cd build`
- `cmake -G Ninja ..`
- `ninja`

The code has been tested only for x86_64 with clang++ 3.2 and g++ 4.8

#Usage#

`build/tools/covenant --help` 

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
`build/tools/covenant test.cfg`, you should obtain:

`Finished after 5 cegar iterations.`   

` UNSAT`

#People#

* [Jorge A. Navas](http://ti.arc.nasa.gov/profile/jorge/)
* Graeme Gange
* Peter Schachte
* Harald Sondergaard
* Peter J. Stuckey