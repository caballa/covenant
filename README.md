# Covenant #

Intersection of Context-Free Languages

#About#

Covenant is a tool for testing emptiness of a set of context free
languages (CFLs). Since this problem is undecidable, Covenant
implements a CEGAR-based schema which might not terminate. Covenant
implements several refinement techiques. One of these refinements is
complete if the CFLs are regularly separable.

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

#Usage#

`build/tools/covenant --help` 

The format of the CFLs is a bit inconvenient for now.  This is an example:

    ; this is a comment

    (  S1 -> [ A B]; 
      A  -> [ \x61 \x61]; 
      A  -> [ \x62 \x62]; 
      A  -> [ \x61 S1 \x61, \x62 S1 \x62 ]; 
      B  -> [ \x61 \x62 B]; 
      B  -> [ \x61 \x62]  
    )
    
    (  S2 -> [ A B]; 
      A  -> [ \x61 \x61]; 
      A  -> [ \x62 \x62]; 
      A  -> [ \x61 S2 \x61, \x62 S2 \x62 ]; 
      B  -> [ \x62 \x61 B];
      B  -> [ \x62 \x61]  
    )  

The non-terminal symbol appearing on the left-hand side in the first
grammar production is considered the start symbol of the grammar. In
the above example the start symbols are `S1` and `S2`, respectively.

A current limitation is that terminal symbols can only ranges from 0
to 127 and they must be described in hexadecimal format. E.g., `\x61` is
the ASCII character "a", `\x62` is "b", and so on.

#People#

* [Jorge A. Navas](http://ti.arc.nasa.gov/profile/jorge/)
* Graeme Gange
* Peter Schachte
* Harald Sondergaard
* Peter J. Stuckey