# AutoMine

Stage 1: Pattern Generation
The compiler in /compiler can generate motifs and cliques of a particular size. This is done using the binary /compiler/bin/compile with the flag --motif <pattern size> or --clique <pattern size>. This first generates the relevant patterns as .pat files in /patterns, then proceeds to Stage 2 automatically.

Stage 2: Pattern Compiler
The compiler in /compiler can either begin scheduling with all the patterns from Stage 1, or with a particular pattern from a .pat file. This is done using the binary /compiler/bin/compile with the flag --pattern <pattern file .pat>. This generates the schedule for counting the relevant pattern(s) as c++ code in a .hpp file named plan_<name>.hpp in /plans.

Stage 3: Host Compiler
The Makefile in / automatically generates binaries in /bin named <name> for each schedule in /plans when running `make` from /

Stage 4: Counting
The binaries in /bin can be run with a flag <prefix> where the files <prefix>.meta.txt, <prefix>.vertex.bin, and <prefix>.edge.bin. The prefixes ending with .ord correspond to graph files where the vertices have been ordered according to the optimization in the GraphZero paper.
