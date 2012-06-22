Introduction to Bioinformatics, Summer 2012.
Programing Project 2
############################################

Operating System
----------------

The project was developed in the Java language, on Ubuntu Linux 12.04, Java version: OpenJdk 1.6.0_24

The content
------------

The 'src' directory contains:
1. Java source code files.
2. Sample sequences in .fa file.
3. The run configuration file, 'run.config'

How to Run
----------

1. Go to the 'src' directory.
2. Compile the java files: "javac *.java"
3. Run the Aline.class file: "java DiscoverMotifs"
** The output is written to a text file with the name 'output.txt' instead of printing to the screen, as the results can be too long.

How to change run configuration
-------------------------------

1. Go to the 'src' directory
2. open 'run.config'
3. Do not remove any lines. Just change the values on the right side of the colon ':'.
4. Only one single input text file is required. The single file should contain the sequences in FASTA format.


Google code project url:
------------------------
http://code.google.com/p/expectation-maximization/


Output files:
-------------
The output files contain the results for each initial random alignment and then finally the maximum log2odd scores over all the iterations of all the initial random alignments, see the headings per each table in the output text files. 
1. 'output6.txt' file contains the output for motif length of 6, initial random count of 50, iterations per each initial random of 500.
2. 'output4.txt' file contains the output for motif length of 4, initial random count of 50, iterations per each initial random of 500.

