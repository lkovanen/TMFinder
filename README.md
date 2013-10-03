TMFinder
========

Temporal Motif Finder can be used to detect and enumerate motifs in temporal networks. For more information, including the formal definition of temporal networks and temporal motifs as well as the algorithms used, see

> Lauri Kovanen, Márton Karsai, Kimmo Kaski, János Kertész, Jari
Saramäki, "Temporal motifs in time-dependent networks",
Journal of Statistical Mechanics: Theory and Experiments. P11005 (2011)
doi:10.1088/1742-5468/2011/11/P11005

The code also implements the null model for identifying differences between different node and event types. This null model was introduced in

> Lauri Kovanen, Kimmo Kaski, János Kertész, Jari Saramäki. "Temporal motifs reveal homophily, gender-specific patterns and group talk in mobile communication networks." arXiv:1302.2563 (2013).

If you use this code in scientific work, please cite the above publication.

TMFinder is published under the GPL v3 licence.

Installation
------------

TMFinder uses [bliss][bliss] to calculate canonical forms of graphs. You need to install bliss separately to use TMFinder:
   
1. [Download bliss][bliss] (version 0.72 should work) and compile it by
	following the instructions included with it.

2. Add the path of the bliss directory to the environment variables
	`LIBRARY_PATH` and `CPLUS_INCLUDE_PATH`.

You should now be able to compile TMFinder by calling `make` in the directory `src`. If you get an error message about `bliss` or `graph.hh`, recheck your installation of bliss and the environment variables pointing to the location of the bliss library.

After the compiling, make sure everything works by running the test script `tests/test_small.sh`. This should produce a single output file, `test_small_output.dat` that contains information about the temporal motifs in the small test data.

Python code for handling temporal motifs
----------------------------------------

The `python` subdirectory contains python scripts that should prove useful for analysing and visualizing temporal motifs in the output file. The code relies on PyBliss for calculating graph isomorphisms.

1. [Download and install PyBliss][bliss], the Python wrapper for bliss.
  	PyBliss uses version 0.50 of bliss (this comes with the PyBliss
 	package).

2. Add the path to your PyBliss directory to the environment
   variable PYTHONPATH.

The scripts also require a number of other Python libraries that also need to be installed: `pylab`, `numpy`, `pygraphviz` and `argparse`. Make sure these are all installed.

Finally, to make sure everything works run `tests/test_plotting.sh`. This script reads the output of `tests/test_small.sh` and visualizes the detected motifs.

Usage instructions
------------------

For documentation about input and output file formats and usage options, call `bin/tmf --help`. The test scripts should also provide an example for getting started.


[bliss]: http://www.tcs.hut.fi/Software/bliss/ "bliss"
