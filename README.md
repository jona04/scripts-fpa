This repository contains all the code implementing the Fixed Point method to solve the Quadratic Convex Separable Knapsack (QSKP) problem, as described in the manuscript

 * A. Alves, J.O.L. Silva, L.C. Matioli, P.S.M. Santos, S.S. Souza. "A fixed-point algorithm for solving quadratic convex separable knapsack problems".

## Important Note

This original source code which this project was based can be found in this page:
http://www.ime.unicamp.br/~pjssilva/research/quadratic_knapsack_source.zip.

The main difference between the original source code and code  presented here is the addiction of the new algorithm Fixed point (FPA) and the removal of SVM experiments, since in this paper we don't present this kind of experiments. 

If you run into problems related to Fixed point algorithm (FPA) do not hesitate to  contact me:

* Jonatas O. L. Silva <jonatas.iw@gmail.com>

For any other doubts you can contact directly the original author:

* Paulo J. S. Silva <pjssilva@ime.unicamp.br>

All the code is licensed under the GPL Version 2 or later, this includes the files in the third_party directory that were originally licensed as such. You can find a full copy of the license in the COPYING file in the current directory.

All original references and comments presented in the original code have been kept.

## Directory organization

* third_party/ : Directory containing the code necessary to implement some third party methods

* root_finding_algorithms/lib/ : Directory containing the implementation of the root finding solvers applied to QSKP problem. It has four files:
    * bissection_method.c (and .h)
    * fixed_point.c (and .h)
    * regula_falsi.c (and .h)
    * secant_method.c (and .h)

* state_of_the_art_algorithms/lib/ : Directory containing the implementation of the state of the art solvers applied to QSKP problem. It has three files:
    * third_party_methods.c (and .h)
    * cont_quad_knapsack.c (and .h)
    * fixed_point.c (and .h)

*  state_of_the_art_algorithms/lib/optmize.py : Run this to generate an optimized version of the methods above
    when the Hessian diagonal is a fixed value (that is, the d vector is a
    constant) or the hyperplane normal is a fixed value (that is, the b vector
    is a constant). You can also use it to turn off compression (variable
    fixing) in Newton method when the problem is expected to be easy (require
    just a few Newton steps). Run the program with --help option to see usage
    information.

## Tests and results
* root_find_algorithms/tests/ and state_of_the_art_algorithms/tests/ : Directory containing the test problems used in the numerical section of the manuscript. To run all test just type in the command line

    
    ``` python run_tests.py``` 


This will run the two set of tests. They can take many hours to complete.

After run the command from the directory root_find_algorithms/tests/, a .tex file inside of test_results/ will be created contained the tables according to Section 4.1 of the manuscript. 

After run the command from the directory state_of_the_art_algorithms/tests/, a .tex file inside of test_results/ will be created contained the tables according to Section 4.2 of the manuscript.

## Compilation 

As already mentioned, to run the tests and the algorithms, you only need to run python run_tests.py.

The file Makefile will then generate all necessary compilation files.

This code was implemented in standard C (version C99) and it can be compiled and used by any standard compliant compiler.

To compile and run the tests you need a unix like shell, make, gcc, gfortran
version 4.4 or later, python (2.6 <= version < 3.9) with scipy and numpy
packages. The code assumes that you also have a BLAS library installed and that
it can be linked with a simple flag like -lblas. IF you do not have BLAS you
will need to edit the files tests/synthetic/Makefile and tests/svm/Makefile and
comment the lines

## Some configs

------------

CFLAGS += -DLAPACK

LIBS += -lblas

------------

You can also change the compilation flags to link your version of BLAS by
editing the second line above.

In a Linux distribution you can install the compiler and python dependecies above easily, example in Ubuntu 12.04 we can try:

sudo apt-get install gcc gfortran make python-numpy python-scipy

In windows I suggest you to install mingw compiler with latest version of
mingw-get-inst from

http://sourceforge.net/projects/mingw/files/Installer/mingw-get-inst/

Run the installer and select both the C and Fortran compilers and
"MSYS Basic System" option (that will give you a unix like
shell). To install Python is suggested to use Python(x,y)

http://code.google.com/p/pythonxy/

I already comes with numpy and scipy.

The installation of BLAS can be a little tricky and depends on the platform.
For Ubuntu I would suggest

`sudo apt-get install libopenblas-base libopenblas-dev`
