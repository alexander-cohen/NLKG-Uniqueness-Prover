# ${1:Uniqueness of excited states to $-\Delta u + u - u^3 = 0$ in three dimensions}
This project consists of code accompanying the paper https://arxiv.org/abs/2101.08356. This is part of joint work done by Alex Cohen, Zhenhao Li, and Wilhelm Schlag during the summer and fall of 2020. 

## Installation
This software relies on VNODE-LP, a validated ODE solver written by Nedialko S. Nedialkov. See https://www.cas.mcmaster.ca/~nedialk/vnodelp/doc/vnode.pdf for the documentation and instructions on installation. Once VNODE-LP is succesfully installed, the code in this project may be build and run, possibly with small modifications to the makefile. 

We note that VNODE-LP can only be built on certain computer architectures, so on many computers a virtual machine will be necessary. 

## Usage
The entire code for this project is contained in the C++ file 'nlkg_uniqueness_prover.cc'. To run the code as is, build the project by running make, and then run './solver N' to prove uniqueness of the first 'N' excited states. 

Some other functionalities are available, e.g. for outputting data from VNODE for graphing purposes. See the C++ file for documentation.

Feel free to reach out to the authors with any questions. Alex Cohen can be reached at alex.cohen@yale.edu.
