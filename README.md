
PolyMage is release under the Apache License, version 2.0. Please see the 
LICENSE file for more details.

**REQUIREMENTS**

1) Python 3.x, numpy, pytest

2) OpenCV 2.4.7 or higher (with QT/GTK support, video codec support for the video demo),  
Python bindings for OpenCV. Install instructions on Ubuntu: https://help.ubuntu.com/community/OpenCV  
If you don't have a GPU on your machine, be sure to call cmake with the option -D WITH_CUDA=OFF  
On a Fedora, these can be installed with 'sudo yum -y install opencv python-opencv'

3) g++ (GNU C++ compiler) version 4.8 or higher or Intel C/C++ compiler (icpc) 12.0 or higher  
(recommended: icpc 14.0 or higher)

4) Python bindings for isl  
islpy http://documen.tician.de/islpy/  
This can be easily installed via python3-pip  
$ sudo yum -y install python3-pip  
$ sudo pip3 install islpy  
C code generation library:  
cgen https://github.com/inducer/cgen.git  

(islpy itself requires ffi development files -- this can be installed by 
installing libffi-devel via yum/apt-get)

5) Python packages numpy, pytest. These can be installed via  
(on a Fedora) $ sudo yum -y install python3-numpy python3-pytest  
(on Ubuntu) $ sudo apt-get install python3-numpy python3-pytest  

**PROJECT STRUCTURE**
sandbox is the main directory of interest and it contains most of the code.  
sandbox/tests is the test directory and has a lot of sample code which you can take a look at.  
You can run the tests by invoking the following command:  
$> py.test-3 test_{name}.py  
For example, the harris corner detection test can be run using the following command from the  
sandbox/tests directory:
.../sandbox/tests$> py.test-3 test_harris.py

Note: The input language does not exactly match that in the ASPLOS 2015 paper. However,  
it is very close.  

sandbox/apps/python : has some benchmark applications written using a python driver code.  
Here we eliminated the need for a C++ driver and manage the pipeline input and output in python.  

sandbox/constructs.py : PolyMage language constructs  

sandbox/pipe.py : high level flow of the optimizer  

sandbox/poly.py : for polyhedral representation of the pipelines  

sandbox/schedule.py : schedule transformation for the computations  

sandbox/codegen.py : code generation for the scheduled pipeline  

sandbox/targetc.py : c++ code generation  

sandbox/tuner.py : autotuning code  

The following repository contains just the base and the best PolyMage optimized codes (for Intel  
Sandybridge) used for experiments in the ASPLOS 2015 paper for all of the benchmarks -- these are  
sufficient if one is purely interested in a final performance comparison without any tweaking/tuning:  
https://github.com/bondhugula/polymage-benchmarks
