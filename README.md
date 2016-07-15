* PolyMage

[![Code 
Climate](https://codeclimate.com/github/bollu/polymage/badges/gpa.svg)](https://codeclimate.com/github/bollu/polymage)

[![Build 
Status](https://travis-ci.org/bollu/polymage.svg?branch=github-integration)](https://travis-ci.org/bollu/polymage)

PolyMage is a domain-specific language and optimizing code generator for 
automatic optimization of image processing pipelines, being developed at the 
Multicore Computing Lab, Indian Institute of Science. PolyMage takes an image 
processing pipeline expressed by the user in a high-level language (embedded in 
Python) and generates an optimized parallelized C++ implementation of the 
pipeline.

* INSTALLATION GUIDE

PolyMage is a pure-python library that requires a C++ compiler for code
generation.


** REQUIREMENTS

1) Python 3.x

2) Python packages numpy, pytest. These can be installed via  
(on a Fedora) $ sudo yum -y install python3-numpy python3-pytest  
(on Ubuntu) $ sudo apt-get install python3-numpy python3-pytest  

3) OpenCV 2.4.7 or higher (with QT/GTK support, video codec support for the video demo),  
Python bindings for OpenCV. Install instructions on Ubuntu: https://help.ubuntu.com/community/OpenCV  
If you don't have a GPU on your machine, be sure to call cmake with the option -D WITH_CUDA=OFF  
On a Fedora, these can be installed with 'sudo yum -y install opencv python-opencv'

4) g++ (GNU C++ compiler) version 4.8 or higher or Intel C/C++ compiler (icpc) 12.0 or higher  
(recommended: icpc 14.0 or higher)

5) Python bindings for isl  
islpy http://documen.tician.de/islpy/  
This can be easily installed via python3-pip  
```
$ sudo yum -y install python3-pip  
$ sudo pip3 install islpy  
```

(islpy itself requires ffi development files -- this can be installed by 
installing libffi-devel via yum/apt-get)


** INSTALLATION

```
$ git clone git@bitbucket.org:udayb/polymage.git

$ cd polymage

$ git submodule update --init

$ cd cgen

$ git am ../patches/0001-ctye-to-dtype-handle-void.patch

$ cd ..
```

** PROJECT STRUCTURE


`polymage`  is the main directory of interest and it contains most of the code.  


`polymage/tests` is the test directory and has a lot of sample code which you can take a look at.  
You can run the tests by invoking the following command:  
```
$ py.test-3 test_{name}.py  
```

For example, the harris corner detection test can be run using the following command from the  
`polymage/tests` directory

```
polymage/polymage/tests$ py.test-3 test_harris.py
```
Note: The input language does not exactly match that in the ASPLOS 2015 paper; however,  
it is very close.  

`python/apps/polymage` : has some benchmark applications written using a python driver code.  
Here we eliminated the need for a C++ driver and manage the pipeline input and output in python.  

`polymage/constructs.py` : PolyMage language constructs  

`polymage/pipe.py` : high level flow of the optimizer  

`polymage/poly.py` : for polyhedral representation of the pipelines  

`polymage/schedule.py` : schedule transformation for the computations  

`polymage/codegen.py` : code generation for the scheduled pipeline  

`polymage/targetc.py` : c++ code generation  

`polymage/tuner.py` : autotuning code  

`polymage/libpluto.py` : FFI access to PLUTO


The following repository contains just the base and the best PolyMage optimized codes (for Intel  
Sandybridge) used for experiments in the ASPLOS 2015 paper for all of the benchmarks -- these are  
sufficient if one is purely interested in a final performance comparison without any tweaking/tuning:  
`https://github.com/bondhugula/polymage-benchmarks`


** LICENSE

PolyMage is available under the Apache License, version 2.0. Please see 
the LICENSE file for details.
