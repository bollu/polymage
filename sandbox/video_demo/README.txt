Compiling the cpp file as a shared library
==========================================
icpc -xhost -openmp -fPIC -shared -o <file>.so <file>.cpp

* Try queuing frames using one thread and processing using the rest
* Add a explanation on how the demo works
* Include the bilateral grid that works for all sizes
  currently the pipeline is built for a specific size.
