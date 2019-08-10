Compile: c++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes` mormon.cpp -o mormon`python3-config --extension-suffix`

Use Wurlitzer to show output in jupyter notebook cells

Python bindings with pybind11