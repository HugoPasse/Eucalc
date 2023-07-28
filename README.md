# Eucalc

Eucalc is a C++ library that implements algorithms to compute hybrid transforms, Radon transform and Euler characteristic transform on cubical complexes. It can be used as an header library for C++ users, or can be compiled into a python module.

## Python module compilation

To use our C++ code in python we use cython. To compile the module you will need the  [`Gudhi`](https://gudhi.inria.fr) library headers. You can find them by downloading the [`latest release of Gudhi`](https://gudhi.inria.fr/release/Gudhi-Release-3.8.0/). Once you've downloaded them, modify the file called `setup.py` :

```python
libs = ["/your_path_to_gudhi"]
```
 Then you can compile the module with :
```
python3 setup.py build_ext --inplace
```

Don't forget to add `~/your_module.so` to your python path (use `export PYTHONPATH=~/your_module.so`) or to put the `your_module.so` file into you python's module directory.
