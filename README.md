# Eucalc

**Authors:** H. Passe, V. Lebovici

`Eucalc` is a `C++`\ `Python` library that implements algorithms to efficiently compute topological integral transforms (Euler characteristic transforms, Radon transforms and hybrid transforms) on cubical complexes. It can be used as a header library for `C++` users, or can be compiled into a `Python` module.

This repository is associated with [our paper](https://arxiv.org/abs/2405.02256) with S. Oudot accepted in *Symposium on Experimental Algorithms (SEA) 2024*.

Please find short usage tutorials in the notebooks of the `tutos` directory.

## Python module compilation

To use our `C++` code in `Python` we use `Cython`. To compile the module, we included the necessary headers from the [`GUDHI`](https://gudhi.inria.fr) library in `src/gudhi`. 

To compile `Eucalc`, open the main directory into a terminal and execute:

```
python3 setup.py build_ext --inplace
```

Then copy-paste the file `eucalc.XXX.so` to the directory you are working in or start your `Python` file or Jupyter notebook with:
```
import sys
sys.path.append('eucalc_directory')
import eucalc
```
where `eucalc_directory` is the main directory of this repository. You can start using `Eucalc`!