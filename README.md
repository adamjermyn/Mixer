Mixer is a software package for computing turbulent transport coefficients. The package is composed of a C++ backend which performs the calculations and an optional Python frontend.

# Building Mixer

In the `C++` directory run `make shared` to make the Python interface along with all of its requirements. The makefile assumes that the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [Cubature](https://github.com/stevengj/cubature) libraries are two directories above the `C++` directory (e.g. one above the root of Mixer).

Mixer is known to compile correctly on both Mac and Linux systems with Eigen v3.3.3 and Cubature v1.0.2, though other versions likely work as well. The Python interface makes use of Matplotlib and Numpy as well as H5Py.

# Example Usage

Please see the `Plot Scripts` directory for examples. Script named with the word 'fit' interface with Mixer and calculate transport coefficients, while those with the word 'plot' generate plots like those found in the methods paper.

# Data

The methods paper provides the results of calculations done with Mixer. The raw data underlying those results may be found in the `Data` directory. Data files are stored in the HDF5 format.

# Referencing Mixer

Mixer is free to use, but if you use it for academic purposes please include a citation to the methods paper:


The BibTex entry for this is included below:

```

@ARTICLE{Mixer2018,
   author = {{Jermyn}, A.~S. and {Lesaffre}, P. and {Tout}, C.~A and {Chitre}, S.~M.},
    title = "{Turbulence Closure for Mixing Length Theories }",
  journal = {\mnras},
 primaryClass = "astro-ph.SR",
 keywords = {convection, stars: evolution, stars: interiors, stars: rotation},
     year = 2018,
    month = jan
}


```

