# OTF - N-Dimensional Optimization Testing Functions

This library contains Fortran implementations of n-dimensional functions
outlined [here](https://www.sfu.ca/~ssurjano/optimization.html). Some 2d-only
functions are also implemented, so beware of those.

Each implementation has a matrix and array input inplementation. These allow
for the evaluation of multiple points at once or a single point respectively.
Matrix implementations are meant to maximize gains from vectorization and
should be used when possible.

## Fortran Package Manager (fpm) Setup

This implementation can be easily used as a dependency using the [Fortran
Package Manager](https://github.com/fortran-lang/fpm) as follows.

Using OTF as a dependency of the main program:

```toml
[dependencies]
OTF = { git ="https://github.com/renatomatz/otf.git" }
```

Using OTF as a dependency of a specific test:

```toml
[[ test ]]
name="test"
source-dir="test/"
main="test.f90"
[test.dependencies]
OTF = { git ="https://github.com/renatomatz/otf.git" }
```

The aboce example works for different test names, soruces etc. depending on
your specific program structure. For more information on using fpm, check out
their [packaging
guide](https://github.com/fortran-lang/fpm/blob/master/PACKAGING.md)

## Example Visualizations

You can visualize the implemented functions using
[gnuplot](http://www.gnuplot.info/) by running the following script:

```shell
fpm test
```
