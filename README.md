# Susa Open Source Project

![example workflow](https://github.com/libsusa/susa/actions/workflows/github-actions.yml/badge.svg) [![codecov][codecov-badge]][codecov-link]

Susa is a mathematics and signal processing C++ [framework](https://en.wikipedia.org/wiki/Software_framework) based on [KISS](https://en.wikipedia.org/wiki/KISS_principle)
 principle. It is stand-alone with a modern architecture. It is designed not to have any dependencies to none standard third
 party libraries. Indeed, a C++17 compiler along with STL is necessary and sufficient in order to compile it. Therefore,
 portability is the key feature of Susa. For example it can be exploited in mobile platforms such as Android NDK (Native
 Development Toolkit) without any restriction. This brings the power and speed of the C++ native code to the user friendly
 Java based mobile applications. Susa is also a simulation framework for the researchers and engineers who design
 computational systems. It has linear algebra, signal processing and common communications blocks.

The *matrix* and *array* template classes i.e. types are the heart of Susa. A *vector* is a single column (or a single row) matrix. They are bundled with a constellation of classes and functions to process their underlying data.

### Highlights
 - Algebraic types such as *matrix* and multi-dimensional *array* (template classes),
 - High and low precision Fixed-Point types (template classes),
 - Linear algebraic operations and analysis (e.g. Solvers, Determinant and SVD),
 - Signal processing operations (e.g. FFT, Filter (FIR/IIR), Convolution and Random Number Generators),
 - Convolutional Forward Error Correction (FEC) blocks: encoder, MLSE (Viterbi) and MAP (BCJR) decoders,
 - Channel equalisers: MLSE (Viterbi) and MAP (BCJR),
 - Automatic memory management i.e. allocation, deallocation, move and copy.

## Build, Test and Install
### Build
To build Susa you need to have a C++17 compiler, Make and [CMake](https://cmake.org) installed.

```
mkdir build
cmake -S . -B build/Debug -D CMAKE_BUILD_TYPE=Debug
cmake --build build/Debug -j
```

### Test
It is highly recommended to run the tests after the build.

```
make test
```
Should you verify which test(s) has/have been failed, run the following for a more detailed report.

```
ctest -V
```
### Install
Once it has been built and tested you are ready to code. Assuming your current path is `build` directory, run
```
make install
```
to be able to build against Susa system-wide. However, you may continue using the local build without installation.
## Examples
In the [examples](https://github.com/libsusa/susa/tree/master/examples) directory
a number of simulation and tutorial source codes have been provided.
## Contribution
This is a non-profit project and it belongs to its users. You can contribute to your project by reporting bugs and extending it by following the provided [guidelines](https://guides.github.com/activities/forking). This paves the way for further improvements and protects the authors' rights.
## History
Susa was born in April 2008 out of a university project course in digital communications.
At the time the libraries that could be used for digital communications simulation had
many dependencies (e.g. LAPACK, BLAS and ATLAS).
Once it took about six hours on a decent PC to compile one of them. On the other hand those
weighty codes had nested bugs that sometimes stemmed from their third party dependencies.
The answer to these problems was Susa that was [released](http://sourceforge.net/projects/susa)
in November 2008.

Later in early 2009, Susa was used for a bandwidth efficient coding scheme design and simulation, namely,
[Faster Than Nyquist (FTN)](http://www.eit.lth.se/fileadmin/eit/courses/eit085f/AndersonFasterThanNyquistSignaling.pdf).
It required performant equalizers to decode up to some twenty taps (compared to the fading channels with few taps).
The simulation of such systems took a long time between an hour to a few days. This library could simulate
an FTN system with thirteen taps using a modified BCJR algorithm (a sub-optimal variant that could outperform
the original algorithm) in about an hour whereas a similar script in a commercial computing software took
at least twelve hours. Non-orthogonal a.k.a. FTN waveforms have been considered for the next generation of radio access networks (6th generation).

[Unearthed tablets from Susa (2000 BC)](https://www.britannica.com/science/pi-mathematics) revealed a rather precise calculation of $\pi = 3.125$ with the fractional part
whereas the other earlier efforts calculated only the integer part.
Since the very first line of code was simply the definition of constant Pi, it has been named Susa.

## License
[Susa](http://libsusa.org) has been released under GNU Lesser General Public License (LGPL).


[codecov-link]:  https://codecov.io/gh/libsusa/susa
[codecov-badge]: https://codecov.io/gh/libsusa/susa/branch/master/graph/badge.svg
[travis-link]: https://travis-ci.com/libsusa/susa
[travis-badge]: https://travis-ci.com/libsusa/susa.svg?branch=master
