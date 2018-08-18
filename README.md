# Susa Open Source Project  [![Build Status](https://travis-ci.com/susalib/susa.svg?branch=master)](https://travis-ci.com/susalib/susa)

Susa is a mathematics and signal processing C++ framework based on [KISS](https://en.wikipedia.org/wiki/KISS_principle)
 principle. It is stand-alone with a modern architecture. It is designed not to have any dependencies to none standard third
 party libraries. Indeed, a C++11 compiler along with STL is necessary and sufficient in order to compile it. Therefore,
 portability is the key feature of Susa. For example it can be exploited in mobile platforms such as Android NDK (Native
 Development Toolkit) without any restriction. This brings the power and speed of the C++ native code to the user friendly
 Java based mobile applications. Susa is also a simulation framework for the researchers and engineers who design
 computational systems. It has linear algebra, signal processing and common communications blocks.

A matrix and array template class are at the heart of Susa with a constellation of classes and functions
that take input matrices and return results as output matrices.

### Highlights
 - Matrix and multi-dimensional array types (template classes)
 - Linear algebraic operations and analysis (e.g. Determinant and SVD).
 - Signal processing operations (e.g. FFT, Filter (FIR/IIR), Convolution and Random Number Generators).
 - Convolutional Forward Error Correction (FEC) blocks: encoder, MLSE (Viterbi) and MAP (BCJR) decoders.
 - Channel equalisers: MLSE (Viterbi) and MAP (BCJR).

## Build and Test Susa
To build Susa you need to have a C++ compiler, Make and [CMake](https://cmake.org) installed.

```
mkdir build
cd build
cmake ..
make
make install
```
It is highly recommended to run the tests after the build.

```
make test
```
Should you verify which test(s) has/have been failed, run the following for a more detailed report.

```
ctest -V
```
## Examples
In the [examples](https://github.com/behrooza/susa/tree/master/examples) directory
a number of simulation and tutorial source codes have been provided.
## History
Susa was born in April 2008 out of a university project course in digital communications.
At the time the libraries that could be used for digital communications simulation had
many dependencies (e.g. LAPACK, BLAS and ATLAS).
Once it took about six hours on a decent PC to compile one of them. On the other hand those
weighty codes had nested bugs that sometimes stemmed from their third party dependencies.
The answer to these problems was Susa that was [released](http://sourceforge.net/projects/susa)
in November 2008.

Later in early 2009, Susa was used for a bandwidth efficient coding scheme, namely,
[Faster Than Nyquist](http://www.eit.lth.se/fileadmin/eit/courses/eit085f/AndersonFasterThanNyquistSignaling.pdf).
It required performant equalizers to decode up to some twenty taps (compared to the fading channel with few taps).
The simulation of such systems took a long time between an hour to a few days. This library could simulate
a FTN system with thirteen taps using a modified BCJR algorithm (a suboptimal variant that could outperform
the original algorithm) in about an hour whereas a similar script in a commercial computing software took
at least twelve hours.

Unearthed tablets from Susa (2000 BC) revealed a rather precise calculation of Pi = 3.125 with the fractional part
whereas the other earlier efforts calculated the integer part.
Since the very first line of code was simply the definition of constant Pi, it has been named Susa.

## License
[Susa](http://susalib.org) has been released under GNU Lesser General Public License (LGPL).
