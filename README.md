# Susa Open Source Project

Susa is a mathematics and signal processing C++ framework based on KISS principle. The mobility is the key feature of Susa. It is a stand-alone and modern signal processing framework that can easily be ported to mobile platforms. It is designed not to have any dependencies to any none standard third party library. A C++11 compiler is necessary and sufficient in order to compile it. Therefore, Susa can be exploited in mobile platforms such as Android NDK (Native Development Toolkit) without restriction. This brings the power and speed of the C++ native code to the user friendly Java based mobile applications. Susa is also a simulation framework for the researchers and engineers who design computational systems. It has linear algebra, signal processing and communications common methods.

### Highlights
 - Linear algebraic operations and analysis (e.g. Determinant and SVD).
 - Signal processing operations (e.g. FFT, Filter, Convolution and Random Number Generators).
 - Convolutional Forward Error Correction (FEC) computational blocks (Encoder, Viterbi and BCJR decoders).
 - Channel equalizers (Viterbi and BCJR decoders).

## Build and Test Susa
To build Susa you need to have a C++ compiler, Make and [CMake](https://cmake.org) installed.

```
mkdir build
cd build
cmake ..
make
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
In the [examples](https://github.com/behrooza/susa/tree/master/examples) directory a number of simulation and tutorial source codes have been provided.
## History
Susa was born in April 2008 out of a university project course in digital communications. At the time the libraries
that could be used for digital communications simulation had many dependencies (e.g. LAPACK, BLAS and ATLAS).
Once it took about six hours on a decent PC to compile one of them. On the other hand those weighty codes had nested
bugs that sometimes stemed from their third party dependencies. The answer to these problems was Susa that was
[released](http://sourceforge.net/projects/susa) in November 2008.

Later in early 2009, Susa was used for a bandwidth efficient coding scheme, namely, [Faster Than Nyquist](http://www.eit.lth.se/fileadmin/eit/courses/eit085f/AndersonFasterThanNyquistSignaling.pdf). It required performant equalizers to decode up to some twenty taps (compared to the fading channel with few taps). The simulation of such systems took a long time between an hour to a few days. This library could simulate an FTN system with thirteen taps using a modified BCJR algorithm (a suboptimal variant that could outperform the original algorithm) in about an hour wheras a similar script in a commercial computing software took at least twelve hours.
## License
[Susa](http://susalib.org) has been released under GNU Lesser General Public License (LGPL).
