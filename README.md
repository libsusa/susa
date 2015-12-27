# Susa Open Source Project

Susa is a mathematics and signal processing C++ framework based on KISS principle. The mobility is the key feature of Susa. It is a stand-alone and modern signal processing framework that can easily be ported to mobile platforms. It is designed not to have any dependencies to any none standard third party library. A C++11 compiler is necessary and sufficient in order to compile it. Therefore, Susa can be exploited in mobile platforms such as Android NDK (Native Development Toolkit) without restriction. This brings the power and speed of the C++ native code to the user friendly Java based mobile applications. Susa is also a simulation framework for the researchers and engineers who design computational systems. It has linear algebra, signal processing and communications common methods.

### Highlights
 - Linear algebraic operations and analysis (e.g. Determinant and SVD).
 - Signal processing operations (e.g. FFT, Filter, Convolution and Random Number Generators).
 - Convolutional Forward Error Correction (FEC) computational blocks (Encoder, Viterbi and BCJR decoders).
 - Channel equalizers (Viterbi and BCJR decoders).

## Build and Test Susa
To build Susa you need to have a C++ compiler, *Make* and [CMake](https://cmake.org) installed.

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
In the *examples* directory a number of simulation and tutorial source codes have been provided.
## License
[Susa](http://susalib.org) has been released under GNU Lesser General Public License (LGPL).
