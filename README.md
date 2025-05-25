# Physics

Reading [**Boudreau, Swanson: Applied Computational Physics**](https://academic.oup.com/book/26392).

## Build

Build with [CMake](https://cmake.org/).

```shell
mkdir build
cd build
cmake ..
cmake --build .
cd ../
```

Binaries are stored in `build/bin`, sorted by chapter.

## Libraries

Libraries are placed in `lib`.
3rd party libraries used:

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)

### Headers

- `futil`: file utilities
- `mutil`: math utilities

## Usage

Sample commands are located in `commands.txt`.
Examples are located in `examples`, sorted by chapter.

## Chapter 03

### Coupled Oscillators

Consider $n$ masses on a wire connected by springs that obey Hooke's law with no forcing or damping.
`src/ch03/oscillation.cpp` computes displacement from equilibrium over time.
Sample input `examples/ch03/oscillation/input.txt` is from Exercise 3.13.
Animation with script `animate.gp` using [gnuplot](https://sourceforge.net/projects/gnuplot/) and [FFmpeg](https://ffmpeg.org/).

### Scattering from potential barriers

Scattering of a particle, incident from the left, hitting piecewise constant potential barriers that
are $0$ outside a finite interval, obeying time-independent Schrödinger equation.
`src/ch03/scattering.cpp` computes wavefunction, probabilities, reflection, and transmission
coefficients.
`funcToStep` uses flag `-g` with a user defined namespace to create a piecewise constant version of
a user defined function, outputted as a **partially valid** input file format.
Flag `-3.9` performs Exercise 3.9.

## Chapter 04

### Interpolation and extrapolation

Function interpolation and extrapolation methods implemented in
`lib/mutil/interp.h`.
