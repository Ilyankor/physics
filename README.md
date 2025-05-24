# Physics

Reading **Boudreau, Swanson: Applied Computational Physics**.

Libraries used: `Eigen`, header `util/util.h`.
Libraries are placed in `lib`.

## Chapter 03

### Coupled Oscillators

Consider $n$ masses on a wire connected by springs that obey Hooke's law with no forcing or damping.
`src/ch03/oscillation.cpp` computes displacement from equilibrium over time.
Sample input `examples/ch03/oscillation/input.txt` from Exercise 3.13.
Animation `animate.gp` using `gnuplot` and `ffmpeg`.

### Scattering from potential barriers

Scattering of a particle, incident from the left, hitting piecewise constant potential barriers that
are $0$ outside a finite interval, obeying time-independent Schrödinger equation.
`src/ch03/scattering.cpp` computes wavefunction, probabilities, reflection, and transmission
coefficients.
`funcToStep` uses flag `-g` with a user defined namespace to create a piecewise constant version of
a user defined function, outputted as a **partially valid** input file format.
Flag `-3.9` performs Exercise 3.9.

## Chapter 04
