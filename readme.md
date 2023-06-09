![spinner](spinner.png)

# Spinner: a spin lattice simulation library

This library allows the creation of *D*-dimensional lattices or more
general structures of spins from *n*-vector models (such as Ising, XY,
Heisenberg, etc), sometimes also called *O(n)* models, and their
simulation using Monte Carlo Markov Chain (MCMC) methods.

> **WARNING**: This library is in early stage of development: it might
produce incorrect results. It comes with absolutely no warranty, and
you should use it at your own risk. Contributions are welcome.

> **WARNING**: So far only the single spin-flip Metropolis and the cluster
flip Wolff dynamics has been implemented, and the only interaction
available is ferromagnetic.

> **NOTE**: Support for the Potts models and more interactions (long range,
spin-glass, etc) is planned.

Read the [documentation](https://bertini97.gitlab.io/spinner/) for
basic usage examples.
