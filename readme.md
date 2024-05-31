![spinner](spinner.png)

# Spinner: a spin system simulation library

This is Spinner, a collection of numerical routines for spin systems simulations. Spinner is free software, you can redistribute it and/or modify it under the terms of the GNU General Public License.

> **WARNING**: This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

> **NOTE**: This library is in early stage of development. Contributions are welcome.

## Introduction

Spin systems, such as the **Ising** or the **Heisenberg** model, are defined on a graph where each vertex $i$ has a (symmetric) interaction $J_{ij}$ with its neighbors $j$. This library allows the creation of a *graph object* that is used by one or more *system objects*. You can then apply a *stepper object* to move the system in its configuration space using a **Monte Carlo Markov Chain** (MCMC) method.

The system object can hold spin variables of different kinds, such as Ising variables ($s_i=\pm1$), or n-vector variables ($\mathbf{s}_i\in\mathbb{R}^n:|s_i|=1$). The graph kind holds the coupling constants $J_{ij}$ that make up the graph. Graph structures can be of the following kinds:

![graph_kinds](graph_kinds.svg)

Read the [documentation](https://bertini97.gitlab.io/spinner/) for a detailed explanation.

## Building and using

Here are the building instructions:
```bash
git clone <spinner repo>
cd spinner
autoreconf -i
./configure
make
sudo make install
```

The default installation prefix is `/usr/local/lib`.

## Usage example

Here is an example code that, at different temperatures, samples:
* an Ising system
* on a cubic lattice
* with ferromagnetic interactions
```c
#include <stdio.h>
#include <spinner/spinner.h>

int main()
{
  size_t N = 10000, n_temps = 30, n_probes = 1000+1, i;
  float temp, temp_start = 0.1, temp_finish=3.0, temp_step, beta;
  float h_mean, m_mean, h_var, m_var;
  FILE * f;
  
  spnr_graph_t *graph = spnr_graph_alloc (spnr_cubic, spnr_ferr, N, 2);
  spnr_sys_t *sys = spnr_sys_alloc (graph, spnr_ising, 2);
  spnr_step_t *step = spnr_step_alloc (spnr_metropolis,
                                       spnr_sys_spin_size(sys));
  spnr_data_t *data = spnr_data_alloc (n_probes);
  
  temp_step = (temp_finish - temp_start) / (float) (n_temps-1);
  f = fopen ("results.txt", "w");
  for (i = 0; i < n_temps; ++i)
    {
      temp = temp_start + i * temp_step;
      beta = 1.0 / temp;
      
      printf ("temp=%f\n", temp);
      spnr_data_run_and_probe (data, sys, step, temp, 10);
      spnr_data_mean_calc (data, &h_mean, &m_mean);
      spnr_data_var_calc (data, &h_var, &m_var);
      fprintf (f, "%f %+f %+f %+f %+f\n", temp, h_mean, m_mean,
               beta * beta * N * h_var, beta * N * m_var);
    }
  fclose (f);
  
  spnr_data_free (data);
  spnr_step_free (step);
  spnr_sys_free (sys);
  spnr_graph_free (graph);
}
```
To compile:
```bash
gcc -I /usr/local/include/ program.c -o program -L /usr/local/lib -lspinner -lm
```

## Planned features:
Here are the planned features for this library. For now, only some of them are in place.
 - Variables
	 - [x] Ising
	 - [x] n-vector (XY, Heisenberg)
	 - [ ] Potts
 - Graphs
	 - [x] Cubic lattice
	 - [ ] Fully connected
	 - [ ] Random graph
 - Interactions
	 - [x] Ferromagnetic
	 - [x] Antiferromagnetic
	 - [ ] Spin glass (bimodal)
 - Steppers
	 - [x] Metropolis
	 - [ ] Heat-Bath
	 - [ ] Wolff
	 - [ ] Swendsen–Wang
 - Utilities
	 - [x] Simulation data object
	 - [ ] Parallelization
	 - [ ] Parallel tempering

