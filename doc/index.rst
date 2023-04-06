.. Spinner documentation master file, created by
   sphinx-quickstart on Tue Apr  4 20:10:48 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

******************************************
Spinner: a spin lattice simulation library
******************************************

This library allows the creation of D-dimensional lattices or more
general structures of spins from `n`-vector models (such as Ising, XY,
Heisenberg, etc), sometimes also called :math:`O(n)` models, and their
simulation using Monte Carlo Markov Chain (MCMC) methods.

.. warning::
   This library is in early stage of development: it might produce
   incorrect results. It comes with absolutely no warranty, and you
   should use it at your own risk. Contributions are welcome.

.. warning::
   So far only the single spin-flip metropolis dynamics has been
   implemented, and the only interaction available is ferromagnetic.

.. note::
   Support for the Potts models and more interactions (long range,
   spin-glass, etc) is planned.

Introduction
==========================================

An `n`-vector model is a simple system of interacting n-vectors
:math:`\bold{s}_i`, called spins, on a structure, usually an ordered
lattice, of size :math:`N`. They interact with an Hamiltonian in the
form:

.. math::
   H(\{\bold{s}_i\}) = -\sum_{i,j}^N J_{ij} \bold{s}_i \cdot \bold{s}_j

The coupling constant :math:`J_{ij}` establishes which spins will
interact.

Non-zero elements of the coupling
------------------------------------------

Which elements of this matrix are nonzero define the extent of the
interaction between the spins, the most common being the
nearest-neighbor (NN) interaction:

.. list-table::
   :header-rows: 1
   
   * - :math:`J_{ij}`
     - Interaction
   * - :math:`\begin{cases}=J,\;|i-j|=1\\=0,\;|i-j|\neq 1\end{cases}`
     - Nearest neighbor interaction
   * - :math:`=J/N\quad \forall i,j`
     - Long range interaction

Magnitude of the couplings
------------------------------------------

Different values for the nonzero coupling elements result in different
models:

.. list-table::
   :header-rows: 1

   * - :math:`J`
     - Interaction
   * - :math:`>0`
     - Ferromagnetic, the spins will try to align to each other.
   * - :math:`<0`
     - Antiferromagnetic, the spins will try to disalign from each other
   * - Random
     - Spin-glass behavior

Creating a lattice
==========================================

In this section are covered the basics of making a lattice object.

The lattice kind
------------------------------------------

All the information about the structure of a particular lattice, the
kinds of spins that populate it and the interation they have with each
other is contained in the datatype :c:type:`spnr_latt_kind_t`.

.. c:type:: spnr_latt_kind_t

   Defines a static structure which holds functions that act on a
   particular lattice. It does `not` hold any data.

   The following lattice types are available.

   .. c:var:: spnr_latt_kind_t *spnr_ising_cubicnn_ferr

      Ising model on a cubic lattice with nearest neighbor ferromagnetic
      interaction.

   .. c:var:: spnr_latt_kind_t *spnr_nvector_cubicnn_ferr

      `n`-vector model on a cubic lattice with nearest neighbor
      ferromagnetic interaction.

The lattice object
------------------------------------------

A lattice is defined using the :c:type:`spnr_latt_t` datatype.

.. c:type:: spnr_latt_t

   This data structure defines a lattice.

   .. c:member:: spnr_latt_kind_t *kind

      Member holding a reference to a datatype
      :c:type:`spnr_latt_kind_t`, which holds all the necessary
      information for that kind of lattice, including functions to
      operate on it.

   .. c:member:: void *priv

      The actual data belonging to the lattice; for example, the array
      of the spin values, among other convenience data. Each
      :c:type:`spnr_latt_kind_t` defines it own internal structure, and
      casts :code:`priv` accordingly.

Creating and destroying a lattice
------------------------------------------

These are the functions for creating and destroying a lattice.

.. c:function:: spnr_latt_t *spnr_latt_alloc (spnr_latt_kind_t *kind, size_t side, size_t n_dims, size_t param)

   :param kind: a pointer to the required lattice kind

   :param side: the size in spins of a side of the lattice

   :param n_dims: the number of dimensions of the lattice

   :param param: an additional parameter which contains the number of vector components in the case of `n`-vector models, or the number of colors in the Potts model. It is ignored in Ising model.

   :returns: a pointer to an allocated :c:type:`spnr_latt_t` variable. It uses the provided :c:type:`spnr_latt_kind_t` to properly allocate the private data.

.. c:function:: void spnr_latt_free (spnr_latt_t *l)

   Frees every pointer allocated by :code:`spnr_latt_alloc`.

Example:

.. code-block:: c
  
   #include "spinner.h"
   
   int main ()
   {
     spnr_latt_t *l = spnr_latt_alloc (spnr_ising_cubicnn_ferr, 3, 3, 0);
     spnr_latt_free (l);
   }

Simulating a lattice
==========================================

In this section is covered the process of simulating a lattice with
MCMC methods.

Simulation data
------------------------------------------

This library provides the following datatype to store the simulation data.

.. c:type:: spnr_data_t

   Contains a size parameter and two arrays, which hold the energy and
   the magnetization for each time step.

These functions can be used to create or destroy one.

.. c:function:: spnr_data_t *spnr_data_alloc (size_t size)

.. c:function:: void spnr_data_free (spnr_data_t *d)

.. c:function:: void spnr_data_write (spnr_data_t *d, char *fname)

   Writes the simulation data in a plain text :file:`.dat` file with
   the provided name.

Running the simulation
------------------------------------------

Running the simulation requires calling the appropriate function.

.. warning:: So far only the single spin flip Metropolis dynamic is
   available. Single spin flip Heat-Bath and cluster flip Wolff are
   being worked on and will be implemented soon

.. c:function:: spnr_data_t *spnr_latt_run_ssf_met (spnr_latt_t *l, size_t n_steps, size_t n_probes, float temp)

   Runs the simulation with a single spin flip metropolis algorhithm.

   :param l: the lattice to sample
   :param n_steps: the number of Monte Carlo steps to run (each MC step is equivalent to N single spin flip steps where N is the size of the lattice)
   :param n_probes: the number of times the lattice is probed for energy and magnetization
   :param temp: the temperature at which the simulation is run
   
   :returns: a pointer to a freshly allocated :c:type:`spnr_data_t` variable which contains the simulation data for each probe in order. This pointer must be freed with :code:`spnr_data_free`.

Example:

.. code-block:: c
  
   #include "spinner.h"
   
   int main ()
   {
     spnr_data_t *run;
     spnr_latt_t *l = spnr_latt_alloc (spnr_ising_cubicnn_ferr, 3, 3, 0);
     
     run = spnr_latt_run_ssf_met (l, 1000, 100, 3.0);
     spnr_data_write (run, "data")
     
     spnr_data_free (run);
     spnr_latt_free (l);
   }

Computing means and correlations
------------------------------------------

This library provides the following functions to study the obtained
data.

.. c:function:: void spnr_data_calc_mean (spnr_data_t *d, float *h, float *m)

   Returns in the pointers provided the mean values for energy and magnetization.

.. c:function:: spnr_data_t *spnr_data_calc_corr (spnr_data_t *d)

   :returns: a pointer to an allocated :c:type:`spnr_data_t` variable holding the temporal correlation for both energy and magnetization. The pointer must be freed with :code:`spnr_data_free`

Example:

.. code-block:: c
  
   #include "spinner.h"
   
   int main ()
   {
     float h = 0, m = 0;
     spnr_data_t *run, *corr;
     spnr_latt_t *l = spnr_latt_alloc (spnr_ising_cubicnn_ferr, 3, 3, 0);
     
     run = spnr_latt_run_ssf_met (l, 1000, 100, 3.0);
     corr = spnr_data_calc_corr (run);
     
     spnr_data_mean_calc (run, &h, &m);
     printf("%f %f\n", h, m);
     
     spnr_data_write (run, "data");
     spnr_data_write (corr, "corr");
     
     spnr_data_free (corr);
     spnr_data_free (run);
     spnr_latt_free (l);
   }
