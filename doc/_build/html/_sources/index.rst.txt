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

Lattice parameters and couplings
==========================================

An `n`-vector model is a simple system of interacting n-vectors
:math:`\bold{s}_i`, called spins, on a structure, usually an ordered
lattice, of size :math:`N`. They interact with an Hamiltonian in the
form:

.. math::
   H(\{\bold{s}_i\}) = -\sum_{i,j}^N J_{ij} \bold{s}_i \cdot \bold{s}_j

The coupling constant :math:`J_{ij}` establishes which spins will
interact.

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

Creating lattice parameter data
==========================================

All the information about the structure of a particular lattice, and the
type of interaction between the spins, is contained in the datatype
:c:type:`spnr_params_t`.

.. c:type:: spnr_params_t

   This datatype holds information about a lattice, and can be shared
   between multiple instances of a lattice.

It can be created and destroyed with

.. c:function:: spnr_params_t *spnr_params_cubicnn_alloc (float (*getter)(), size_t side, size_t n_dims, size_t param)
   
   :param getter: pointer to function that returns the desired interaction; available getters are:
   
      .. c:function:: float spnr_ferr ()
      
         ferromagnetic interaction, returns :code:`+SPNR_J`
         
      .. c:function:: float spnr_antiferr ()
      
         antiferromagnetic interaction, returns :code:`-SPNR_J`
      
      .. c:function:: float spnr_bim ()
      
         bimodal interaction, returns :code:`+SPNR_J` or :code:`-SPNR_J` with equal probability
   
   :param side: the side length in spins of the lattice structure
   
   :param n_dims: the number of dimensions of the lattice structure
   
   :param param: a generic parameter that lattice kinds interpret differently: in the n-vector model it represent n, the number of components of the vectors, while in the Potts model it's q, the number of colors
      
   :returns: a pointer to a newly allocated :c:type:`spnr_param_t` object with cubic structure and nearest neighbor interaction

.. c:function:: spnr_params_t *spnr_params_longrange_alloc (float (*getter)(), size_t side, size_t n_dims, size_t param)
   
   :returns: a pointer to a newly allocated :c:type:`spnr_param_t` object with long range interactions; only :code:`size = pow (side, n_dims)` is needed for the computations, and these parameters are otherwise ignored

.. c:function:: void spnr_params_nn_free (spnr_params_t *p)

   Frees the :c:type:`spnr_params_t` onject with nearest neighbor interaction

.. c:function:: void spnr_params_lr_free (spnr_params_t *p)

   Frees the :c:type:`spnr_params_t` onject with long range interaction

Example:

.. code-block:: c
  
   #include "spinner.h"
   
   int main ()
   {
     spnr_params_t *params = spnr_params_cubicnn_alloc (spnr_ferr, 16, 2, 0);
     spnr_prams_nn_free (params);
   }

Creating lattices from the parameter type
==========================================

In this section are covered the basics of making lattice objects from
a single shared set of parameters.

All the information about the kinds of spins that populate the lattice
and the functions needed to operate on them, are contained in the
datatype :c:type:`spnr_latt_kind_t`.

.. c:type:: spnr_kind_t

   Defines a static structure which holds functions that act on a
   particular kind of lattice. It does `not` hold any data, and it's
   shared between lattices.

   The following lattice types are available.

   .. c:var:: spnr_kind_t *spnr_ising_cubicnn

      Ising model on a cubic lattice with nearest neighbor interaction
   
   .. c:var:: spnr_kind_t *spnr_ising_longrange

      Ising model on a graph with long range interaction

   .. c:var:: spnr_kind_t *spnr_nvector_cubicnn_ferr

      `n`-vector model on a cubic lattice with nearest neighbor
      interaction

A lattice is defined using the :c:type:`spnr_latt_t` datatype

.. c:type:: spnr_latt_t

   This datatype represents a lattice object with a defined kind of
   spins and it can share its :c:type:`spnr_param_t` structure with
   other lattice instances

These are the functions for creating and destroying a lattice.

.. c:function:: spnr_latt_t *spnr_latt_alloc (spnr_latt_kind_t *kind, size_t side, size_t n_dims, size_t param)

   :param kind: a pointer to the requested lattice kind

   :param params: a pointer to the requested lattice params

   :returns: a pointer to an allocated :c:type:`spnr_latt_t` variable

.. c:function:: void spnr_latt_free (spnr_latt_t *l)

   Frees every pointer allocated by :code:`spnr_latt_alloc`.

Example:

.. code-block:: c
  
   #include "spinner.h"
   
   int main ()
   {
     spnr_params_t *p = spnr_params_cubicnn_alloc (spnr_ferr, 8, 2, 0);
     spnr_latt_t *l = spnr_latt_alloc (spnr_ising_cubicnn, p);
     spnr_latt_free (l);
     spnr_params_free (p);
   }

Simulating a lattice
==========================================

In this section is covered the process of simulating a lattice with
MCMC methods.

This library provides the following datatype to store the simulation
data.

.. c:type:: spnr_data_t

   Contains a size parameter and two arrays, which hold the energy and
   the magnetization for each time step.

These functions can be used to create or destroy one.

.. c:function:: spnr_data_t *spnr_data_alloc (size_t size)

.. c:function:: void spnr_data_free (spnr_data_t *d)

.. c:function:: void spnr_data_write (spnr_data_t *d, char *fname)

   Writes the simulation data in a plain text :file:`.dat` file with
   the provided name.

Running the simulation requires calling the appropriate function.

.. c:type:: spnr_func_getter_t

   Datatype that holds a pointer to a function that grabs the required
   algorithm from the 

.. c:function:: spnr_data_t *spnr_latt_run (spnr_func_getter_t *getter, spnr_latt_t *l, size_t n_steps, size_t n_probes, float temp)

   Runs the simulation with a single spin flip metropolis algorhithm.

   :param getter: a pointer to the requested algorithm getter, to be chosen from:
   
      .. c:var:: spnr_func_getter_t spnr_metr
      
         gets the single spin-flip Metropolis dynamics
         
      .. c:var:: spnr_func_getter_t spnr_wolff
      
         gets the clister-flip Wolff dynamics
   
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
     spnr_params_t *p = spnr_params_cubicnn_alloc (spnr_ferr, 8, 2, 0);
     spnr_latt_t *l = spnr_latt_alloc (spnr_ising_cubicnn, params);
     
     run = spnr_latt_run (spnr_metr, l, 1000, 100, 3.0);
     corr = spnr_data_calc_corr (run);
     
     spnr_data_mean_calc (run, &h, &m);
     printf("%f %f\n", h, m);
     
     spnr_data_write (run, "data");
     spnr_data_write (corr, "corr");
     
     spnr_data_free (corr);
     spnr_data_free (run);
     spnr_latt_free (l);
     spnr_params_nn_free (p);
   }
