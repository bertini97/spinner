/* spinner.h
 * 
 * Copyright (C) 2023 L. Bertini
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef SPNR_H
#define SPNR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#undef BEGIN_C_DECLS
#undef END_C_DECLS
#ifdef __cplusplus
# define BEGIN_C_DECLS extern "C" {
# define END_C_DECLS }
#else
# define BEGIN_C_DECLS /* empty */
# define END_C_DECLS /* empty */
#endif

#define SPNR_FALSE          0
#define SPNR_TRUE           1

#define SPNR_PI             3.1415926536
#define SPNR_J              1
#define SPNR_BUFFER_SIZE    100
#define SPNR_SIDE_MAX       512
#define SPNR_DIMS_MAX       8
#define SPNR_COMPS_MAX      8

BEGIN_C_DECLS

/* structure holding a coupling matrix */
typedef struct
{
  size_t side;
  size_t n_dims;
  size_t param;
  size_t size;
  size_t n_inters;
  
  float * coups;
  size_t * nbors;
} spnr_graph_t;

/* Opaque structure representing a lattice kind, containing the functions
 * necessary to create it, simulate it with MCMC methods and
 * measure the relevant quantities
 */
typedef struct
{
  char const * name;
  void * (*priv_alloc) (spnr_graph_t * graph);
  void (*priv_free) (void *priv);
  void (*spins_set_up) (spnr_graph_t *graph, void *priv);
  void (*spins_set_rand) (spnr_graph_t *graph, void *priv);
  void (*spins_print_2d) (spnr_graph_t *graph, void *priv);
  void (*spins_print_3d) (spnr_graph_t *graph, void *priv);
  float (*calc_h) (spnr_graph_t *graph, void *priv);
  float (*calc_m) (spnr_graph_t *graph, void *priv);
  void  (*mcstep_metr) (spnr_graph_t *graph, void *priv, float beta);
  void  (*mcstep_wolff) (spnr_graph_t *graph, void *priv, float beta);
} spnr_kind_t;

/* Structure representing a lattice containing a lattice type structure
 * and a pointer to the private data, that has a different structure
 * for every lattice kind
 */
typedef struct
{
  spnr_kind_t const * kind;
  spnr_graph_t * graph;
  void * priv;
} spnr_system_t;


/* Structure for a double vector of data, holding the energy and
 * magnetization for every time step
 */
typedef struct
{
  size_t size;
  float *h;
  float *m;
} spnr_data_t;

/* lord forgive me */

typedef void spnr_func_t (spnr_graph_t *graph, void *priv, float beta);
typedef spnr_func_t * spnr_func_getter_t (spnr_kind_t const *kind);

/* available lattice kinds */

extern spnr_kind_t const *spnr_ising_cubicnn;
extern spnr_kind_t const *spnr_ising_longrange;

extern spnr_kind_t const *spnr_nvector_cubicnn;

/* available interactions */

extern float spnr_ferr ();
extern float spnr_antiferr ();
extern float spnr_bim ();

/* available func methods */

extern spnr_func_t *spnr_metr (spnr_kind_t const * kind);
extern spnr_func_t *spnr_wolff (spnr_kind_t const * kind);

/* coupling matrix methods */

spnr_graph_t *
spnr_graph_cubicnn_alloc (float (*getter)(), size_t side, size_t n_dims, size_t param);

spnr_graph_t *
spnr_graph_longrange_alloc (float (*getter)(), size_t side, size_t n_dims, size_t param);

void
spnr_graph_nn_free (spnr_graph_t *graph);

void
spnr_graph_lr_free (spnr_graph_t *graph);

/* lattice structure methods */

spnr_system_t *
spnr_system_alloc (spnr_kind_t const * kind, spnr_graph_t *graph);

void
spnr_system_free (spnr_system_t * const s);

float
spnr_system_calc_h (spnr_system_t const * const s);

float
spnr_system_calc_m (spnr_system_t const * const s);

spnr_data_t *
spnr_latt_run (spnr_func_getter_t *getter, spnr_system_t * s,
               size_t n_steps, size_t n_probes, float temp);

void spnr_system_spins_set_up (spnr_system_t * const s);
void spnr_system_spins_set_rand (spnr_system_t * const s);
void spnr_system_spins_print_2d (spnr_system_t const * const s);
void spnr_system_spins_print_3d (spnr_system_t const * const s);

/*data structure methods */

spnr_data_t *
spnr_data_alloc (size_t const size);

void
spnr_data_free (spnr_data_t *d);

void
spnr_data_write (spnr_data_t const * const d,
                 char * const fname);

spnr_data_t *
spnr_data_corr_calc (spnr_data_t const * const data);

void
spnr_data_mean_calc (spnr_data_t const * const d,
                     float *h,
                     float *m);

END_C_DECLS
#endif
