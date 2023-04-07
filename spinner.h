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

typedef void spnr_func_t (void *priv, float beta);

/* Opaque structure representing a lattice kind, containing the functions
 * necessary to create it, simulate it with MCMC methods and
 * measure the relevant quantities
 */
typedef struct
{
  char const * name;
  void *(*alloc) (size_t const side, size_t n_dims, size_t const param);
  void (*free) (void *priv);
  void (*print_spins_2d) (void *priv);
  void (*print_spins_3d) (void *priv);
  float (*calc_h) (void *priv);
  float (*calc_m) (void *priv);
  void  (*mcstep_metr) (void *priv, float beta);
  void  (*mcstep_wolff) (void *priv, float beta);
} spnr_latt_kind_t;

/* Structure representing a lattice containing a lattice type structure
 * and a pointer to the private data, that has a different structure
 * for every lattice kind
 */
typedef struct
{
  spnr_latt_kind_t const * kind;
  void * priv;
} spnr_latt_t;


/* Structure for a double vector of data, holding the energy and
 * magnetization for every time step
 */
typedef struct
{
  size_t size;
  float *h;
  float *m;
} spnr_data_t;


/* Available lattice kinds */

extern spnr_latt_kind_t const *spnr_ising_cubicnn_ferr;
extern spnr_latt_kind_t const *spnr_nvector_cubicnn_ferr;

/* Available func methods */
extern spnr_func_t *spnr_metr (spnr_latt_t * const l);
extern spnr_func_t *spnr_wolff (spnr_latt_t * const l);

/*lattice structure methods */

spnr_latt_t *
spnr_latt_alloc (spnr_latt_kind_t const * const kind,
                 size_t const side,
                 size_t const n_dims,
                 size_t const param);

void
spnr_latt_free (spnr_latt_t * const l);

float
spnr_latt_calc_h (spnr_latt_t const * const l);

float
spnr_latt_calc_m (spnr_latt_t const * const l);

spnr_data_t *
spnr_latt_run (spnr_func_t * (*getter) (spnr_latt_t *l),
               spnr_latt_t * const l, size_t const n_steps,
               size_t const n_probes, float const beta);

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
