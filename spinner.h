/* spinner.h
 * 
 * Copyright (C) 2024 L. Bertini
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

#include <stdlib.h>
#include <stdio.h>

#undef BEGIN_C_DECLS
#undef END_C_DECLS
#ifdef __cplusplus
# define BEGIN_C_DECLS extern "C" {
# define END_C_DECLS }
#else
# define BEGIN_C_DECLS /* empty */
# define END_C_DECLS /* empty */
#endif

#define SPNR_FALSE 0
#define SPNR_TRUE 1

BEGIN_C_DECLS

/* struct forward declaration */

typedef struct spnr_graph_struct spnr_graph_t;
typedef struct spnr_sys_struct spnr_sys_t;
typedef struct spnr_step_struct spnr_step_t;

/* System object
 *
 * Opaque object representing a spin system
 * TODO: comment architecture
 */

typedef struct
{
  char const * name;
  void * (*priv_alloc) (size_t N, size_t param);
  void (*priv_free) (void *priv);
  size_t (*spin_size) (void *priv);
  void (*set_up) (void *priv, size_t N);
  void (*set_rand) (void *priv, size_t N);
  
  void (*fill_prop) (void const *priv, void *prop, size_t k);
  void (*accept_prop) (void *priv, void const *prop, size_t k);
  
  float (*calc_delta_h_binary) (void const *priv, size_t n_sites,
                                float const *J, size_t const *sites,
                                void const *prop, size_t k);
  float (*calc_part_h_binary) (void const* priv, size_t n_sites,
                               float const *J, size_t const *sites,
                               size_t k);
  
  float (*calc_phi) (void const *priv, size_t N);
} spnr_sys_kind_t;

struct spnr_sys_struct
{
  spnr_sys_kind_t const * kind;
  void *priv;
  size_t N;
  spnr_graph_t *graph;
};

/* Available system kinds */

extern spnr_sys_kind_t const *spnr_ising;
extern spnr_sys_kind_t const *spnr_nvector;

/* System object methods */

spnr_sys_t * spnr_sys_alloc (spnr_graph_t *graph, spnr_sys_kind_t const * kind, size_t param);
void spnr_sys_free (spnr_sys_t * sys);
float spnr_sys_spin_size (spnr_sys_t * sys);
float spnr_sys_calc_h (spnr_sys_t * sys);
float spnr_sys_calc_phi (spnr_sys_t * sys);

/* Graph object
 *
 * Opaque object representing a graph
 * TODO: comment architecture
 */

typedef struct
{
  char const * name;
  void * (*priv_alloc) (float (*getter)(), size_t N, size_t param);
  void (*priv_free) (void *priv);
  float (*calc_delta_h) (void const *priv, spnr_sys_t const *sys, void const *prop, size_t k);
  float (*calc_h) (void const *priv, size_t N, spnr_sys_t const *sys);
} spnr_graph_kind_t;

struct spnr_graph_struct
{
  spnr_graph_kind_t const * kind;
  void * priv;
  size_t N;
};

/* Available graph kinds */

extern spnr_graph_kind_t const *spnr_cubic;

/* System object methods */

spnr_graph_t * spnr_graph_alloc (spnr_graph_kind_t const *kind,
                                 float (*getter)(),
                                 size_t N, size_t param);
void spnr_graph_free (spnr_graph_t *graph);

/* Stepper object
 *
 * Opaque object representing a stepper
 * TODO: comment architecture
 */

typedef struct {
  char const * name;
  void * (*priv_alloc) (size_t param);
  void (*priv_free) (void *priv);
  void (*apply) (void *priv, spnr_sys_t const *sys, float beta);
} spnr_step_kind_t;

struct spnr_step_struct {
  spnr_step_kind_t const * kind;
  void *priv;
};

/* Available stepper kinds */

extern spnr_step_kind_t const *spnr_metropolis;

/* System object methods */

spnr_step_t * spnr_step_alloc (spnr_step_kind_t const *kind, size_t param);
void spnr_step_free (spnr_step_t *step);
void spnr_step_apply (spnr_step_t const *step, spnr_sys_t const *sys,
                      float beta);

/* Graph couplings getter functions */

float spnr_ferr ();


/* Data struct
 *
 * Convenience struct to hold simulation data
 * TODO: comment architecture
 */

typedef struct
{
  size_t size;
  float * h;
  float * phi;
} spnr_data_t;

/* Data object methods */

spnr_data_t * spnr_data_alloc (size_t const size);
void spnr_data_free (spnr_data_t *data);
void spnr_data_write (spnr_data_t const *data, char const *fname);
void spnr_data_mean_calc (spnr_data_t const *data,
                          float *h_mean, float * const phi_mean);
void spnr_data_var_calc (spnr_data_t const * const data,
                         float * const h_var, float * const phi_var);
void spnr_data_corr_calc (spnr_data_t *corr, spnr_data_t const *data);
void spnr_data_run_and_probe (spnr_data_t *data, spnr_sys_t *sys,
                              spnr_step_t const *step,
                              float temp, size_t n_steps_before_probe);

END_C_DECLS

#endif