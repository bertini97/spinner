/* lattice.c
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

#include "spinner.h"
#include "utils.h"
#include "error.h"

spnr_system_t *
spnr_system_alloc (spnr_kind_t const * kind,
                   spnr_graph_t * graph)
{
  spnr_system_t * s = malloc_err (sizeof (spnr_system_t));
  
  if (!kind || !graph)
    spnr_err (SPNR_ERROR_ARG_NULL, "requested argument is NULL");
  s->kind = kind;
  s->graph = graph;
  s->priv = kind->priv_alloc (graph);
  
  return s;
}

void
spnr_system_free (spnr_system_t * const s)
{
  s->kind->priv_free (s->priv);
  free (s);
}

void
spnr_system_spins_set_up (spnr_system_t * const s)
{
  s->kind->spins_set_up (s->graph, s->priv);
}

void
spnr_system_spins_set_rand (spnr_system_t * const s)
{
  s->kind->spins_set_rand (s->graph, s->priv);
}

void
spnr_system_spins_print_2d (spnr_system_t const * const s)
{
  s->kind->spins_print_2d (s->graph, s->priv);
}

void
spnr_system_spins_print_3d (spnr_system_t const * const s)
{
  s->kind->spins_print_3d (s->graph, s->priv);
}

float
spnr_system_calc_h (spnr_system_t const * const s)
{
  float h;
  h = s->kind->calc_h (s->graph, s->priv);
  return h; 
}

float
spnr_system_calc_m (spnr_system_t const * const s)
{
  float h;
  h = s->kind->calc_m (s->graph, s->priv);
  return h; 
}

spnr_func_t *
spnr_metr (spnr_kind_t const * kind)
{
  if (!kind->mcstep_metr)
    spnr_err (SPNR_ERROR_FUNC_NULL,
              "ERROR: func not found for this lattice kind");
  return kind->mcstep_metr;
}

spnr_func_t *
spnr_wolff (spnr_kind_t const * kind)
{
  if (!kind->mcstep_wolff)
    spnr_err (SPNR_ERROR_FUNC_NULL,
              "ERROR: func not found for this lattice kind");
  return kind->mcstep_wolff;
}

spnr_data_t *
spnr_system_run (spnr_func_getter_t *getter, spnr_system_t * s,
               size_t n_steps, size_t n_probes, float temp)
{
  size_t i, j;
  size_t const n_steps_bef_probe = n_steps / n_probes;
  float const beta = 1.0 / temp;
  
  spnr_graph_t *graph = s->graph;
  void * priv = s->priv;
  float (*calc_h) (spnr_graph_t *graph, void *priv) = s->kind->calc_h;
  float (*calc_m) (spnr_graph_t *graph, void *priv) = s->kind->calc_m;
  void (*func) (spnr_graph_t *graph, void *priv, float beta) = getter(s->kind);
  
  spnr_data_t * const data = spnr_data_alloc (n_probes + 1);
  
  srand (time (0));
  
  data->h[0] = calc_h (graph, priv);
  data->m[0] = calc_m (graph, priv);
  
  for (i = 1; i <= n_probes; ++i)
    {
      for (j = 0; j < n_steps_bef_probe; ++j)
          func (s->graph, priv, beta);
      data->h[i] = calc_h (graph, priv);
      data->m[i] = calc_m (graph, priv);
    }
  
  return data;
}
