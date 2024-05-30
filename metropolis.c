/* metropolis.c
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

#include <math.h>

#include "spinner.h"
#include "error.h"

static void *
priv_alloc (size_t const spin_size)
{
  return malloc_err (spin_size);
}

static void
priv_free (void * priv)
{
  free (priv);
}

static int
metr_prop_accept (float const delta_h, float const beta)
{
  if (delta_h <= 0)
    return SPNR_TRUE;
  else
    {
      float acc_ratio = exp (- beta * delta_h);
      float rand_num = (float) rand () / RAND_MAX;
      if (rand_num < acc_ratio)
        return SPNR_TRUE;
      else
        return SPNR_FALSE;
  }
}


/* here priv is not a struct but holds the proposed value */
void
apply (void * const priv, spnr_sys_t const * const sys,
       float const beta)
{
  spnr_graph_t const * const graph = sys->graph;
  size_t i, k;
  size_t const N = graph->N;
  void *prop;
  float delta_h;
  
  for (i = 0; i < N; ++i)
    {
      k = rand() % N;
      sys->kind->fill_prop (sys->priv, priv, k);
      delta_h = graph->kind->calc_delta_h (graph->priv, sys, priv, k);
      
      if (metr_prop_accept (delta_h, beta))
        sys->kind->accept_prop (sys->priv, priv, k);
    }
}

static const spnr_step_kind_t metropolis_kind =
{
  "metropolis",
  &priv_alloc,
  &priv_free,
  &apply
};

const spnr_step_kind_t *spnr_metropolis = &metropolis_kind;