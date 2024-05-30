/* cubic.c
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

#define SPNR_DIMS_MAX 8

typedef struct
{
  size_t L;
  size_t D;
  
  float *J;
  size_t *neighbors;
} cubic_priv_t;

static void *
priv_alloc (float (*getter)(), size_t const N, size_t const D)
{
  size_t i, j;
  size_t unit, row, last_row, left;
  size_t slices[SPNR_DIMS_MAX];
  
  /*if (D <= 0 || D >= SPNR_DIMS_MAX)
    spnr_err (SPNR_ERROR_PARAM_OOB, "graph parameters out of bounds");*/
  
  cubic_priv_t *priv = malloc_err (sizeof(cubic_priv_t));
  
  priv->L = nearbyintf (pow (N, 1.0/D));
  priv->D = D;
  priv->J = malloc_err (N*2*D * sizeof (float));
  priv->neighbors = malloc (N*2*D * sizeof (size_t));

  for (i = 0; i <= D; ++i)
    slices[i] = pow (priv->L, i);
  
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < D; ++j)
        {
          unit = slices[j];
          row  = slices[j+1];
          last_row = row - unit;
          left = ((i % row) < unit) ? (i + last_row) : (i - unit);
      
          priv->J[i*2*D + j] = getter ();
          priv->J[left*2*D + j + D] = priv->J[i*2*D + j];
          
          priv->neighbors[i*2*D + j] = left;
          priv->neighbors[left*2*D + j + D] = i;
        }
    }
  
  return priv;
}

static void
priv_free (void * const priv)
{
  cubic_priv_t *priv_cast = (cubic_priv_t*)priv;
  free (priv_cast->J);
  free (priv_cast->neighbors);
  free (priv_cast);
}

static float
calc_delta_h (void const * const priv,
              spnr_sys_t const * const sys,
              void const * const prop,
              size_t const k)
{
  cubic_priv_t const * const priv_ = (cubic_priv_t *) priv;
  size_t const stride = 2 * priv_->D;
  
  return sys->kind->calc_delta_h_binary (sys->priv, stride,
                                         priv_->J + k * stride,
                                         priv_->neighbors + k * stride,
                                         prop, k);
}

static float
calc_h (void const * const priv, size_t const N, spnr_sys_t const * const sys)
{
  cubic_priv_t const * const priv_ = (cubic_priv_t *) priv;
  size_t i, stride = 2 * priv_->D;
  float h = 0;
  
  for (i = 0; i < N; ++i)
    h += sys->kind->calc_part_h_binary (sys->priv, stride,
                                        priv_->J + i * stride,
                                        priv_->neighbors + i * stride,
                                        i);
  
  return h / (2.0*N);
}

static const spnr_graph_kind_t cubic_kind =
{
  "cubic",
  &priv_alloc,
  &priv_free,
  &calc_delta_h,
  &calc_h
};

const spnr_graph_kind_t *spnr_cubic = &cubic_kind;