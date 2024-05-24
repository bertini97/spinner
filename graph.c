/* graph.c
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
#include "error.h"
#include "utils.h"

spnr_graph_t *
spnr_graph_cubicnn_alloc (float (*getter)(),
                           size_t side, size_t n_dims, size_t param)
{
  size_t i, j, unit, row, last_row, left, slices[SPNR_DIMS_MAX];
  size_t size = pow (side, n_dims), n_inters = 2 * n_dims;
  spnr_graph_t *g = malloc (sizeof (spnr_graph_t));
  
  if (side <= 0 || side >= SPNR_SIDE_MAX
      || n_dims <= 0 || n_dims >= SPNR_DIMS_MAX)
    spnr_err (SPNR_ERROR_PARAM_OOB, "lattice parameters out of bounds");
  
  g->side = side;
  g->n_dims = n_dims;
  g->param = param;
  g->size = size;
  g->n_inters = n_inters;
  
  g->coups = malloc_err (g->size * g->n_inters * sizeof (float));
  g->nbors = malloc_err (g->size * g->n_inters * sizeof (size_t));

  for (i = 0; i <= n_dims; ++i)
    slices[i] = pow (side, i);
  
  for (i = 0; i < size; ++i)
    {
      for (j = 0; j < n_dims; ++j)
        {
          unit = slices[j];
          row  = slices[j+1];
          last_row = row - unit;
          left = ((i % row) < unit) ? (i + last_row) : (i - unit);
      
          g->coups[i * n_inters + j] = getter ();
          g->coups[left * n_inters + j + n_dims] = g->coups[i * n_inters + j];
          
          g->nbors[i * n_inters + j] = left;
          g->nbors[left * n_inters + j + n_dims] = i;
        }
    }
  
  return g;
}

spnr_graph_t *
spnr_graph_longrange_alloc (float (*getter)(),
                             size_t side, size_t n_dims, size_t param)
{
  size_t i, j, size = pow (side, n_dims);
  spnr_graph_t *g = malloc_err (sizeof (spnr_graph_t));
  
  if (side <= 0 || side >= SPNR_SIDE_MAX
      || n_dims <= 0 || n_dims >= SPNR_DIMS_MAX)
    spnr_err (SPNR_ERROR_PARAM_OOB, "lattice parameters out of bounds");
  
  g->side = side;
  g->n_dims = n_dims;
  g->param = param;
  
  g->size = size;
  g->n_inters = size;
  
  g->coups = malloc_err (g->size * g->n_inters * sizeof (float));
  
  for (i = 0; i < size; ++i)
    {
      for (j = i; j < size; ++j)
        {
          g->coups[i * size + j] = getter ();
          g->coups[j * size + i] = g->coups[i * size + j];
        }
    }
  
  return g;
}

void
spnr_graph_nn_free (spnr_graph_t *g)
{
  free (g->coups);
  free (g->nbors);
  free (g);
}

void
spnr_graph_lr_free (spnr_graph_t *g)
{
  free (g->coups);
  free (g);
}
