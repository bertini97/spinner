/* params.c
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

spnr_params_t *
spnr_params_cubicnn_alloc (float (*getter)(),
                           size_t side, size_t n_dims, size_t param)
{
  size_t i, j, unit, row, last_row, left, slices[SPNR_DIMS_MAX];
  size_t size = pow (side, n_dims), n_inters = 2 * n_dims;
  spnr_params_t *p = malloc (sizeof (spnr_params_t));
  
  if (side <= 0 || side >= SPNR_SIDE_MAX
      || n_dims <= 0 || n_dims >= SPNR_DIMS_MAX)
    spnr_err (SPNR_ERROR_PARAM_OOB, "lattice parameters out of bounds");
  
  p->side = side;
  p->n_dims = n_dims;
  p->param = param;
  p->size = size;
  p->n_inters = n_inters;
  
  p->coups = malloc_err (p->size * p->n_inters * sizeof (float));
  p->nbors = malloc_err (p->size * p->n_inters * sizeof (size_t));

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
      
          p->coups[i * n_inters + j] = getter ();
          p->coups[left * n_inters + j + n_dims] = p->coups[i * n_inters + j];
          
          p->nbors[i * n_inters + j] = left;
          p->nbors[left * n_inters + j + n_dims] = i;
        }
    }
  
  return p;
}

spnr_params_t *
spnr_params_longrange_alloc (float (*getter)(),
                             size_t side, size_t n_dims, size_t param)
{
  size_t i, j, size = pow (side, n_dims);
  spnr_params_t *p = malloc_err (sizeof (spnr_params_t));
  
  if (side <= 0 || side >= SPNR_SIDE_MAX
      || n_dims <= 0 || n_dims >= SPNR_DIMS_MAX)
    spnr_err (SPNR_ERROR_PARAM_OOB, "lattice parameters out of bounds");
  
  p->side = side;
  p->n_dims = n_dims;
  p->param = param;
  
  p->size = size;
  p->n_inters = size;
  
  p->coups = malloc_err (p->size * p->n_inters * sizeof (float));
  
  for (i = 0; i < size; ++i)
    {
      for (j = i; j < size; ++j)
        {
          p->coups[i * size + j] = getter ();
          p->coups[j * size + i] = p->coups[i * size + j];
        }
    }
  
  return p;
}

void
spnr_params_nn_free (spnr_params_t *p)
{
  free (p->coups);
  free (p->nbors);
  free (p);
}

void
spnr_params_lr_free (spnr_params_t *p)
{
  free (p->coups);
  free (p);
}
