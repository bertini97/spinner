/* utils.c
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

#include "utils.h"

#define NAME  ferr
#include "utils_source.c"
#undef  NAME

/*#define	NAME	bim
#include "utils_source.c"
#undef	NAME*/

void
nbors_cubicnn_set (size_t * const nbors,
                   size_t const side,
                   size_t const n_dims)
{
  size_t i, j;
  size_t slices[SPNR_DIMS_MAX], size;
  size_t unit, row, last_row, left;
  
  size_t const n_inters = n_dims * 2;

  for (i = 0; i <= n_dims; ++i)
    slices[i] = pow (side, i);
  size = slices[n_dims];
  
  for (i = 0; i < size; ++i)
    {
      for (j = 0; j < n_dims; ++j)
        {
          unit = slices[j];
          row  = slices[j+1];
          last_row = row - unit;
          left = ((i % row) < unit) ? (i + last_row) : (i - unit);
      
          nbors[i * n_inters + j] = left;
          nbors[left * n_inters + j + n_dims] = i;
        }
    }
}


int
metr_prop_accept (float const h_delta,
                  float const beta)
{
  if (h_delta <= 0)
    return SPNR_TRUE;
  else
    {
      float acc_ratio = exp (- beta * h_delta);
      float rand_num = (float) rand () / RAND_MAX;
      if (rand_num < acc_ratio)
        return SPNR_TRUE;
      else
        return SPNR_FALSE;
  }
}
