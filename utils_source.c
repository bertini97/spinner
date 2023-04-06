/* utils_source.c
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

#include "macros.h"
#include "couplings.h"

void
FUNC (coups_cubicnn_set, NAME) (float * const coups,
                                size_t * const nbors,
                                size_t const side,
                                size_t const n_dims)
{
  size_t i, j, left;
  size_t const n_inters = n_dims * 2;
  size_t const size = pow (side, n_dims);
  
  for (i = 0; i < size; ++i)
    {
      for (j = 0; j < n_dims; ++j)
        {
          left = nbors[i * n_inters + j];
          coups[i * n_inters + j] = FUNC (coups_get, NAME) ();
          coups[left * n_inters + j + n_dims] = coups[i * n_inters + j];
        }
    }
}

void
FUNC (coups_longrange_set, NAME) (float * const coups,
                                  size_t const size)
{
  size_t i, j;
  
  for (i = 0; i < size; ++i)
    {
      for (j = 0; j < size; ++j)
        {
          coups[i * size + j] = +SPNR_J;
          coups[j * size + i] = coups[i * size + j];
        }
    }
}
