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
#include "error.h"

void *
malloc_err (size_t const size)
{
  void * p = malloc (size);
  if (!p)
    spnr_err (SPNR_ERROR_ALLOC, "malloc returned NULL");
  return p;
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
