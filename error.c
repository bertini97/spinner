/* error.c
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

#include <stdio.h>
#include <stdlib.h>

#include "error.h"

void
spnr_err (int err, char const *mess)
{
  fprintf (stderr, "ERROR: %s\n", mess);
  exit (err);
}

void
spnr_warn (int warn, char const *mess)
{
  fprintf (stderr, "WARNING: %s\n", mess);
}

void *
malloc_err (size_t const size)
{
  void * const p = malloc (size);
  if (!p)
    spnr_err (SPNR_ERROR_ALLOC, "malloc returned NULL");
  return p;
}