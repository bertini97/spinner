/* step.c
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

#include <time.h>

#include "spinner.h"
#include "error.h"

spnr_step_t *
spnr_step_alloc (spnr_step_kind_t const * const kind, size_t const param)
{
  srand (time(0));
  spnr_step_t * step = malloc_err (sizeof (spnr_step_t));
  step->kind = kind;
  step->priv = kind->priv_alloc (param);
  
  return step;
}

void
spnr_step_free (spnr_step_t * const step)
{
  step->kind->priv_free (step->priv);
  free (step);
}

void
spnr_step_apply (spnr_step_t const * const step,
                 spnr_sys_t const * const sys,
                 float const beta)
{
  step->kind->apply(step->priv, sys, beta);
}