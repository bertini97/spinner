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

spnr_latt_t *
spnr_latt_alloc (spnr_kind_t const * kind,
                 spnr_params_t * params)
{
  spnr_latt_t * l = malloc_err (sizeof (spnr_latt_t));
  
  if (!kind || !params)
    spnr_err (SPNR_ERROR_ARG_NULL, "requested argument is NULL");
  l->kind = kind;
  l->params = params;
  l->priv = kind->priv_alloc (params);
  
  return l;
}

void
spnr_latt_free (spnr_latt_t * const l)
{
  l->kind->priv_free (l->priv);
  free (l);
}

void
spnr_latt_spins_set_up (spnr_latt_t * const l)
{
  l->kind->spins_set_up (l->params, l->priv);
}

void
spnr_latt_spins_set_rand (spnr_latt_t * const l)
{
  l->kind->spins_set_rand (l->params, l->priv);
}

void
spnr_latt_spins_print_2d (spnr_latt_t const * const l)
{
  l->kind->spins_print_2d (l->params, l->priv);
}

void
spnr_latt_spins_print_3d (spnr_latt_t const * const l)
{
  l->kind->spins_print_3d (l->params, l->priv);
}

float
spnr_latt_calc_h (spnr_latt_t const * const l)
{
  float h;
  h = l->kind->calc_h (l->params, l->priv);
  return h; 
}

float
spnr_latt_calc_m (spnr_latt_t const * const l)
{
  float h;
  h = l->kind->calc_m (l->params, l->priv);
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
spnr_latt_run (spnr_func_getter_t *getter, spnr_latt_t * l,
               size_t n_steps, size_t n_probes, float temp)
{
  size_t i, j;
  size_t const n_steps_bef_probe = n_steps / n_probes;
  float const beta = 1.0 / temp;
  
  spnr_params_t *params = l->params;
  void * priv = l->priv;
  float (*calc_h) (spnr_params_t *params, void *priv) = l->kind->calc_h;
  float (*calc_m) (spnr_params_t *params, void *priv) = l->kind->calc_m;
  void (*func) (spnr_params_t *params, void *priv, float beta) = getter(l->kind);
  
  spnr_data_t * const data = spnr_data_alloc (n_probes + 1);
  
  srand (time (0));
  
  data->h[0] = calc_h (params, priv);
  data->m[0] = calc_m (params, priv);
  
  for (i = 1; i <= n_probes; ++i)
    {
      for (j = 0; j < n_steps_bef_probe; ++j)
          func (l->params, priv, beta);
      data->h[i] = calc_h (params, priv);
      data->m[i] = calc_m (params, priv);
    }
  
  return data;
}
