/* spinner.c
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
#include "common.h"

spnr_latt_t *
spnr_latt_alloc (spnr_latt_kind_t const * const kind,
                 size_t const side,
                 size_t const n_dims,
                 size_t const param)
{
  spnr_latt_t * l = spnr_malloc (sizeof (spnr_latt_t));
  
  l->kind = kind;
  l->priv = kind->alloc (side, n_dims, param);
  
  return l;
}

void
spnr_latt_free (spnr_latt_t * const l)
{
  l->kind->free (l->priv);
  free (l);
}


void
spnr_latt_print_spins_2d (spnr_latt_t const * const l)
{
  l->kind->print_spins_2d (l->priv);
}

void
spnr_latt_print_spins_3d (spnr_latt_t const * const l)
{
  l->kind->print_spins_3d (l->priv);
}

float
spnr_latt_calc_h (spnr_latt_t const * const l)
{
  float h;
  h = l->kind->calc_h (l->priv);
  return h; 
}

float
spnr_latt_calc_m (spnr_latt_t const * const l)
{
  float h;
  h = l->kind->calc_m (l->priv);
  return h; 
}

spnr_func_t *
spnr_metr (spnr_latt_t * const l)
{
  if (!l->kind->mcstep_metr)
    spnr_err (SPNR_ERROR_FUNC_NULL,
              "ERROR: func not found for this lattice kind");
  return l->kind->mcstep_metr;
}

spnr_func_t *
spnr_wolff (spnr_latt_t * const l)
{
  if (!l->kind->mcstep_wolff)
    spnr_err (SPNR_ERROR_FUNC_NULL,
              "ERROR: func not found for this lattice kind");
  return l->kind->mcstep_wolff;
}


spnr_data_t *
spnr_latt_run (spnr_func_t * (*getter) (spnr_latt_t *l),
               spnr_latt_t * const l, size_t const n_steps,
               size_t const n_probes, float const temp)
{
  size_t i, j;
  size_t const n_steps_bef_probe = n_steps / n_probes;
  float const beta = 1.0 / temp;
  
  void * priv = l->priv;
  float (*calc_h) (void *priv) = l->kind->calc_h;
  float (*calc_m) (void *priv) = l->kind->calc_m;
  void (*func) (void *priv, float beta) = getter(l);
  
  spnr_data_t * const data = spnr_data_alloc (n_probes + 1);
  
  srand (time (0));
  
  data->h[0] = calc_h (priv);
  data->m[0] = calc_m (priv);
  
  for (i = 1; i <= n_probes; ++i)
    {
      for (j = 0; j < n_steps_bef_probe; ++j)
          func (priv, beta);
      data->h[i] = calc_h (priv);
      data->m[i] = calc_m (priv);
    }
  
  return data;
}

spnr_data_t *
spnr_data_alloc (size_t const size)
{
  spnr_data_t *d = spnr_malloc (sizeof (spnr_data_t));
  
  d->size = size;
  d->h = malloc (d->size * sizeof (float));
  d->m = malloc (d->size * sizeof (float));
  
  return d;
}

void
spnr_data_free (spnr_data_t *d)
{
  free (d->h);
  free (d->m);
  free (d);
}

void
spnr_data_write (spnr_data_t const * const d,
                 char * const fname)
{
  FILE *f;
  size_t i;
  size_t const size = d->size;
  float * const h = d->h;
  float * const m = d->m;
  
  f = fopen (fname, "w");
  if (!f)
    exit (EXIT_FAILURE);
  
  for (i = 0; i < size; ++i)
    fprintf (f, "%05ld %+.3f %+.3f\n", i, h[i], m[i]);
  
  fclose (f);
}

void
spnr_data_mean_calc (spnr_data_t const * const d,
                     float *h,
                     float *m)
{
  size_t i;
  size_t const size = d->size;
  float *dh = d->h, *dm = d->m;
  float tmp_h, tmp_m;
  
  tmp_h = 0;
  tmp_m = 0;
  for (i = 0; i < size; ++i)
    {
      tmp_h += dh[i];
      tmp_m += dm[i];
    }
  
  *h = tmp_h / size;
  *m = tmp_m / size;
}

static void
spnr_corr_calc (float * const corr,
                float const * const arr,
                size_t const size)
{
  size_t i, j;
  size_t max;
  float a, b, c;
  
  for (i = 0; i < size; ++i)
    {
      max = size - i;
      a = 0;
      b = 0;
      c = 0;
      
      for (j = 0; j < max; ++j)
        {
          a += arr[j] * arr[j + i];
          b += arr[j];
          c += arr[j + i];
        }
      
      corr[i] = (a - b * c / max) / max;
  }
}

spnr_data_t *
spnr_data_corr_calc (spnr_data_t const * const data)
{
  size_t size = data->size;
  spnr_data_t * corr = spnr_data_alloc (size);
  spnr_corr_calc (corr->h, data->h, size);
  spnr_corr_calc (corr->m, data->m, size);
  
  return corr;
}
