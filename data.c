/* data.c
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

#include "spinner.h"

spnr_data_t *
spnr_data_alloc (size_t const size)
{
  spnr_data_t * const data = malloc (sizeof (spnr_data_t));
  
  data->size = size;
  data->h = malloc (data->size * sizeof (float));
  data->phi = malloc (data->size * sizeof (float));
  
  return data;
}

void
spnr_data_free (spnr_data_t * const data)
{
  free (data->h);
  free (data->phi);
  free (data);
}

void
spnr_data_write (spnr_data_t const * const data,
                 char const * const fname)
{
  FILE *f;
  size_t i;
  size_t const size = data->size;
  float * const h = data->h;
  float * const phi = data->phi;
  
  f = fopen (fname, "w");
  if (!f)
    exit (EXIT_FAILURE);
  
  for (i = 0; i < size; ++i)
    fprintf (f, "%05ld %+.3f %+.3f\n", i, h[i], phi[i]);
  
  fclose (f);
}

void
spnr_data_mean_calc (spnr_data_t const * const data,
                     float * const h_mean,
                     float * const phi_mean)
{
  size_t i;
  size_t const size = data->size;
  float * const h = data->h;
  float * const phi = data->phi;
  float sum_h = 0;
  float sum_phi = 0;
  
  for (i = 0; i < size; ++i)
    {
      sum_h += h[i];
      sum_phi += phi[i];
    }
  
  *h_mean = sum_h / size;
  *phi_mean = sum_phi / size;
}

void
spnr_data_var_calc (spnr_data_t const * const data,
                    float * const h_var,
                    float * const phi_var)
{
  size_t i;
  size_t const size = data->size;
  float * const h = data->h;
  float * const phi = data->phi;
  float sum_sq_h = 0;
  float sum_sq_phi = 0;
  float sum_h = 0;
  float sum_phi = 0;
  
  for (i = 0; i < size; ++i)
    {
      sum_sq_h += h[i] * h[i];
      sum_sq_phi += phi[i] * phi[i];
      sum_h += h[i];
      sum_phi += phi[i];
    }
  
  *h_var = sum_sq_h / size - (sum_h * sum_h) / (size * size);
  *phi_var = sum_sq_phi / size - (sum_phi * sum_phi) / (size * size);
}

static void
spnr_corr_calc (float * const corr, float const * const arr, size_t const size)
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

void
spnr_data_corr_calc (spnr_data_t * const corr,
                     spnr_data_t const * const data)
{
  if (corr->size != data->size)
    exit (EXIT_FAILURE);
  spnr_corr_calc (corr->h, data->h, corr->size);
  spnr_corr_calc (corr->phi, data->phi, corr->size);
}

void
spnr_data_run_and_probe (spnr_data_t * const data,
                         spnr_sys_t * const sys,
                         spnr_step_t const * const step,
                         float const temp,
                         size_t const n_steps_before_probe)
{
  size_t i, j;
  size_t const n_probes = data->size;
  float const beta = 1.0 / temp;
  float * const h = data->h;
  float * const phi = data->phi;
  
  h[0] = spnr_sys_calc_h (sys);
  phi[0] = spnr_sys_calc_phi (sys);
  for (i = 1; i < n_probes; ++i)
  {
    for (j = 0; j < n_steps_before_probe; ++j)
      spnr_step_apply (step, sys, beta);
    h[i] = spnr_sys_calc_h (sys);
    phi[i] = spnr_sys_calc_phi (sys);
  }
}