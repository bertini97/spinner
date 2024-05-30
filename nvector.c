/* nvector.c
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

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "spinner.h"
#include "error.h"

#define SPNR_NVECTOR_D_MAX  8
#define SPNR_PI             3.1415926536

typedef float spin_t;

typedef struct
{
  size_t n;
  spin_t *spins;
}
nvector_priv_t;

static float
rand_gauss ()
{
  static int has_prev = 0;
  static float prev_rand;
  float u, v, x, y;

  if (has_prev)
    {
      has_prev = 0;
      return prev_rand;
    }
  else
    {
      u = (float) rand () / RAND_MAX;
      v = (float) rand () / RAND_MAX;
      x = sqrt (-2. * log (u)) * sin (2 * SPNR_PI * v);
      y = sqrt (-2. * log (u)) * cos (2 * SPNR_PI * v);

      has_prev = 1;
      prev_rand = x;
      return y;
    }
}

static void
spin_rand (spin_t * const u, const size_t n)
{
  size_t i;
  float rand_g, norm = 0;

  for (i = 0; i < n; ++i)
    {
      rand_g = rand_gauss ();
      u[i] = rand_g;
      norm += rand_g * rand_g;
    }
  norm = sqrt(norm);

  for (i = 0; i < n; ++i)
    u[i] /= norm;
}

static float
spin_sprod (spin_t * const u, spin_t * const v, size_t const n_comps)
{
  size_t i;
  float prod;

  prod = 0;
  for (i = 0; i < n_comps; ++i)
    prod += u[i] * v[i];

  return prod;
}

static void
spin_sub (spin_t *s, spin_t * const u, spin_t * const v, size_t const n)
{
  size_t i;
  
  for (i = 0; i < n; ++i)
    s[i] = u[i] - v[i];
}

static float
spin_mod (spin_t * const u, size_t const n)
{
  size_t i;
  float mag = 0;
  
  for (i = 0; i < n; ++i)
    mag += u[i] * u[i];
  
  return sqrt(mag);
}

static void
spin_print (spin_t * const u, size_t const n)
{
  size_t i;
  for (i = 0; i < n; ++i)
    printf ("%f ", u[i]);
}

static void
set_up (void * const priv, size_t const N)
{
  nvector_priv_t * const priv_ = (nvector_priv_t *)priv;
  spin_t *spins = priv_->spins;
  size_t i, j, n = priv_->n;
  
  for (i = 0; i < N; ++i)
    {
      spins[i * n] = +1.0;
      for (j = 1; j < n; ++j)
        spins[i * n + j] = 0.0;
    }
}

static void
set_rand (void * const priv, size_t const N)
{
  nvector_priv_t * const priv_ = (nvector_priv_t *)priv;
  spin_t *spins = priv_->spins;
  size_t i, j, n = priv_->n;
  
  for (i = 0; i < N; ++i)
    spin_rand (spins + i * n, n);
}

static void *
priv_alloc (size_t const N, size_t const n)
{
  if (n <= 0)
    spnr_err (SPNR_ERROR_PARAM_OOB, "number of components must be positive");
  
  nvector_priv_t * const priv = malloc (sizeof (nvector_priv_t));
  priv->n = n;
  priv->spins = malloc (N * n * sizeof (spin_t));
  
  set_up (priv, N);
  
  return priv;
}

static void
priv_free (void *priv)
{
  nvector_priv_t * const priv_ = (nvector_priv_t *)priv;
  free (priv_->spins);
  free (priv);
}

static size_t
spin_size (void *priv)
{
  nvector_priv_t *priv_ = (nvector_priv_t *)priv;
  return sizeof (priv_->n * sizeof(spin_t));
}

static void
fill_prop (void const * const priv, void * const prop, size_t const k)
{
  nvector_priv_t const * const priv_ = (nvector_priv_t*) priv;
  spin_t * const prop_ = (spin_t*) prop;
  spin_rand (prop, priv_->n);
}

static void
accept_prop (void * const priv, void const * const prop, size_t k)
{
  nvector_priv_t * const priv_ = (nvector_priv_t*) priv;
  spin_t const * const prop_ = (spin_t*) prop;
  size_t const n = priv_->n;
  spin_t * const  spin_k = priv_->spins + k * n;
  
  memcpy (spin_k, prop_, n * sizeof (spin_t));
}

static float
calc_delta_h_binary (void const * const priv,
                     size_t const n_sites,
                     float const * const J,
                     size_t const * const sites,
                     void const * const prop,
                     size_t const k)
{
  nvector_priv_t const * const priv_ = (nvector_priv_t*) priv;
  spin_t const * const spins = priv_->spins;
  size_t const n = priv_->n;
  spin_t const * const prop_ = (spin_t*) prop;
  spin_t const * const spin_k = spins + k * n;
  size_t i, j;
  float h;
  spin_t sum[SPNR_NVECTOR_D_MAX];
  
  memset (sum, 0, sizeof (sum));
  for (i = 0; i < n_sites; ++i)
    for (j = 0; j < n; ++j)
      sum[j] += J[i] * spins[sites[i] * n + j];
  
  h = 0;
  for (j = 0; j < n; ++j)
    h += (spin_k[j] - prop_[j]) * sum[j];
  
  return h;
}

static float
calc_part_h_binary (void const * const priv, size_t const n_sites,
                    float const * const J, size_t const * const sites,
                    size_t const k)
{
  nvector_priv_t const * const priv_ = (nvector_priv_t*) priv;
  spin_t *spins = priv_->spins;
  size_t i, j, n = priv_->n;
  spin_t *spin_k = spins + k * n;
  float h;
  spin_t sum[SPNR_NVECTOR_D_MAX];
  
  memset (sum, 0, sizeof (sum));
  for (i = 0; i < n_sites; ++i)
    for (j = 0; j < n; ++j)
      sum[j] += J[i] * spins[sites[i] * n + j];
  
  h = 0;
  for (j = 0; j < n; ++j)
    h += spin_k[j] * sum[j];
  
  return -h;
}

static float
calc_phi (void const * const priv, size_t const N)
{
  nvector_priv_t const * const priv_ = (nvector_priv_t*) priv;
  spin_t *spins = priv_->spins;
  size_t i, j, n = priv_->n;
  float m;
  spin_t sum[SPNR_NVECTOR_D_MAX];
  
  memset (sum, 0, sizeof (sum));
  for (i = 0; i < N; ++i)
    for (j = 0; j < n; ++j)
      sum[j] += spins[i * n + j];
  
  m = 0;
  for (j = 0; j < n; ++j)
    m += sum[j] * sum[j];
  
  return sqrt(m) / (float) N;
}

static const spnr_sys_kind_t nvector_kind =
{
  "nvector",
  &priv_alloc,
  &priv_free,
  &spin_size,
  &set_up,
  &set_rand,
  &fill_prop,
  &accept_prop,
  &calc_delta_h_binary,
  &calc_part_h_binary,
  &calc_phi
};

const spnr_sys_kind_t *spnr_nvector = &nvector_kind;