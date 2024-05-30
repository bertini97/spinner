/* ising.c
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

#include "spinner.h"
#include "error.h"

typedef int spin_t;

static void
set_up (void *priv, size_t N)
{
  spin_t *spins = (spin_t *)priv;
  size_t i;
  for (i = 0; i < N; ++i)
    spins[i] = +1.0;
}

static void
set_rand (void *priv, size_t N)
{
  spin_t *spins = (spin_t *)priv;
  size_t i;
  for (i = 0; i < N; ++i)
      spins[i] = ((float) rand() / RAND_MAX < 0.5) ? +1 : -1;
}

static void *
priv_alloc (size_t N, size_t _)
{
  spin_t *spins = malloc_err (N * sizeof(spin_t));
  set_up (spins, N);
  
  return spins;
}

static void
priv_free (void *priv)
{
  free (priv);
}

static size_t
spin_size (void *priv)
{
  return sizeof (spin_t);
}

static void
fill_prop (void const * const priv, void * const prop, size_t const k)
{
  spin_t const * const spins = (spin_t*) priv;
  spin_t * const prop_ = (spin_t*) prop;
  *prop_ = -spins[k];
}

static void
accept_prop (void * const priv, void const * const prop, size_t const k)
{
  spin_t * const spins = (spin_t*) priv;
  spin_t const * const prop_ = (spin_t*) prop;
  spins[k] = *prop_;
}

static float
calc_delta_h_binary (void const * const priv,
                     size_t const n_sites,
                     float const * const J,
                     size_t const * const sites,
                     void const * const prop,
                     size_t const k)
{
  spin_t const * const spins = (spin_t*) priv;
  float h = 0;
  size_t i;
  
  for (i = 0; i < n_sites; ++i)
    h += J[i] * spins[sites[i]];
  
  //printf ("-h * -2 = %f\n", -h * -2);
  return -2 * spins[k] * -h;
}

static float
calc_part_h_binary (void const * const priv, size_t const n_sites,
                    float const * const J, size_t const * const sites,
                    size_t const k)
{
  spin_t const * const spins = (spin_t*) priv;
  float h = 0;
  size_t i;
  
  for (i = 0; i < n_sites; ++i)
      h += J[i] * spins[k] * spins[sites[i]];
  
  return -h;
}

static float
calc_phi (void const * const priv, size_t const N)
{
  spin_t const * const spins = (spin_t*) priv;
  size_t i;
  float m = 0;
  
  for (i = 0; i < N; ++i)
    m += spins[i];
  
  return m / (float) N;
}

static const spnr_sys_kind_t ising_kind =
{
  "ising",
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

const spnr_sys_kind_t *spnr_ising = &ising_kind;