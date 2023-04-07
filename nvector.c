/* ising.c
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
#include "utils.h"

typedef float spin_t;

typedef struct
{
  size_t n_dims;
  size_t side;
  
  size_t size;
  size_t n_comps;
  spin_t *spins;
  
  size_t n_inters;
  size_t *nbors;
  float *coups;
} nvector_priv_t;

static void
spins_set_up (spin_t *spins, size_t size, size_t n_comps)
{
  size_t i, j;
  
  for (i = 0; i < size; ++i)
    {
      spins[i * n_comps + 0] = 1.0;
      for (j = 1; j < n_comps; ++j)
        spins[i * n_comps + j] = 0.0;
    }
}

static void *
cnf_alloc (size_t const side, size_t const n_dims,
                         size_t const param)
{
  nvector_priv_t *priv = spnr_malloc (sizeof (nvector_priv_t));
  
  if (side <= 0 || side >= SPNR_SIDE_MAX
      || n_dims <= 0 || n_dims >= SPNR_DIMS_MAX
      || param <= 0 || param >= SPNR_COMPS_MAX)
    spnr_err (SPNR_ERROR_PARAM_OOB, "lattice parameters out of bounds");
  
  priv->n_dims = n_dims;
  priv->side = side;
  priv->n_comps = param;
  priv->size = pow (side, n_dims);
  priv->n_inters = n_dims * 2;
  
  priv->spins = spnr_malloc (priv->size * priv->n_comps * sizeof (spin_t));
  priv->nbors = spnr_malloc (priv->size * priv->n_inters * sizeof (size_t));
  priv->coups = spnr_malloc (priv->size * priv->n_inters * sizeof (float));
  
  spins_set_up (priv->spins, priv->size, priv->n_comps);
  nbors_cubicnn_set (priv->nbors, side, n_dims);
  coups_cubicnn_set_ferr (priv->coups, priv->nbors, side, n_dims);
  
  return priv;
}

static void
cnf_free (void *priv)
{
  nvector_priv_t *priv_c = (nvector_priv_t *) priv;
  
  free (priv_c->spins);
  free (priv_c->nbors);
  free (priv_c->coups);
  
  free (priv_c);
}

static float
rand_gauss ()
{
  static int has_prev = 0;
  static float prev_rand;

  if (has_prev)
    {
      has_prev = 0;
      return prev_rand;
    }
  else
    {
      float const u = (float) rand () / RAND_MAX;
      float const v = (float) rand () / RAND_MAX;
      float const x = sqrt (-2. * log (u)) * sin (2 * SPNR_PI * v);
      float const y = sqrt (-2. * log (u)) * cos (2 * SPNR_PI * v);

      has_prev = 1;
      prev_rand = x;
      return y;
    }
}

static void
spin_rand (spin_t *u, const size_t n_comps)
{
  size_t i;
  float rand_g, norm = 0;

  for (i = 0; i < n_comps; ++i)
    {
      rand_g = rand_gauss ();
      u[i] = rand_g;
      norm += rand_g * rand_g;
    }
  norm = sqrt(norm);

  for (i = 0; i < n_comps; ++i)
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
spin_sub (spin_t *s, spin_t * const u, spin_t * const v,
          size_t const n_comps)
{
  size_t i;
  
  for (i = 0; i < n_comps; ++i)
    s[i] = u[i] - v[i];
}

static float
spin_mod (spin_t * const u, size_t const n_comps)
{
  size_t i;
  float mag = 0;
  
  for (i = 0; i < n_comps; ++i)
    mag += u[i] * u[i];
  
  return sqrt(mag);
}

static void
spin_print (spin_t * const u, size_t const n_comps)
{
  size_t i;
  
  printf ("[");
  for (i = 0; i < n_comps; ++i)
    printf ("%+.1f ", u[i]);
  printf ("]");
}

static void
print_spins_2d (void *priv)
{
  size_t i, j, index;
  nvector_priv_t *priv_c = priv;
  spin_t *spins = priv_c->spins;
  size_t const side = priv_c->side;
  size_t const n_comps = priv_c->n_comps;
  
  for (i = 0; i < side; i++)
    {
      for (j = 0; j < side; j++)
        {
          index = i * side + j;
          spin_print (spins + index * n_comps, n_comps);
        }
      printf ("\n");
    }
  printf ("\n");
}

static void
print_spins_3d (void *priv)
{
  size_t i, j, k, index;
  nvector_priv_t *priv_c = priv;
  spin_t *spins = priv_c->spins;
  size_t const side = priv_c->side;
  size_t const n_comps = priv_c->n_comps;
  
  for (i = 0; i < side; i++)
    {
      for (j = 0; j < side; j++)
        {
          for (k = 0; k < side; k++)
            {
              index = i * side * side + j * side + k;
              spin_print (spins + index * n_comps, n_comps);
            }
          printf ("\n");
        }
      printf ("\n");
    }
}

static float
cnf_calc_h (void *priv)
{
  size_t i, j, index;
  nvector_priv_t *priv_c = (nvector_priv_t *) priv;
  
  size_t const size = priv_c->size;
  size_t const n_dims = priv_c->n_dims;
  size_t const n_comps = priv_c->n_comps;
  size_t const n_inters = priv_c->n_inters;
  size_t * const nbors = priv_c->nbors, *nbors_i;
  float * const coups = priv_c->coups, *coups_i;
  spin_t * const spins = priv_c->spins, *spins_i;
  
  float h = 0;
  
  for (i = 0; i < size; ++i)
    {
      spins_i = spins + i * n_comps;
      index = i * n_inters;
      nbors_i = nbors + index;
      coups_i = coups + index;
      
      for (j = 0; j < n_dims; ++j)
        h += coups_i[j] *
          spin_sprod(spins_i, spins + nbors_i[j] * n_comps, n_comps);
    }
  h *= - 2.0 / (float) size;
  
  return h;
}

static float
calc_m (void *priv)
{
  size_t i, j;
  nvector_priv_t *priv_c = (nvector_priv_t *) priv;
  
  size_t const size = priv_c->size;
  size_t const n_comps = priv_c->n_comps;
  spin_t * const spins = priv_c->spins, *spins_i;
  
  spin_t m[SPNR_COMPS_MAX] = {0};
  
  for (i = 0; i < size; ++i)
    {
      spins_i = spins + i * n_comps;
      for (j = 0; j < n_comps; ++j)
        m[j] += spins_i[j];
    }
  
  return spin_mod (m, n_comps) / (float) size;
}

static float
h_delta_calc (size_t * const nbors_k, float * const coups_k,
              spin_t * const spin_delta, spin_t * const spins,
              size_t const n_inters, size_t const n_comps)
{
  size_t i;
  float h = 0;
  
  for (i = 0; i < n_inters; ++i)
    h += coups_k[i] *
      spin_sprod (spin_delta, spins + nbors_k[i] * n_comps, n_comps);
  h *= -1.0;
  
  return h;
}

static void
cnf_mcstep_metr (void *priv, float const beta)
{
  size_t i, k, index;
  nvector_priv_t *priv_c = (nvector_priv_t *) priv;
  
  size_t const size = priv_c->size;
  size_t const n_comps = priv_c->n_comps;
  size_t const n_inters = priv_c->n_inters;
  size_t * const nbors = priv_c->nbors;
  float * const coups = priv_c->coups;
  spin_t * const spins = priv_c->spins, *spins_k;
  
  spin_t spin_prop[SPNR_COMPS_MAX], spin_delta[SPNR_COMPS_MAX];
  float h_delta;
  
  for (i = 0; i < size; ++i)
    {
      k = rand() % size;
      spins_k = spins + k * n_comps;
      
      spin_rand (spin_prop, n_comps);
      spin_sub (spin_delta, spin_prop, spins_k, n_comps);
      
      index = k * n_inters;
      h_delta = h_delta_calc (nbors + index, coups + index,
                              spin_delta, spins,
                              n_inters, n_comps);
      
      if (metr_prop_accept (h_delta, beta))
        memcpy (spins_k, spin_prop, n_comps * sizeof (spin_t));
    }
}

static void
cnf_mcstep_wolff (void *priv, float beta)
{
  /* TODO */
}

static const spnr_latt_kind_t cnf_kind =
{
  "ising cubic nn ferr",
  &cnf_alloc,
  &cnf_free,
  &print_spins_2d,
  &print_spins_3d,
  &cnf_calc_h,
  &calc_m,
  &cnf_mcstep_metr,
  NULL
};

const spnr_latt_kind_t *spnr_nvector_cubicnn_ferr = &cnf_kind;
