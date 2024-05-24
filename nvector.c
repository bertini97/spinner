/* nvector.c
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

typedef float spin_t;

typedef struct
{
  size_t n_comps;
  spin_t *spins;
} nvector_priv_t;

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
spins_set_up (spnr_graph_t *graph, void *priv)
{
  size_t i, j, size = graph->size, n_comps = graph->param;
  spin_t *spins = ((nvector_priv_t *) priv)->spins;
  
  for (i = 0; i < size; ++i)
    {
      spins[i * n_comps + 0] = 1.0;
      for (j = 1; j < n_comps; ++j)
        spins[i * n_comps + j] = 0.0;
    }
}

static void
spins_set_rand (spnr_graph_t *graph, void *priv)
{
  size_t i, j, size = graph->size, n_comps = graph->param;
  spin_t *spins = ((nvector_priv_t *) priv)->spins;
  
  for (i = 0; i < size; ++i)
      spin_rand (spins + i * n_comps, n_comps);
}

static void *
priv_alloc (spnr_graph_t *graph)
{
  nvector_priv_t *priv = malloc_err (sizeof (nvector_priv_t));
  
  if (!graph->param)
    spnr_err (SPNR_ERROR_PARAM_OOB, "param cannot be zero for nvector");
  priv->n_comps = graph->param;
  priv->spins = malloc_err (graph->size * priv->n_comps * sizeof (spin_t));
  
  spins_set_up (graph, priv);
  
  return priv;
}

static void
priv_free (void *priv)
{
  nvector_priv_t *priv_c = (nvector_priv_t *) priv;
  
  free (priv_c->spins);
  free (priv_c);
}

static void
spins_print_2d (spnr_graph_t *graph, void *priv)
{
  size_t i, j, index;
  nvector_priv_t *priv_c = priv;
  
  size_t const side = graph->side;
  spin_t *spins = priv_c->spins;
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
spins_print_3d (spnr_graph_t *graph, void *priv)
{
  size_t i, j, k, index;
  nvector_priv_t *priv_c = priv;
  size_t const side = graph->side;
  size_t const n_comps = priv_c->n_comps;
  spin_t *spins = priv_c->spins;
  
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
cn_calc_h (spnr_graph_t *graph, void *priv)
{
  size_t i, j, index;
  nvector_priv_t *priv_c = (nvector_priv_t *) priv;
  
  size_t const size = graph->size;
  size_t const n_dims = graph->n_dims;
  size_t const n_inters = graph->n_inters;
  size_t * const nbors = graph->nbors, *nbors_i;
  float * const coups = graph->coups, *coups_i;
  
  size_t const n_comps = priv_c->n_comps;
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
calc_m (spnr_graph_t *graph, void *priv)
{
  size_t i, j;
  nvector_priv_t *priv_c = (nvector_priv_t *) priv;
  
  size_t const size = graph->size;
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
cn_mcstep_metr (spnr_graph_t *graph, void *priv, float const beta)
{
  size_t i, k, index;
  nvector_priv_t *priv_c = (nvector_priv_t *) priv;
  
  size_t const size = graph->size;
  size_t const n_inters = graph->n_inters;
  size_t * const nbors = graph->nbors;
  float * const coups = graph->coups;
  
  size_t const n_comps = priv_c->n_comps;
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

static const spnr_kind_t cn_kind =
{
  "ising cubic nn ferr",
  &priv_alloc,
  &priv_free,
  &spins_set_up,
  &spins_set_rand,
  &spins_print_2d,
  &spins_print_3d,
  &cn_calc_h,
  &calc_m,
  &cn_mcstep_metr,
  NULL
};

const spnr_kind_t *spnr_nvector_cubicnn = &cn_kind;
