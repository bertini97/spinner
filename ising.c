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
#include "utils.h"

typedef int_fast8_t spin_t;

typedef struct
{
  spin_t *spins;
  
  size_t *buffer;
  int *in_cluster;
} ising_priv_t;

static void
spins_set_up (spnr_graph_t *graph, void *priv)
{
  spin_t *spins = ((ising_priv_t *) priv)->spins;
  memset (spins, +1, graph->size * sizeof (spin_t));
}

static void
spins_set_rand (spnr_graph_t *graph, void *priv)
{
  size_t i, size = graph->size;
  spin_t *spins = ((ising_priv_t *) priv)->spins;
  
  for (i = 0; i < size; ++i)
      spins[i] = ((float) rand() / RAND_MAX < 0.5) ? +1 : -1;
}

static void *
priv_alloc (spnr_graph_t * const graph)
{
  ising_priv_t *priv = malloc_err (sizeof (ising_priv_t));
  priv->spins = malloc_err (graph->size * sizeof (spin_t));
  spins_set_up (graph, priv);
  
  priv->buffer = malloc_err (graph->size * sizeof (size_t));
  priv->in_cluster = malloc_err (graph->size * sizeof (int));
  
  return priv;
}

static void
priv_free (void *priv)
{
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  free (priv_c->spins);
  
  free (priv_c->buffer);
  free (priv_c->in_cluster);
  
  free (priv_c);
}

static void
spins_print_2d (spnr_graph_t *graph, void *priv)
{
  size_t i, j;
  ising_priv_t *priv_c = priv;
  spin_t *spins = priv_c->spins;
  size_t const side = graph->side;
  
  for (i = 0; i < side; i++)
    {
      for (j = 0; j < side; j++)
        printf ("%+d ", spins[i * side + j]);
      printf ("\n");
    }
  printf ("\n");
}


static void
spins_print_3d (spnr_graph_t *graph, void *priv)
{
  size_t i, j, k;
  ising_priv_t *priv_c = priv;
  spin_t *spins = priv_c->spins;
  size_t const side = graph->side;
  
  for (i = 0; i < side; i++)
    {
      for (j = 0; j < side; j++)
        {
          for (k = 0; k < side; k++)
            printf ("%+d ", spins[i * side * side + j * side + k]);
          printf ("\n");
        }
      printf ("\n");
    }
}


static float
cn_calc_h (spnr_graph_t *graph, void *priv)
{
  size_t i, j, i_inters;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = graph->size;
  size_t const n_dims = graph->n_dims;
  size_t const n_inters = graph->n_inters;
  size_t * const nbors = graph->nbors, *nbors_i;
  float * const coups = graph->coups, *coups_i;
  spin_t * const spins = priv_c->spins;
  
  float h = 0;
  
  for (i = 0; i < size; ++i)
    {
      i_inters = i * n_inters;
      nbors_i = nbors + i_inters;
      coups_i = coups + i_inters;
      for (j = 0; j < n_dims; ++j)
        h += coups_i[j] * spins[i] * spins[nbors_i[j]];
    }
  h *= - 2.0 / size;
  
  return h;
}


static float
lr_calc_h (spnr_graph_t *graph, void *priv)
{
  size_t i, j;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = graph->size;
  float * const coups = graph->coups, *coups_i;
  spin_t * const spins = priv_c->spins;
  
  float h = 0;
  
  for (i = 0; i < size; ++i)
    {
      coups_i = coups + i * size;
      for (j = 0; j < size; ++j)
        {
          if (j == i)
            continue;
          h += coups_i [j] * spins[i] * spins[j];
        }
    }
  h *= - 2.0 / (size * sqrt(size));
  
  return h;
}


static float
calc_m (spnr_graph_t *graph, void *priv)
{
  size_t i;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = graph->size;
  spin_t * const spins = priv_c->spins;
  
  float m = 0;
  
  for (i = 0; i < size; ++i)
    m += spins[i];
  m /= size;
  
  return m;
}

static void
cn_mcstep_metr (spnr_graph_t *graph, void *priv, float const beta)
{
  size_t i, j, k, k_inters;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = graph->size;
  size_t const n_inters = graph->n_inters;
  size_t * const nbors = graph->nbors, *nbors_k;
  float * const coups = graph->coups, *coups_k;
  spin_t * const spins = priv_c->spins, *spins_k;
  
  float h_delta;
  
  for (i = 0; i < size; ++i)
    {
      k = rand() % size;
      k_inters = k * n_inters;
      spins_k = spins + k;
      nbors_k = nbors + k_inters;
      coups_k = coups + k_inters;
      
      h_delta = 0;
      for (j = 0; j < n_inters; ++j)
        h_delta += coups_k[j] * spins[nbors_k[j]];
      h_delta *= 2 * *spins_k;
      
      if (metr_prop_accept (h_delta, beta))
        *spins_k *= -1;
    }
}

static void
lr_mcstep_metr (spnr_graph_t *graph, void *priv, float const beta)
{
  size_t i, j, k;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = graph->size;
  float const one_over_sqsize = 1.0 / sqrt(size);
  float * const coups = graph->coups, *coups_k;
  spin_t * const spins = priv_c->spins, *spins_k;
  
  float h_delta;
  
  for (i = 0; i < size; ++i)
    {
      k = rand() % size;
      spins_k = spins + k;
      coups_k = coups + k * size;
      
      h_delta = 0;
      for (j = 0; j < k; ++j)
        h_delta += coups_k[j] * spins[j];
      for (j = k + 1; j < size; ++j)
        h_delta += coups_k[j] * spins[j];
      h_delta *= 2 * *spins_k * one_over_sqsize;
      
      if (metr_prop_accept (h_delta, beta))
        *spins_k *= -1;
    }
}

void
cn_mcstep_wolff (spnr_graph_t *graph, void *priv, float const beta)
{
  size_t i, k, head, index_cur, index_nbor;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = graph->size;
  size_t const n_inters = graph->n_inters;
  size_t * const nbors = graph->nbors, *nbors_cur;
  float * const coups = graph->coups, *coups_cur;
  
  spin_t * const spins = priv_c->spins, spin_cur;
  size_t * const buffer = priv_c->buffer;
  int * const in_cluster = priv_c->in_cluster;
  
  float p_add, rand_num;
  
  head = 0;
  memset (in_cluster, SPNR_FALSE, size * sizeof (int));
  
  k = rand() % size;
  spin_cur = spins[k];
  nbors_cur = nbors + k * n_inters;
  
  buffer[head++] = k;
  in_cluster[k] = SPNR_TRUE;
  
  while (head)
    {
      index_cur = buffer[--head];
      spin_cur = spins[index_cur];
      nbors_cur = nbors + index_cur * n_inters;
      coups_cur = coups + index_cur * n_inters;
      
      for (i = 0; i < n_inters; ++i)
        {
          index_nbor = nbors_cur[i];
          if (!in_cluster[index_nbor] && spin_cur == spins[index_nbor])
            {
              p_add = 1.0 - exp (-2.0 * beta * coups_cur[i]
                                 * spins[nbors_cur[i]] * spin_cur);
              rand_num = (double) rand () / RAND_MAX;
              if (rand_num < p_add)
                {
                  buffer[head++] = index_nbor;
                  in_cluster[index_nbor] = SPNR_TRUE;
                }
            }
        }
      spins[index_cur] *= -1;
    }
}

static const spnr_kind_t cn_kind =
{
  "icn",
  &priv_alloc,
  &priv_free,
  &spins_set_up,
  &spins_set_rand,
  &spins_print_2d,
  &spins_print_3d,
  &cn_calc_h,
  &calc_m,
  &cn_mcstep_metr,
  &cn_mcstep_wolff,
};

const spnr_kind_t *spnr_ising_cubicnn = &cn_kind;

static const spnr_kind_t lr_kind =
{
  "ilr",
  &priv_alloc,
  &priv_free,
  &spins_set_up,
  &spins_set_rand,
  &spins_print_2d,
  &spins_print_3d,
  &lr_calc_h,
  &calc_m,
  &lr_mcstep_metr,
  NULL
};

const spnr_kind_t *spnr_ising_longrange = &lr_kind;
