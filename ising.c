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

typedef int_fast8_t spin_t;

typedef struct
{
  size_t n_dims;
  size_t side;
  
  size_t size;
  spin_t *spins;
  
  size_t n_inters;
  size_t *nbors;
  float *coups;
  
  size_t *buffer;
  int *in_cluster;
} ising_priv_t;

static void
spins_set_up (spin_t *spins, size_t size)
{
  memset (spins, +1, size * sizeof (spin_t));
}

static void *
cnf_alloc (size_t const side, size_t const n_dims,
           __attribute__((unused)) size_t const param)
{
  ising_priv_t *priv = spnr_malloc (sizeof (ising_priv_t));
  
  if (side <= 0 || side >= SPNR_SIDE_MAX
      || n_dims <= 0 || n_dims >= SPNR_DIMS_MAX)
    spnr_err (SPNR_ERROR_PARAM_OOB, "lattice parameters out of bounds");
  
  priv->n_dims = n_dims;
  priv->side = side;
  priv->size = pow (side, n_dims);
  priv->n_inters = n_dims * 2;
  
  priv->spins = spnr_malloc (priv->size * sizeof (spin_t));
  priv->nbors = spnr_malloc (priv->size * priv->n_inters * sizeof (size_t));
  priv->coups = spnr_malloc (priv->size * priv->n_inters * sizeof (float));
  
  spins_set_up (priv->spins, priv->size);
  nbors_cubicnn_set (priv->nbors, side, n_dims);
  coups_cubicnn_set_ferr (priv->coups, priv->nbors, side, n_dims);
  
  priv->buffer = spnr_malloc (priv->size * sizeof (size_t));
  priv->in_cluster = spnr_malloc (priv->size * sizeof (int));
  
  return priv;
}


/*static void *
ising_priv_ilrf_alloc (size_t const side, size_t const n_dims,
                       __attribute__((unused)) size_t const param)
{
  ising_priv_t *priv = spnr_malloc (sizeof (ising_priv_t));
  
  priv->size = pow (side, n_dims);
  priv->spins = malloc (priv->size * sizeof (spin_t));
  priv->coups = malloc (priv->size * priv->size * sizeof (float));
  
  spins_set_up (priv->spins, priv->size);
  coups_longrange_set_ferr (priv->coups, priv->size);
  
  return priv;
}*/


static void
cnf_free (void *priv)
{
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  free (priv_c->spins);
  free (priv_c->nbors);
  free (priv_c->coups);
  
  free (priv_c);
}


/*static void
lrf_free (void *priv)
{
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  free (priv_c->spins);
  free (priv_c->coups);
  
  free (priv_c);
}*/


static void
print_spins_2d (void *priv)
{
  size_t i, j;
  ising_priv_t *priv_c = priv;
  spin_t *spins = priv_c->spins;
  size_t const side = priv_c->side;
  
  for (i = 0; i < side; i++)
    {
      for (j = 0; j < side; j++)
        printf ("%+d ", spins[i * side + j]);
      printf ("\n");
    }
  printf ("\n");
}


static void
print_spins_3d (void *priv)
{
  size_t i, j, k;
  ising_priv_t *priv_c = priv;
  spin_t *spins = priv_c->spins;
  size_t const side = priv_c->side;
  
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
cnf_calc_h (void *priv)
{
  size_t i, j;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = priv_c->size;
  size_t const n_dims = priv_c->n_dims;
  size_t const n_inters = priv_c->n_inters;
  size_t * const nbors = priv_c->nbors;
  spin_t * const spins = priv_c->spins;
  
  float h = 0;
  
  for (i = 0; i < size; ++i)
    {
      for (j = 0; j < n_dims; ++j)
        h += spins[i] * spins[nbors[i * n_inters + j]];
    }
  h *= - 2.0 / size;
  
  return h;
}


/*static float
ising_priv_ilrf_calc_h (void *priv)
{
  size_t i, j;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = priv_c->size;
  spin_t * const spins = priv_c->spins;
  
  float h = 0;
  
  for (i = 0; i < size; ++i)
    {
      for (j = 0; j < size; ++j)
        {
          if (j == i)
            continue;
          h += spins[i] * spins[j];
        }
    }
  h *= - 2.0 / size;
  
  return h;
}*/


static float
calc_m (void *priv)
{
  size_t i;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = priv_c->size;
  spin_t * const spins = priv_c->spins;
  
  float m = 0;
  
  for (i = 0; i < size; ++i)
    m += spins[i];
  m /= size;
  
  return m;
}


static float
h_delta_calc (size_t * const nbors_k, float * const coups_k,
              spin_t * const spins, spin_t * const spins_k,
              size_t const n_inters)
{
  size_t i;
  float h = 0;
  
  h = 0;
  for (i = 0; i < n_inters; ++i)
    h += coups_k[i] * spins[nbors_k[i]];
  h *= 2 * *spins_k;
  
  return h;
}


static void
cnf_mcstep_metr (void *priv, float const beta)
{
  size_t i, k, index;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = priv_c->size;
  size_t const n_inters = priv_c->n_inters;
  size_t * const nbors = priv_c->nbors;
  float * const coups = priv_c->coups;
  spin_t * const spins = priv_c->spins, *spins_k;
  
  float h_delta;
  
  for (i = 0; i < size; ++i)
    {
      k = rand() % size;
      spins_k = spins + k;
      
      index = k * n_inters;
      h_delta = h_delta_calc (nbors + index, coups + index,
                  spins, spins_k, n_inters);
      
      if (metr_prop_accept (h_delta, beta))
        *spins_k *= -1;
    }
}

void
cnf_mcstep_wolff (void *priv, float const beta)
{
  size_t i, k, head, index_cur, index_nbor;
  ising_priv_t *priv_c = (ising_priv_t *) priv;
  
  size_t const size = priv_c->size;
  size_t const n_inters = priv_c->n_inters;
  spin_t * const spins = priv_c->spins, spin_cur;
  size_t * const nbors = priv_c->nbors, *nbors_cur;
  float * const coups = priv_c->coups, *coups_cur;
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
  
  /*printf ("p_add=%f\n", p_add);
  printf ("k=%ld\n", k);*/
  
  while (head)
    {
      index_cur = buffer[--head];
      spin_cur = spins[index_cur];
      nbors_cur = nbors + index_cur * n_inters;
      coups_cur = coups + index_cur * n_inters;
      /*spins[index_cur] *= -1;*/
      
      /*printf ("index_cur=%ld spin_cur=%d\n", index_cur, spin_cur);*/
      
      for (i = 0; i < n_inters; ++i)
        {
          index_nbor = nbors_cur[i];
          /*printf ("index_nbor=%ld spin_nbor=%d\n", index_nbor, spins[index_nbor]);*/
          if (!in_cluster[index_nbor] && spin_cur == spins[index_nbor])
            {
              p_add = 1.0 - exp (-2.0 * beta * coups_cur[i] * spins[nbors_cur[i]] * spin_cur);
              rand_num = (double) rand () / RAND_MAX;
              /*printf ("rand_num=%f\n", rand_num);*/
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

static const spnr_latt_kind_t cnf_kind =
{
  "icnf",
  &cnf_alloc,
  &cnf_free,
  &print_spins_2d,
  &print_spins_3d,
  &cnf_calc_h,
  &calc_m,
  &cnf_mcstep_metr,
  &cnf_mcstep_wolff,
};

const spnr_latt_kind_t *spnr_ising_cubicnn_ferr = &cnf_kind;

/*static const spnr_latt_kind_t ilrf_kind =
{
	"icnf",
	&ising_priv_ilrf_alloc,
	&ising_priv_ilrf_free,
	&ising_priv_print_spins_2d,
	&ising_priv_print_spins_3d,
	&ising_priv_ilrf_calc_h,
	&ising_priv_calc_m,
	NULL
};

const spnr_latt_kind_t *spnr_ising_longrange_ferr = &ilrf_kind;*/
