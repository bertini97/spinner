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

#include "spinner.h"
#include "error.h"

spnr_graph_t *
spnr_graph_alloc (spnr_graph_kind_t const * const kind,
                  float const (*getter)(),
                  size_t const N,
                  size_t const param)
{
  spnr_graph_t *graph = malloc_err (sizeof(spnr_graph_t));
  graph->N = N;
  graph->kind = kind;
  graph->priv = kind->priv_alloc(getter, N, param);
  return graph;
}

void
spnr_graph_free (spnr_graph_t * const graph)
{
  graph->kind->priv_free (graph->priv);
  free (graph);
}