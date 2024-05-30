#include "spinner.h"

spnr_sys_t *
spnr_sys_alloc (spnr_graph_t * graph, spnr_sys_kind_t const * kind, size_t param)
{
  spnr_sys_t * sys = malloc (sizeof (spnr_sys_t));
  sys->graph = graph;
  sys->kind = kind;
  sys->priv = kind->priv_alloc (graph->N, param);
  
  return sys;
}

void
spnr_sys_free (spnr_sys_t * sys)
{
  sys->kind->priv_free(sys->priv);
  free (sys);
}

float
spnr_sys_spin_size (spnr_sys_t * sys)
{
  return sys->kind->spin_size(sys->priv);
}

float spnr_sys_calc_h (spnr_sys_t * sys)
{
  spnr_graph_t *g = sys->graph;
  return g->kind->calc_h (g->priv, g->N, sys);
}

float spnr_sys_calc_phi (spnr_sys_t * sys)
{
  return sys->kind->calc_phi (sys->priv, sys->graph->N);
}