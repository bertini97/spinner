/* couplings.h
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

#ifndef COUPLINGS_H
#define COUPLINGS_H

float
coups_get_ferr ()
{
  return +SPNR_J;
}


float
coups_get_antiferr ()
{
  return -SPNR_J;
}


float
coups_get_bim ()
{
  float const randnum = (float) rand () / RAND_MAX;
  return (randnum < 0.5) ? (+SPNR_J) : (-SPNR_J);
}

#endif
