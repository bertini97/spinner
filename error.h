/* error.h
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

#ifndef ERROR_H
#define ERROR_H

#undef BEGIN_C_DECLS
#undef END_C_DECLS
#ifdef __cplusplus
# define BEGIN_C_DECLS extern "C" {
# define END_C_DECLS }
#else
# define BEGIN_C_DECLS /* empty */
# define END_C_DECLS /* empty */
#endif

BEGIN_C_DECLS

enum
{
  SPNR_SUCCESS          = 0, /* success */
  SPNR_FAILURE          = 1, /* generic failure */
  SPNR_ERROR_PARAM_OOB  = 2, /* given parameters out of bounds */
  SPNR_ERROR_ALLOC      = 3, /* memory allocation failure */
  SPNR_ERROR_FUNC_NULL  = 4, /* func is NULL for that lattice kind*/
  SPNR_ERROR_ARG_NULL   = 5  /* function argument is NULL */
};

extern void spnr_warn (int warn, char const *mess);
extern void spnr_err (int err, char const *mess);

extern void * malloc_err (size_t size);

END_C_DECLS

#endif
