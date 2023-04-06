/* macros.h
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

#ifndef MACROS_H
#define MACROS_H

#define CONCAT2x(A, B)      A ## _ ## B
#define CONCAT2(A, B)       CONCAT2x(A, B)

#define CONCAT3x(A, B, C)   A ## _ ## B ## _ ## C
#define CONCAT3(A, B, C)    CONCAT3x(A, B, C)

#define FUNC(NAME, SUFFIX)  CONCAT2(NAME, SUFFIX)

#endif
