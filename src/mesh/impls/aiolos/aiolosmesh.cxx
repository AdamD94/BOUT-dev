/**************************************************************************
 * Implementation of the Mesh class, handling input files compatible with
 * BOUT++.
 *
 * Changelog
 * ---------
 *
 * 2016-09 David Schwörer
 *           based on the BoutMesh
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include "aiolosmesh.hxx"

#include <output.hxx>

AiolosMesh::AiolosMesh(GridDataSource *s, Options *options)
    : BoutMesh(s, options), is_x_uniform(2), is_y_uniform(2), is_z_uniform(1),
      stencil_x_CtoL(this), stencil_x_LtoC(this), stencil_y_CtoL(this),
      stencil_y_LtoC(this) {
  output_info.write("  Using Aiolos Mesh!\n");
  derivs_init(options);
}
AiolosMesh::~AiolosMesh() {}
