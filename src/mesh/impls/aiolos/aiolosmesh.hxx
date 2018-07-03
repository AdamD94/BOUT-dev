/**************************************************************************
 * Implementation of the Mesh class. Based on the BoutMesh, but with
 * hopefuly faster derivatives.
 *
 * Changelog
 * ---------
 *
 * 2016..2018 David Schwörer
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

#pragma once

#include "../bout/boutmesh.hxx"

#include <stencils.hxx>

#include <cmath>
#include <list>
#include <vector>

using std::list;
using std::vector;

class AiolosMesh : public BoutMesh {
public:
  AiolosMesh(GridDataSource *s, Options *options = NULL);
  ~AiolosMesh();

  struct cart_diff_lookup_table {
    Mesh::deriv_func func; // Single-argument differencing function
    deriv_func norm;
    deriv_func on;
    deriv_func off;
  };

#include "aiolos_header.hxx"

  // virtual const Field3D interp_to(const Field3D &var, CELL_LOC loc) const;

  virtual const Field3D interp_to(const Field3D &f , CELL_LOC loc, REGION region) override {
    return ((const AiolosMesh*)this)->interp_to(f, loc, region);
  }
  virtual const Field3D interp_to(const Field3D &f , CELL_LOC loc, REGION region) const {
    ASSERT2(f.getMesh() == this);
    if (loc == f.getLocation() || loc == CELL_DEFAULT) {
      return f;
    } else {
      return interp_to_do(f, loc, region);
    }
  }


  virtual void derivs_init(Options *option) override;

// to check in debugger we have the right mesh
#if CHECK > 1
  bool isAiolos = true;
#endif
private:
  int is_x_uniform;
  int is_y_uniform;
  int is_z_uniform;
  const Field3D interp_to_do(const Field3D &f, CELL_LOC loc, REGION region) const;

  template <typename T>
  struct Stencil {
    T c1, c2, c3, c4;
    bool isSet;
    Stencil() : isSet(false) {};
    Stencil(T c1, T c2, T c3, T c4) : c1(c1), c2(c2), c3(c3), c4(c4), isSet(true) {};
    template<typename ... Args>
    Stencil(Args && ... args) : c1(args ...), c2(args ...), c3(args ...), c4(args ...), isSet(false) {};
  };

  Stencil<BoutReal> calc_interp_to_stencil(BoutReal a, BoutReal b, BoutReal c,
                                       BoutReal d) const {
    auto ab = a * b;
    auto cd = c * d;
    auto bma = b - a;
    auto cma = c - a;
    auto cmb = c - b;
    auto dma = d - a;
    auto bmd = b - d;
    auto dmc = d - c;
    auto s1 = bma * cma * dma;
    auto s2 = bma * cmb * bmd;
    auto s3 = cma * cmb * dmc;
    auto s4 = dma * bmd * dmc;
    auto e = cd * b / s1;
    auto f = a * cd / s2;
    auto g = ab * d / s3;
    auto h = ab * c / s4;
    return {e,f,g,h};
  }

  Stencil<Field2D> stencil_x_CtoL, stencil_x_LtoC, stencil_y_CtoL, stencil_y_LtoC;

#include "aiolos_derivs.hxx"

#include "aiolos_stencils.hxx"

#include "aiolos_interp_to.hxx"

};
