
#include "aiolosmesh.hxx"
#include <interpolation.hxx>

/********************************************************
 * BOUT++ Library - Write fluid simulations in curviilinear geometry
 * Copyright (C) 2016, 2017, 2018 David Schwörer
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
 *****************************************************************/

void AiolosMesh::interp_to_CtoL_Field3D_x(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 2; x < LocalNx - 1; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            -static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
            static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];
      }
    }
  }
  int x = 1;
  for (int y = 0; y < LocalNy; ++y) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x + 2) * LocalNy + y) * LocalNz + z];
      result_ptr[((x - 1) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(21. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 2) * LocalNy + y) * LocalNz + z];
    }
  }
  x = LocalNx - 1;
  for (int y = 0; y < LocalNy; ++y) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x - 3) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z];
    }
  }
}
void AiolosMesh::interp_to_CtoL_Field3D_y(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 2; y < LocalNy - 1; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            -static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
            static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];
      }
    }
  }
  int y = 1;
  for (int x = 0; x < LocalNx; ++x) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
          +static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] +
          static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x)*LocalNy + y + 2) * LocalNz + z];
      result_ptr[((x)*LocalNy + y - 1) * LocalNz + z] =
          +static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
          static_cast<const BoutReal>(21. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 2) * LocalNz + z];
    }
  }
  y = LocalNy - 1;
  for (int x = 0; x < LocalNx; ++x) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
          +static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x)*LocalNy + y - 3) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z];
    }
  }
}
void AiolosMesh::interp_to_CtoL_Field3D_z(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
        }
        {
          int z = 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 2 + LocalNz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
        }
        for (int z = 2; z < LocalNz - 1; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1];
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 2] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz];
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz +
                         +((z - 2 + 2 * LocalNz) % LocalNz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz +
                         +((z - 1 + 1 * LocalNz) % LocalNz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)];
        }
      }
    }
  }
}
void AiolosMesh::interp_to_LtoC_Field3D_x(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (xstart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in X direction!");
  }
#endif
  for (int x = 1; x < LocalNx - 2; ++x) {
    for (int y = 0; y < LocalNy; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
            -static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
            static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x + 2) * LocalNy + y) * LocalNz + z];
      }
    }
  }
  int x = 0;
  for (int y = 0; y < LocalNy; ++y) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 2) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x + 3) * LocalNy + y) * LocalNz + z];
    }
  }
  x = LocalNx - 2;
  for (int y = 0; y < LocalNy; ++y) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x + 0) * LocalNy + y) * LocalNz + z] =
          +static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];
      result_ptr[((x + 1) * LocalNy + y) * LocalNz + z] =
          -static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x - 2) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(21. / 16.) *
              in_ptr[((x - 1) * LocalNy + y) * LocalNz + z] -
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x + 0) * LocalNy + y) * LocalNz + z] +
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x + 1) * LocalNy + y) * LocalNz + z];
    }
  }
}
void AiolosMesh::interp_to_LtoC_Field3D_y(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
#if CHECK > 0
  if (ystart < 2) {
    throw BoutException(
        "Cannot compute derivative - need at least 2 guard cells in Y direction!");
  }
#endif
  for (int x = 0; x < LocalNx; ++x) {
    for (int y = 1; y < LocalNy - 2; ++y) {
      for (int z = 0; z < LocalNz; ++z) {
        result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
            -static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
            static_cast<const BoutReal>(9. / 16.) *
                in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
            static_cast<const BoutReal>(1. / 16.) *
                in_ptr[((x)*LocalNy + y + 2) * LocalNz + z];
      }
    }
  }
  int y = 0;
  for (int x = 0; x < LocalNx; ++x) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
          +static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 2) * LocalNz + z] +
          static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x)*LocalNy + y + 3) * LocalNz + z];
    }
  }
  y = LocalNy - 2;
  for (int x = 0; x < LocalNx; ++x) {
    for (int z = 0; z < LocalNz; ++z) {
      result_ptr[((x)*LocalNy + y + 0) * LocalNz + z] =
          +static_cast<const BoutReal>(1. / 16.) *
              in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] -
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] +
          static_cast<const BoutReal>(15. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
          static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];
      result_ptr[((x)*LocalNy + y + 1) * LocalNz + z] =
          -static_cast<const BoutReal>(5. / 16.) *
              in_ptr[((x)*LocalNy + y - 2) * LocalNz + z] +
          static_cast<const BoutReal>(21. / 16.) *
              in_ptr[((x)*LocalNy + y - 1) * LocalNz + z] -
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x)*LocalNy + y + 0) * LocalNz + z] +
          static_cast<const BoutReal>(35. / 16.) *
              in_ptr[((x)*LocalNy + y + 1) * LocalNz + z];
    }
  }
}
void AiolosMesh::interp_to_LtoC_Field3D_z(BoutReal *__restrict__ result_ptr,
                                          const BoutReal *__restrict__ in_ptr) const {
  if (LocalNz > 3) {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        {
          int z = 0;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1 + LocalNz] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 2];
        }
        for (int z = 1; z < LocalNz - 2; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 2];
        }
        {
          int z = LocalNz - 2;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz];
        }
        {
          int z = LocalNz - 1;

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z - 1] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 0] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 1 - LocalNz] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + z + 2 - LocalNz];
        }
      }
    }
  } else {
    for (int x = 0; x < LocalNx; ++x) {
      for (int y = 0; y < LocalNy; ++y) {
        for (int z = 0; z < LocalNz; ++z) {

          result_ptr[((x)*LocalNy + y) * LocalNz + z + 0] =
              -static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz +
                         +((z - 1 + 1 * LocalNz) % LocalNz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 0) % LocalNz)] +
              static_cast<const BoutReal>(9. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 1) % LocalNz)] -
              static_cast<const BoutReal>(1. / 16.) *
                  in_ptr[((x)*LocalNy + y) * LocalNz + +((z + 2) % LocalNz)];
        }
      }
    }
  }
}

Field3D AiolosMesh::interp_to_do_non_uniform_x_CtoL(const Field3D &in) const {
  Field3D result(const_cast<AiolosMesh *>(this));
  result.allocate();

  if (!stencil_x_CtoL.isSet) {
    const Field2D &dy = coords->dx;
    const_cast<AiolosMesh *>(this)->stencil_x_CtoL =
        Stencil<Field2D>(const_cast<AiolosMesh *>(this));
    const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c1.allocate();
    const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c2.allocate();
    const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c3.allocate();
    const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c4.allocate();
    for (DataIterator di{0, 1 - 1, 0, LocalNy - 1, 0, 0}; !di.done(); di++) {
      auto a = 0.5 * dy[di];
      auto b = dy[di] + 0.5 * dy[di.xp()];
      auto c = dy[di] + dy[di.xp()] + 0.5 * dy[di.xpp()];
      auto d = dy[di] + dy[di.xp()] + dy[di.xp(2)] + 0.5 * dy[di.xp(3)];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c4[di] = sten.c4;
    }
    for (DataIterator di{1, 2 - 1, 0, LocalNy - 1, 0, 0}; !di.done(); di++) {
      auto a = -0.5 * dy[di.xm()];
      auto b = 0.5 * dy[di];
      auto c = dy[di] + 0.5 * dy[di.xp()];
      auto d = dy[di] + dy[di.xp()] + 0.5 * dy[di.xpp()];
      ;
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c4[di] = sten.c4;
    }
    for (DataIterator di{LocalNx - 1, LocalNx + 0 - 1, 0, LocalNy - 1, 0, 0}; !di.done();
         di++) {
      auto a = -0.5 * dy[di.xm(3)] - dy[di.xm(2)] - dy[di.xm()];
      auto b = -0.5 * dy[di.xmm()] - dy[di.xm()];
      auto c = -0.5 * dy[di.xm()];
      auto d = 0.5 * dy[di];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c4[di] = sten.c4;
    }
    for (DataIterator di{2, LocalNx - 1 - 1, 0, LocalNy - 1, 0, 0}; !di.done(); ++di) {
      auto a = -0.5 * dy[di.xmm()] - dy[di.xm()];
      auto b = -0.5 * dy[di.xm()];
      auto c = 0.5 * dy[di];
      auto d = dy[di] + 0.5 * dy[di.xp()];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_x_CtoL.c4[di] = sten.c4;
    }
  }
  for (DataIterator di{0, 2 - 1, 0, LocalNy - 1, 0, LocalNz - 1}; !di.done(); ++di) {
    Indices c1, c2, c3, c4{di.x, di.y, di.z};
    c1 = c2 = c3 = c4;
    c1.x = 0;
    c2.x = 1;
    c3.x = 2;
    c4.x = 3;
    result[di] = +stencil_x_CtoL.c1[di] * in[c1] + stencil_x_CtoL.c2[di] * in[c2] +
                 stencil_x_CtoL.c3[di] * in[c3] + stencil_x_CtoL.c4[di] * in[c4];
  }
  for (DataIterator di{2, LocalNx - 1 - 1, 0, LocalNy - 1, 0, LocalNz - 1}; !di.done();
       ++di) {
    result[di] = +stencil_x_CtoL.c1[di] * in[di.xmm()] +
                 stencil_x_CtoL.c2[di] * in[di.xm()] + stencil_x_CtoL.c3[di] * in[di] +
                 stencil_x_CtoL.c4[di] * in[di.xp()];
  }
  for (DataIterator di{LocalNx - 2, LocalNx + 0 - 1, 0, LocalNy - 1, 0, LocalNz - 1};
       !di.done(); ++di) {
    Indices c1, c2, c3, c4{di.x, di.y, di.z};
    c1 = c2 = c3 = c4;
    c1.x = LocalNx - 4;
    c2.x = LocalNx - 3;
    c3.x = LocalNx - 2;
    c4.x = LocalNx - 1;
    result[di] = +stencil_x_CtoL.c1[di] * in[c1] + stencil_x_CtoL.c2[di] * in[c2] +
                 stencil_x_CtoL.c3[di] * in[c3] + stencil_x_CtoL.c4[di] * in[c4];
  }
  result.setLocation(CELL_XLOW);
  return result;
}

Field3D AiolosMesh::interp_to_do_non_uniform_x_LtoC(const Field3D &in) const {
  Field3D result(const_cast<AiolosMesh *>(this));
  result.allocate();

  if (!stencil_x_LtoC.isSet) {
    const Field2D &dy = coords->dx;
    const_cast<AiolosMesh *>(this)->stencil_x_LtoC =
        Stencil<Field2D>(const_cast<AiolosMesh *>(this));
    const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c1.allocate();
    const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c2.allocate();
    const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c3.allocate();
    const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c4.allocate();
    //////////////////////// L -> C
    for (DataIterator di{0, 1 - 1, 0, LocalNy - 1, 0, 0}; !di.done(); di++) {
      auto b = 0.5 * dy[di];
      auto a = -b;
      auto c = b + dy[di.xp()];
      auto d = c + dy[di.xp(2)];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c4[di] = sten.c4;
    }
    for (DataIterator di{LocalNx - 2, LocalNx - 1 - 1, 0, LocalNy - 1, 0, 0}; !di.done();
         di++) {
      auto d = 0.5 * dy[di];
      auto c = -d;
      auto b = c - dy[di.xm()];
      auto a = b - dy[di.xm(2)];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c4[di] = sten.c4;
    }
    for (DataIterator di{LocalNx - 1, LocalNx + 0 - 1, 0, LocalNy - 1, 0, 0}; !di.done();
         di++) {
      auto d = -0.5 * dy[di];
      auto c = d - dy[di.xm()];
      auto b = c - dy[di.xm(2)];
      auto a = b - dy[di.xm(3)];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c4[di] = sten.c4;
    }
    for (DataIterator di{1, LocalNx - 2 - 1, 0, LocalNy - 1, 0, 0}; !di.done(); ++di) {
      auto c = 0.5 * dy[di];
      auto b = -c;
      auto a = b - dy[di.xm()];
      auto d = c + dy[di.xp()];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_x_LtoC.c4[di] = sten.c4;
    }
  }
  for (DataIterator di{0, 2 - 1, 0, LocalNy - 1, 0, LocalNz - 1}; !di.done(); ++di) {
    Indices c1, c2, c3, c4{di.x, di.y, di.z};
    c1 = c2 = c3 = c4;
    c1.x = 0;
    c2.x = 1;
    c3.x = 2;
    c4.x = 3;
    result[di] = +stencil_x_LtoC.c1[di] * in[c1] + stencil_x_LtoC.c2[di] * in[c2] +
                 stencil_x_LtoC.c3[di] * in[c3] + stencil_x_LtoC.c4[di] * in[c4];
  }
  for (DataIterator di{1, LocalNx - 2 - 1, 0, LocalNy - 1, 0, LocalNz - 1}; !di.done();
       ++di) {
    result[di] = +stencil_x_LtoC.c1[di] * in[di.xm()] + stencil_x_LtoC.c2[di] * in[di] +
                 stencil_x_LtoC.c3[di] * in[di.xp()] +
                 stencil_x_LtoC.c4[di] * in[di.xpp()];
  }
  for (DataIterator di{LocalNx - 2, LocalNx + 0 - 1, 0, LocalNy - 1, 0, LocalNz - 1};
       !di.done(); ++di) {
    Indices c1, c2, c3, c4{di.x, di.y, di.z};
    c1 = c2 = c3 = c4;
    c1.x = LocalNx - 4;
    c2.x = LocalNx - 3;
    c3.x = LocalNx - 2;
    c4.x = LocalNx - 1;
    result[di] = +stencil_x_LtoC.c1[di] * in[c1] + stencil_x_LtoC.c2[di] * in[c2] +
                 stencil_x_LtoC.c3[di] * in[c3] + stencil_x_LtoC.c4[di] * in[c4];
  }
  result.setLocation(CELL_CENTRE);
  return result;
}

Field3D AiolosMesh::interp_to_do_non_uniform_y_CtoL(const Field3D &in) const {
  Field3D result(const_cast<AiolosMesh *>(this));
  result.allocate();

  if (!stencil_y_CtoL.isSet) {
    const Field2D &dy = coords->dy;
    const_cast<AiolosMesh *>(this)->stencil_y_CtoL =
        Stencil<Field2D>(const_cast<AiolosMesh *>(this));
    const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c1.allocate();
    const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c2.allocate();
    const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c3.allocate();
    const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c4.allocate();
    for (DataIterator di{0, LocalNx - 1, 0, 1 - 1, 0, 0}; !di.done(); di++) {
      auto a = 0.5 * dy[di];
      auto b = dy[di] + 0.5 * dy[di.yp()];
      auto c = dy[di] + dy[di.yp()] + 0.5 * dy[di.ypp()];
      auto d = dy[di] + dy[di.yp()] + dy[di.yp(2)] + 0.5 * dy[di.yp(3)];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c4[di] = sten.c4;
    }
    for (DataIterator di{0, LocalNx - 1, 1, 2 - 1, 0, 0}; !di.done(); di++) {
      auto a = -0.5 * dy[di.ym()];
      auto b = 0.5 * dy[di];
      auto c = dy[di] + 0.5 * dy[di.yp()];
      auto d = dy[di] + dy[di.yp()] + 0.5 * dy[di.ypp()];
      ;
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c4[di] = sten.c4;
    }
    for (DataIterator di{0, LocalNx - 1, LocalNy - 1, LocalNy + 0 - 1, 0, 0}; !di.done();
         di++) {
      auto a = -0.5 * dy[di.ym(3)] - dy[di.ym(2)] - dy[di.ym()];
      auto b = -0.5 * dy[di.ymm()] - dy[di.ym()];
      auto c = -0.5 * dy[di.ym()];
      auto d = 0.5 * dy[di];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c4[di] = sten.c4;
    }
    for (DataIterator di{0, LocalNx - 1, 2, LocalNy - 1 - 1, 0, 0}; !di.done(); ++di) {
      auto a = -0.5 * dy[di.ymm()] - dy[di.ym()];
      auto b = -0.5 * dy[di.ym()];
      auto c = 0.5 * dy[di];
      auto d = dy[di] + 0.5 * dy[di.yp()];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_y_CtoL.c4[di] = sten.c4;
    }
  }
  for (DataIterator di{0, LocalNx - 1, 0, 2 - 1, 0, LocalNz - 1}; !di.done(); ++di) {
    Indices c1, c2, c3, c4{di.x, di.y, di.z};
    c1 = c2 = c3 = c4;
    c1.y = 0;
    c2.y = 1;
    c3.y = 2;
    c4.y = 3;
    result[di] = +stencil_y_CtoL.c1[di] * in[c1] + stencil_y_CtoL.c2[di] * in[c2] +
                 stencil_y_CtoL.c3[di] * in[c3] + stencil_y_CtoL.c4[di] * in[c4];
  }
  for (DataIterator di{0, LocalNx - 1, 2, LocalNy - 1 - 1, 0, LocalNz - 1}; !di.done();
       ++di) {
    result[di] = +stencil_y_CtoL.c1[di] * in[di.ymm()] +
                 stencil_y_CtoL.c2[di] * in[di.ym()] + stencil_y_CtoL.c3[di] * in[di] +
                 stencil_y_CtoL.c4[di] * in[di.yp()];
  }
  for (DataIterator di{0, LocalNx - 1, LocalNy - 2, LocalNy + 0 - 1, 0, LocalNz - 1};
       !di.done(); ++di) {
    Indices c1, c2, c3, c4{di.x, di.y, di.z};
    c1 = c2 = c3 = c4;
    c1.y = LocalNy - 4;
    c2.y = LocalNy - 3;
    c3.y = LocalNy - 2;
    c4.y = LocalNy - 1;
    result[di] = +stencil_y_CtoL.c1[di] * in[c1] + stencil_y_CtoL.c2[di] * in[c2] +
                 stencil_y_CtoL.c3[di] * in[c3] + stencil_y_CtoL.c4[di] * in[c4];
  }
  result.setLocation(CELL_YLOW);
  return result;
}

Field3D AiolosMesh::interp_to_do_non_uniform_y_LtoC(const Field3D &in) const {
  Field3D result(const_cast<AiolosMesh *>(this));
  result.allocate();

  if (!stencil_y_LtoC.isSet) {
    const Field2D &dy = coords->dy;
    const_cast<AiolosMesh *>(this)->stencil_y_LtoC =
        Stencil<Field2D>(const_cast<AiolosMesh *>(this));
    const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c1.allocate();
    const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c2.allocate();
    const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c3.allocate();
    const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c4.allocate();
    //////////////////////// L -> C
    for (DataIterator di{0, LocalNx - 1, 0, 1 - 1, 0, 0}; !di.done(); di++) {
      auto b = 0.5 * dy[di];
      auto a = -b;
      auto c = b + dy[di.yp()];
      auto d = c + dy[di.yp(2)];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c4[di] = sten.c4;
    }
    for (DataIterator di{0, LocalNx - 1, LocalNy - 2, LocalNy - 1 - 1, 0, 0}; !di.done();
         di++) {
      auto d = 0.5 * dy[di];
      auto c = -d;
      auto b = c - dy[di.ym()];
      auto a = b - dy[di.ym(2)];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c4[di] = sten.c4;
    }
    for (DataIterator di{0, LocalNx - 1, LocalNy - 1, LocalNy + 0 - 1, 0, 0}; !di.done();
         di++) {
      auto d = -0.5 * dy[di];
      auto c = d - dy[di.ym()];
      auto b = c - dy[di.ym(2)];
      auto a = b - dy[di.ym(3)];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c4[di] = sten.c4;
    }
    for (DataIterator di{0, LocalNx - 1, 1, LocalNy - 2 - 1, 0, 0}; !di.done(); ++di) {
      auto c = 0.5 * dy[di];
      auto b = -c;
      auto a = b - dy[di.ym()];
      auto d = c + dy[di.yp()];
      auto sten = calc_interp_to_stencil(a, b, c, d);
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c1[di] = sten.c1;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c2[di] = sten.c2;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c3[di] = sten.c3;
      const_cast<AiolosMesh *>(this)->stencil_y_LtoC.c4[di] = sten.c4;
    }
  }
  for (DataIterator di{0, LocalNx - 1, 0, 2 - 1, 0, LocalNz - 1}; !di.done(); ++di) {
    Indices c1, c2, c3, c4{di.x, di.y, di.z};
    c1 = c2 = c3 = c4;
    c1.y = 0;
    c2.y = 1;
    c3.y = 2;
    c4.y = 3;
    result[di] = +stencil_y_LtoC.c1[di] * in[c1] + stencil_y_LtoC.c2[di] * in[c2] +
                 stencil_y_LtoC.c3[di] * in[c3] + stencil_y_LtoC.c4[di] * in[c4];
  }
  for (DataIterator di{0, LocalNx - 1, 1, LocalNy - 2 - 1, 0, LocalNz - 1}; !di.done();
       ++di) {
    result[di] = +stencil_y_LtoC.c1[di] * in[di.ym()] + stencil_y_LtoC.c2[di] * in[di] +
                 stencil_y_LtoC.c3[di] * in[di.yp()] +
                 stencil_y_LtoC.c4[di] * in[di.ypp()];
  }
  for (DataIterator di{0, LocalNx - 1, LocalNy - 2, LocalNy + 0 - 1, 0, LocalNz - 1};
       !di.done(); ++di) {
    Indices c1, c2, c3, c4{di.x, di.y, di.z};
    c1 = c2 = c3 = c4;
    c1.y = LocalNy - 4;
    c2.y = LocalNy - 3;
    c3.y = LocalNy - 2;
    c4.y = LocalNy - 1;
    result[di] = +stencil_y_LtoC.c1[di] * in[c1] + stencil_y_LtoC.c2[di] * in[c2] +
                 stencil_y_LtoC.c3[di] * in[c3] + stencil_y_LtoC.c4[di] * in[c4];
  }
  result.setLocation(CELL_CENTRE);
  return result;
}

const Field3D AiolosMesh::interp_to_do(const Field3D &f, CELL_LOC loc,
                                       REGION region) const {
  Mesh *localmesh = const_cast<AiolosMesh *>(this);
  CELL_LOC dir = f.getLocation() == CELL_CENTRE ? loc : f.getLocation();
  switch (dir) {
  case CELL_XLOW:
    if (is_x_uniform == 2) {
      const_cast<AiolosMesh *>(this)->is_x_uniform =
          min(localmesh->coordinates()->dx, false, RGN_ALL) ==
          max(localmesh->coordinates()->dx, false, RGN_ALL);
    }
    if (is_x_uniform) {
      Field3D result(localmesh);
      result.allocate();
      Indices i0{0, 0, 0};
      if (f.getLocation() != CELL_CENTRE) {
        // we first need to go back to centre before we can go anywhere else
        interp_to_LtoC_Field3D_x(&result[i0], &f[i0]);
        result.setLocation(CELL_CENTRE);
        return interp_to(result, loc, region);
      } else {
        interp_to_CtoL_Field3D_x(&result[i0], &f[i0]);
        result.setLocation(CELL_XLOW);
        return result;
      }
    } else {
      if (f.getLocation() != CELL_CENTRE) {
        return interp_to(interp_to_do_non_uniform_x_LtoC(f), loc, region);
      } else {
        return interp_to_do_non_uniform_x_CtoL(f);
      }
    }
    break;
  case CELL_YLOW:
    if (is_y_uniform == 2) {
      const_cast<AiolosMesh *>(this)->is_y_uniform =
          min(localmesh->coordinates()->dy, false, RGN_ALL) ==
          max(localmesh->coordinates()->dy, false, RGN_ALL);
    }
    if (is_y_uniform) {
      Field3D result(localmesh);
      result.allocate();
      Indices i0{0, 0, 0};
      if (f.getLocation() != CELL_CENTRE) {
        // we first need to go back to centre before we can go anywhere else
        interp_to_LtoC_Field3D_y(&result[i0], &f[i0]);
        result.setLocation(CELL_CENTRE);
        return interp_to(result, loc, region);
      } else {
        interp_to_CtoL_Field3D_y(&result[i0], &f[i0]);
        result.setLocation(CELL_YLOW);
        return result;
      }
    } else {
      if (f.getLocation() != CELL_CENTRE) {
        return interp_to(interp_to_do_non_uniform_y_LtoC(f), loc, region);
      } else {
        return interp_to_do_non_uniform_y_CtoL(f);
      }
    }
    break;
  case CELL_ZLOW:
    if (is_z_uniform) {
      Field3D result(localmesh);
      result.allocate();
      Indices i0{0, 0, 0};
      if (f.getLocation() != CELL_CENTRE) {
        // we first need to go back to centre before we can go anywhere else
        interp_to_LtoC_Field3D_z(&result[i0], &f[i0]);
        result.setLocation(CELL_CENTRE);
        return interp_to(result, loc, region);
      } else {
        interp_to_CtoL_Field3D_z(&result[i0], &f[i0]);
        result.setLocation(CELL_ZLOW);
        return result;
      }
    }
    break;
  default:
    throw BoutException("AiolosMesh::interp_to: Cannot interpolate to %s!",
                        strLocation(loc));
  }
  // gcc thinks we might end up here. So be save and throw.
  throw BoutException("AiolosMesh::interp_to - failed to return result!");
}
