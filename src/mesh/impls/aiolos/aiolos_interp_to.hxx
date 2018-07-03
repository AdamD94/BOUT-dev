
// to be included from aiolosmesh.hxx

void interp_to_CtoL_Field3D_x(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_CtoL_Field3D_y(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_CtoL_Field3D_z(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_LtoC_Field3D_x(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_LtoC_Field3D_y(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
void interp_to_LtoC_Field3D_z(BoutReal *__restrict__ result_ptr,
                              const BoutReal *__restrict__ in_ptr) const;
Field3D interp_to_do_non_uniform_x_CtoL(const Field3D &in) const;
Field3D interp_to_do_non_uniform_x_LtoC(const Field3D &in) const;
Field3D interp_to_do_non_uniform_y_CtoL(const Field3D &in) const;
Field3D interp_to_do_non_uniform_y_LtoC(const Field3D &in) const;
