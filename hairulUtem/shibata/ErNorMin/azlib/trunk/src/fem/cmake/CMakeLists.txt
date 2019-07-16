#####################################################################
# CMakeLists.txt(fem)
#
# Copyright (C) 2016 AzLib Developers Group
#
# Written by
#  <Takahiro OKADA>
#
# Refer to log documents for details.
#####################################################################

cmake_minimum_required(VERSION 3.0)
project(fem)

include_directories(../../base)
include_directories(../../mathematics)

add_library(${PROJECT_NAME} STATIC
#	../adv_component.c
	../ansys_component.c
	../avs_utl.c
	../contact.c
	../dxf_utl.c
	../eigen_solver.c
	../element_volume.c
	../eps_utl.c
	../extract_surface.c
	../face_edge_table.c
	../fem_geometry.c
	../fem_renumber.c
	../fem_search_common.c
	../fem_solver.c
	../fem_solver_6dof.c
	../fem_solver_plate_shell.c
	../fem_solver_shell.c
	../fem_struct.c
	../fem_struct_binio.c
	../fem_struct_print.c
	../fem_utl.c
	../flow_solver.c
	../integration.c
	../interpolation.c
	../matrix_db.c
	../matrix_db_amls.c
	../nl_solver.c
	../nl_solver_shell.c
	../nst_blk_component.c
	../nst_blk_llio.c
	../nst_f06_component.c
	../nst_pch_component.c
	../optimization_utl.c
	../pressure_force.c
	../sky_cholesky.c
	../sky_crout.c
	../stl_utl.c
	../sysnoise_component.c
)

# target_link_libraries(${PROJECT_NAME} m) # -lm
# target_link_libraries(${PROJECT_NAME} base)
# target_link_libraries(${PROJECT_NAME} mathematics)
