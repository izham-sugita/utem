#####################################################################
# CMakeLists.txt(mathematics)
#
# Copyright (C) 2016 AzLib Developers Group
#
# Written by
#  <Takahiro OKADA>
#
# Refer to log documents for details.
#####################################################################

cmake_minimum_required(VERSION 3.0)
project(mathematics)

include_directories(../../base)

add_library(${PROJECT_NAME} STATIC
	../complex_utl.c
	../crs_complex_utl.c
	../geometry.c
	../lq.c
	../lu.c
	../math_utl.c
	../matrix.c
	../nonzero_cg.c
	../operation.c
	../pokecom2.c
	../qr.c
	../random.c
	../vect3d.c
	../vect3d_complex.c
	../vect_storage.c
)

target_link_libraries(${PROJECT_NAME} m) # -lm
target_link_libraries(${PROJECT_NAME} base)
