#####################################################################
# CMakeLists.txt(image)
#
# Copyright (C) 2016 AzLib Developers Group
#
# Written by
#  <Takahiro OKADA>
#
# Refer to log documents for details.
#####################################################################

cmake_minimum_required(VERSION 3.0)
project(image)

include_directories(../../base)
include_directories(../../fem)
include_directories(../../mathematics)

add_library(${PROJECT_NAME} STATIC
	../avs_field_utl.c
	../bmp_utl.c
	../dicom_utl.c
	../graphic_binio.c
	../graphic_utl.c
	../tiff_utl.c
	../voxel_utl.c
)

# target_link_libraries(${PROJECT_NAME} m) # -lm
# target_link_libraries(${PROJECT_NAME} base)
# target_link_libraries(${PROJECT_NAME} fem)
# target_link_libraries(${PROJECT_NAME} mathematics)
