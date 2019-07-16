#####################################################################
# CMakeLists.txt(gui)
#
# Copyright (C) 2016 AzLib Developers Group
#
# Written by
#  <Takahiro OKADA>
#
# Refer to log documents for details.
#####################################################################

cmake_minimum_required(VERSION 3.0)
project(gui)

find_package(PkgConfig)

pkg_check_modules(GTK2 REQUIRED gtk+-2.0)
pkg_check_modules(GTKGLEXT REQUIRED gtkglext-1.0)

include_directories(${GTK2_INCLUDE_DIRS} ${GTKGLEXT_INCLUDE_DIRS})
include_directories(../../base)
include_directories(../../fem)
include_directories(../../image)
include_directories(../../mathematics)

add_library(${PROJECT_NAME} STATIC
	../glut_utl.c
	../gtk_control_panel.c
	../gtk_model_view.c
	../gtk_toolbar.c
	../gtk_utl.c
	../gtkgl_article.c
	../gtkgl_fem.c
	../gtkgl_image.c
	../gtkgl_utl.c
)

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

target_link_libraries(${PROJECT_NAME} ${GTK2_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${GTKGLEXT_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${GLUT_LIBRARIES})

# target_link_libraries(${PROJECT_NAME} m) # -lm
# target_link_libraries(${PROJECT_NAME} base)
# target_link_libraries(${PROJECT_NAME} fem)
# target_link_libraries(${PROJECT_NAME} image)
# target_link_libraries(${PROJECT_NAME} mathematics)
