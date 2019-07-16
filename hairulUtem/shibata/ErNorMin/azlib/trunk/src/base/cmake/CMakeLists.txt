#####################################################################
# CMakeLists.txt(base)
#
# Copyright (C) 2016 AzLib Developers Group
#
# Written by
#  <Takahiro OKADA>
#
# Refer to log documents for details.
#####################################################################

cmake_minimum_required(VERSION 3.0)
project(base)

add_library(${PROJECT_NAME} STATIC
	../bin_io.c
	../data_model.c
	../idc.c
	../list.c
	../log_printf.c
	../mail_utl.c
	../memory_manager.c
	../rc.c
	../scratch_io.c
	../string_utl.c
	../time_utl.c
	../wrapper_64.c
)
