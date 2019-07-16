#!/bin/sh
make clean
make
time ./matrix_db_test /dev/shm
#time ./test_search_matrix_db
