#!/bin/sh

#EXE_FILE="2geom_contact"
EXE_FILE="node_to_elem"

make clean
make

if [ -x $EXE_FILE ]
then
	echo "############# RUN $EXE_FILE #############"
	./$EXE_FILE Al_Steel_target_2.NAS Al_Steel_hit_2.NAS
#	./$EXE_FILE knuckle.bdf knuckle.bdf
	echo "############# END $EXE_FILE #############"
fi
