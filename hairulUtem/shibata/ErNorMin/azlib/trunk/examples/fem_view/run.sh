#!/bin/sh

EXE_FILE=main
MODEL_FILE=bar.bdf

echo "############# BUILD $EXE_FILE START #############"
make
echo "############# BUILD $EXE_FILE END   #############"


if [ -x $EXE_FILE ]
then
	echo "############# RUN $EXE_FILE START #############"
	echo $EXE_FILE $MODEL_FILE
	./$EXE_FILE $MODEL_FILE
	echo "############# END $EXE_FILE END   #############"
fi
