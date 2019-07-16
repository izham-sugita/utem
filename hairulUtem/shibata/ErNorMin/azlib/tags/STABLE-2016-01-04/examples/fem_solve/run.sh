#!/bin/sh

EXE_FILE=fem_solve
MODEL_FILE=bar.bdf
RESULT_NAME=result

echo ############# BUILD %EXE_FILE% START #############
make
echo ############# BUILD %EXE_FILE% END #############


if [ -x $EXE_FILE ]
then
	echo "############# RUN $EXE_FILE #############"
	echo $EXE_FILE $MODEL_FILE $RESULT_NAME
	./$EXE_FILE $MODEL_FILE $RESULT_NAME
	echo "############# END $EXE_FILE #############"
fi
