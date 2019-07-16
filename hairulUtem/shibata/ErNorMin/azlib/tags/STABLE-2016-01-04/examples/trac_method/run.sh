#!/bin/sh

EXE_FILE=trac_method
PARAM_FILE=PARAM_bodyless
MODEL_FILE=bar.bdf
REST_FILE=bar_rest.bdf
RESULT_NAME=bar_bodyless

echo ############# BUILD %EXE_FILE% START #############
make
echo ############# BUILD %EXE_FILE% END #############


if [ -x $EXE_FILE ]
then
	echo "############# RUN $EXE_FILE #############"
	echo $EXE_FILE $PARAM_FILE $MODEL_FILE $REST_FILE $RESULT_NAME
	./$EXE_FILE $PARAM_FILE $MODEL_FILE $REST_FILE $RESULT_NAME
	echo "############# END $EXE_FILE #############"
fi
