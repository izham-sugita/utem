@echo off
setlocal


rem <- �R�����g�A�E�g�� rem
rem ���̕ϐ��͓K�X�C��
set EXE_FILE=fem_solve.exe
set MODEL_FILE=bar.bdf
set RESULT_NAME=bar

del %EXE_FILE%

echo ############# BUILD %EXE_FILE% START #############
nmake -f Makefile.win
echo ############# BUILD %EXE_FILE% END #############

rem ���s�t�@�C��������Ύ��s
IF EXIST %EXE_FILE% (
	rem �r���h�����D���s��
	goto run
) ELSE (
	rem �r���h���s�D�I��
	echo ############# BUILD ERROR %EXE_FILE% #############
	endlocal
	goto :EOF
)

:run
	echo ############# RUN %EXE_FILE% #############
	echo %EXE_FILE% %MODEL_FILE% %RESULT_NAME%
	%EXE_FILE% %MODEL_FILE% %RESULT_NAME%
	echo ############# END %EXE_FILE% #############

endlocal

