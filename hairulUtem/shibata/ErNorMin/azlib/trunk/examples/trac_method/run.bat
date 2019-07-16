@echo off
setlocal


rem <- コメントアウトは rem
rem この変数は適宜修正
set EXE_FILE=trac_method.exe
set PARAM_FILE=PARAM_bodyless
set MODEL_FILE=bar.bdf
set REST_FILE=bar_rest.bdf
set RESULT_NAME=bar_bodyless

del %EXE_FILE%

echo ############# BUILD %EXE_FILE% START #############
nmake -f Makefile.win
echo ############# BUILD %EXE_FILE% END #############

rem 実行ファイルがあれば実行
IF EXIST %EXE_FILE% (
	rem ビルド成功．実行へ
	goto run
) ELSE (
	rem ビルド失敗．終了
	echo ############# BUILD ERROR %EXE_FILE% #############
	endlocal
	goto :EOF
)

:run
	echo ############# RUN %EXE_FILE% #############
	echo %EXE_FILE% %PARAM_FILE% %MODEL_FILE% %REST_FILE% %RESULT_NAME%
	%EXE_FILE% %PARAM_FILE% %MODEL_FILE% %REST_FILE% %RESULT_NAME%
	echo ############# END %EXE_FILE% #############

endlocal

