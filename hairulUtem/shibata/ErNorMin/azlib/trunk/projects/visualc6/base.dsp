# Microsoft Developer Studio Project File - Name="base" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** 編集しないでください **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=base - Win32 Debug
!MESSAGE これは有効なﾒｲｸﾌｧｲﾙではありません。 このﾌﾟﾛｼﾞｪｸﾄをﾋﾞﾙﾄﾞするためには NMAKE を使用してください。
!MESSAGE [ﾒｲｸﾌｧｲﾙのｴｸｽﾎﾟｰﾄ] ｺﾏﾝﾄﾞを使用して実行してください
!MESSAGE 
!MESSAGE NMAKE /f "base.mak".
!MESSAGE 
!MESSAGE NMAKE の実行時に構成を指定できます
!MESSAGE ｺﾏﾝﾄﾞ ﾗｲﾝ上でﾏｸﾛの設定を定義します。例:
!MESSAGE 
!MESSAGE NMAKE /f "base.mak" CFG="base - Win32 Debug"
!MESSAGE 
!MESSAGE 選択可能なﾋﾞﾙﾄﾞ ﾓｰﾄﾞ:
!MESSAGE 
!MESSAGE "base - Win32 Release" ("Win32 (x86) Static Library" 用)
!MESSAGE "base - Win32 Debug" ("Win32 (x86) Static Library" 用)
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "base - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x411 /d "NDEBUG"
# ADD RSC /l 0x411 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "base - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x411 /d "_DEBUG"
# ADD RSC /l 0x411 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"base.lib"

!ENDIF 

# Begin Target

# Name "base - Win32 Release"
# Name "base - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\base\bin_io.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\data_model.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\idc.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\list.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\log_printf.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\mail_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\memory_manager.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\rc.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\scratch_io.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\string_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\time_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\base\wrapper_64.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\src\base\base.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\bin_io.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\data_model.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\list.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\log_printf.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\macros.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\mail_utl.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\memory_manager.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\messages.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\rc.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\scratch_io.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\string_utl.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\time_utl.h
# End Source File
# Begin Source File

SOURCE=..\..\src\base\wrapper_64.h
# End Source File
# End Group
# End Target
# End Project
