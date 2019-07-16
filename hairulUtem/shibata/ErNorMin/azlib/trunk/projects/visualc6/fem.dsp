# Microsoft Developer Studio Project File - Name="fem" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** 編集しないでください **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=fem - Win32 Debug
!MESSAGE これは有効なﾒｲｸﾌｧｲﾙではありません。 このﾌﾟﾛｼﾞｪｸﾄをﾋﾞﾙﾄﾞするためには NMAKE を使用してください。
!MESSAGE [ﾒｲｸﾌｧｲﾙのｴｸｽﾎﾟｰﾄ] ｺﾏﾝﾄﾞを使用して実行してください
!MESSAGE 
!MESSAGE NMAKE /f "fem.mak".
!MESSAGE 
!MESSAGE NMAKE の実行時に構成を指定できます
!MESSAGE ｺﾏﾝﾄﾞ ﾗｲﾝ上でﾏｸﾛの設定を定義します。例:
!MESSAGE 
!MESSAGE NMAKE /f "fem.mak" CFG="fem - Win32 Debug"
!MESSAGE 
!MESSAGE 選択可能なﾋﾞﾙﾄﾞ ﾓｰﾄﾞ:
!MESSAGE 
!MESSAGE "fem - Win32 Release" ("Win32 (x86) Static Library" 用)
!MESSAGE "fem - Win32 Debug" ("Win32 (x86) Static Library" 用)
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "fem - Win32 Release"

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

!ELSEIF  "$(CFG)" == "fem - Win32 Debug"

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
# ADD RSC /l 0x411 /i "..\..\src\base" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"fem.lib"

!ENDIF 

# Begin Target

# Name "fem - Win32 Release"
# Name "fem - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\src\fem\ansys_component.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\avs_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\contact.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\dxf_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\eigen_solver.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\element_volume.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\eps_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\extract_surface.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\face_edge_table.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_geometry.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_renumber.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_search_common.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_solver.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_solver_6dof.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_struct.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_struct_binio.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_struct_print.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\integration.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\interpolation.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\nl_solver.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\nst_blk_component.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\nst_blk_llio.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\nst_f06_component.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\nst_pch_component.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\optimization_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\pressure_force.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\sky_cholesky.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\sky_crout.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\stl_utl.c
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\sysnoise_component.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\src\fem\ansys_component.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\avs_utl.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\dxf_utl.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\eps_utl.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_solver.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\fem_struct.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\nst_component.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\optimization_utl.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\sky_cholesky.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\sky_crout.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\stl_utl.h
# End Source File
# Begin Source File

SOURCE=..\..\src\fem\sysnoise_component.h
# End Source File
# End Group
# End Target
# End Project
