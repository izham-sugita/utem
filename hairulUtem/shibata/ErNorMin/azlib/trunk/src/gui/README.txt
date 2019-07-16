# azlib/gui

 $Id: README.txt 1066 2016-05-19 09:53:33Z hayashi $

--------------------------------------------------------------------------------

 * AzLibのGUIを OpenGL および GTK+ で実装している．
 * GTK+上でOpenGLを描画するためのライブラリが GtkGLExt である．

 * gtku_*() 系関数は，このライブラリで独自に定義されているものである．


OpenGL
------
 * glut_utl.c/h: OpenGL用．下記 GTK+ では一切使われていない模様．


GTK+
----
 * gtk_control_panel.c : コントロールパネル
 * gtk_icons.h         : アイコン画像
 * gtk_model_view.c    : モデルビューワーの基礎となる関数
 * gtk_toolbar.c       : ツールバー
 * gtk_utl.c/h


GtkGLExt
--------
 * gtkgl_article.c
 * gtkgl_fem.c
 * gtkgl_image.c
 * gtkgl_utl.c


