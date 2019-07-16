/*********************************************************************
 * gtk_graph.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: gtk_graph.c,v 1.2 2003/07/22 08:02:11 sasaoka Exp $ */

#include <stdio.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include "gtk_utl.h"

GtkWidget *gtku_drow_graph()
{	
	GtkWidget *vbox;
	GtkWidget *frame;


	/* Create main Vbox */
	vbox = gtk_vbox_new(FALSE, 2);
	gtk_widget_show(vbox);

	/* Create new FRAME */
	frame = gtk_frame_new("GRAPH AREA");
	gtk_frame_set_label_align(GTK_FRAME(frame), 0.5, 0);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_OUT);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(frame), TRUE, TRUE, 0);
	gtk_widget_show(frame);

	return(vbox);
}

