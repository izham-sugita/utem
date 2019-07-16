/*********************************************************************
 * gtk_utl.c
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

/* $Id: gtk_utl.c,v 1.7 2003/12/04 07:26:08 sasaoka Exp $ */

#include <stdio.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include "gtk_utl.h"


gboolean gtku_window_destory(GtkWidget *widget)
{
	gtk_main_quit();
	return(FALSE);
}


static void dialog_yes (GtkWidget *widget, gint *flag)
{
	*flag = FALSE;
}


static void dialog_no (GtkWidget *widget, gint *flag)
{
	*flag = TRUE;
}


gboolean gtku_window_destory_ask(GtkWidget *widget, GdkEvent *event,
                                 GtkWidget *parent)
{
	GtkWidget *dialog_window = NULL;
	gint flag = TRUE;
	GtkWidget *label;
	GtkWidget *button_y;
	GtkWidget *button_n;
	int dialog_size_x = 250, dialog_size_y = 110;
	gint x, y, w, h;

	if(dialog_window)
		return(FALSE);

	/* Dialog Window */
	dialog_window = gtk_dialog_new();
	g_signal_connect(G_OBJECT(dialog_window), "destroy",
	                 G_CALLBACK(gtk_main_quit), NULL);
	gtk_window_set_title(GTK_WINDOW(dialog_window), "Dialog -- Exit");
	gtk_widget_set_usize(dialog_window, dialog_size_x, dialog_size_y);

	/* Messeage */
	label = gtk_label_new ("Realy Quit ?");
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog_window)->vbox), 
	                   label, TRUE, TRUE, 0);

	/* YES Botton */
	button_y = gtk_button_new_with_label("YES");
	g_signal_connect(G_OBJECT(button_y), "clicked",
 	                 G_CALLBACK(dialog_yes), &flag);
	g_signal_connect_swapped(G_OBJECT(button_y), "clicked",
	                         G_CALLBACK(gtk_widget_destroy), dialog_window);
	GTK_WIDGET_SET_FLAGS(button_y, GTK_CAN_DEFAULT);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog_window)->action_area), 
	                   button_y, TRUE, TRUE, 0);
	gtk_widget_grab_default(button_y);

	/* NO Botton */
	button_n = gtk_button_new_with_label ("NO");
	g_signal_connect(G_OBJECT (button_n), "clicked",
	                 G_CALLBACK(dialog_no), &flag);
	g_signal_connect_swapped(G_OBJECT(button_n), "clicked",
	                         G_CALLBACK(gtk_widget_destroy), dialog_window);
	GTK_WIDGET_SET_FLAGS(button_n, GTK_CAN_DEFAULT);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog_window)->action_area),
	                   button_n, TRUE, TRUE, 0);

	gtk_window_set_modal(GTK_WINDOW(dialog_window), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(dialog_window), GTK_WINDOW(parent));
	gdk_window_get_position(parent->window, &x, &y);
	gdk_window_get_size(parent->window, &w, &h);
	gtk_widget_set_uposition(dialog_window, x+(w-dialog_size_x)/2,
	                                        y+(h-dialog_size_y)/2);

	gtk_widget_show_all(dialog_window);
	gtk_main();

	return flag;
}


GtkWidget *gtku_create_gtkimage_from_xpm(GtkWidget *widget, gchar **xpm_data)
{
	GdkBitmap *mask;
	GdkPixmap *pixmap_data;
	GtkWidget *pixmap_widget;

	g_return_val_if_fail(widget->window, NULL);

	pixmap_data = gdk_pixmap_create_from_xpm_d(
	               widget->window, &mask,
	               &widget->style->bg[GTK_STATE_NORMAL],
	               (gchar **)xpm_data);
	pixmap_widget = gtk_image_new_from_pixmap(pixmap_data, mask);
	gtk_widget_show(pixmap_widget);

	return(pixmap_widget);
}


void create_gtk_optionmenu_item(GtkWidget *menu, GSList *group,
                                GtkSignalFunc func, const gchar *menu_text,
                                gpointer cb_data)
{
	GtkWidget*  menuitem;
	
	menuitem = gtk_radio_menu_item_new_with_label(group, menu_text);
	group = gtk_radio_menu_item_group(GTK_RADIO_MENU_ITEM(menuitem));
	gtk_menu_append(GTK_MENU(menu), menuitem);
	gtk_signal_connect(GTK_OBJECT(menuitem), "activate", GTK_SIGNAL_FUNC(func),
	                   (gpointer)cb_data);
	gtk_widget_show(menuitem);
}

