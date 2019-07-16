/*********************************************************************
 * gtk_analysis.c
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

/* $Id: gtk_analysis.c,v 1.4 2003/07/22 08:02:11 sasaoka Exp $ */

#include <stdio.h>
#include <math.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <gtkgl/gtkglarea.h>
#include <GL/gl.h>
#include "gtk_utl.h"
#include "fem_struct.h"
#include "nst_component.h"
#include "rc.h"

