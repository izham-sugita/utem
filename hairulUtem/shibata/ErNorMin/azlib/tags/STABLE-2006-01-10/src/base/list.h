/*********************************************************************
 * list.h
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: list.h 261 2004-05-07 04:05:42Z sasaoka $ */

#ifndef LIST_H
#define LIST_H

#include <stdio.h>
#include <stdlib.h>


typedef struct list{
	void *data;
	char *tag;
	/* <private> */
	int id;
	struct list *prev;
	struct list *next;
} LIST;


LIST *list_add_tag(LIST *list, const char *tag);
LIST *list_new(const char *tag, void *const data); 

LIST *list_append(LIST *list, const char *tag, void *const data); 
LIST *list_prepend(LIST *list, const char *tag, void *const data); 
LIST *list_insert_after(LIST *list, const char *tag, void *const data);
LIST *list_insert_before(LIST *list, const char *tag, void *const data);

LIST *list_last(LIST *list);
LIST *list_first(LIST *list);
LIST *list_next(LIST *list);
LIST *list_prev(LIST *list);
LIST *list_nth(LIST *list, int n);
int list_length(LIST *list);

LIST *list_recent(LIST *list);
LIST *list_old(LIST *list);

LIST *list_remove(LIST *list); 
LIST *list_remove_all(LIST *list); 
LIST *list_remove_tag(LIST *list, const char *tag); 

LIST *list_search_tag(LIST *list, const char *tag);
LIST *list_search(LIST *list, void *data, int (*matchfunc)(LIST *, void *));

LIST *list_sort_id(LIST *list);
LIST *list_sort(LIST *list, int (*compfunc)(void *, void *));

LIST *list_id_renum(LIST *list);

LIST *list_cat(LIST *list1, LIST *list2);

/* RING LIST */
LIST *ring_list_new(int size);

void list_print_info(LIST *list);
void list_print_info_all(LIST *list);


#endif /* LIST_H */
