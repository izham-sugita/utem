/*********************************************************************
 * list.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: list.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"
#include "rc.h"
#include "memory_manager.h"

static int tag_match(LIST *list, void *tag);
static LIST *sort_part(LIST *first, int num, int (*compfunc)(void *, void *));
static int compar_id(void *id1, void *id2);
static LIST *move_before(LIST *list, LIST *new);


LIST *
list_add_tag (LIST *list, const char *tag)
{
	size_t length;

	if(!list) return(NULL);

	if(!tag){
		log_printf(3, "No tag name. (NULL)\n");
		list->tag = NULL;
		return(list);
	}

	length = strlen(tag) + 1;
	if(list->tag){
		RC rc = mm_free(list->tag); 
		if(rc != NORMAL_RC){
			error_printf(99999);
			return(NULL);
		}
	}
	list->tag = mm_alloc(sizeof(char)*length);
	if(!list->tag){
		error_printf(9001, sizeof(char)*length);
		return(NULL);
	}
	strncpy(list->tag, tag, length);
	list->tag[length-1] = '\0';

	return(list);
}


LIST *
list_new (const char *tag, void *const data)
{
	LIST *new = (LIST *)mm_alloc(sizeof(LIST));
	if(new){
		new->data = data;
		if(!(list_add_tag(new, tag))){
			warning_printf(0, "Can't add tag\n");
			return(list_remove(new));
		}
		new->id = 0;
		new->next = NULL;
		new->prev = NULL;
	}else{
		error_printf(9001, sizeof(LIST));
	}
	return(new);
}


/* list : arbitrary position */
LIST * 
list_append (LIST *list, const char *tag, void *const data)
{
	LIST *recent, *last = NULL;
	LIST *new = list_new(tag, data);

	if(!new){
		error_printf(0, "Can't create new list\n");
		return(NULL);
	}

	if(list){
		last = list_last(list);
		if(!last){
			error_printf(0, "Can't find tail of the list\n");
			return(list_remove(new));
		}
		last->next = new;

		recent = list_recent(list);
		if(!recent){
			warning_printf(0, "Can't find recent list\n");
			return(list_remove(new));
		}
		new->id = recent->id + 1;
	}

	new->prev = last;
	new->next = NULL;

	return(new);
}


/* list : arbitrary position */
LIST * 
list_prepend (LIST *list, const char *tag, void *const data)
{
	LIST *recent, *first = NULL;
	LIST *new = list_new(tag, data);

	if(!new){
		error_printf(0, "Can't create new list\n");
		return(NULL);
	}

	if(list){
		first = list_first(list);
		if(!first){
			error_printf(0, "Can't find tail of the list\n");
			return(list_remove(new));
		}
		first->prev = new;

		recent = list_recent(list);
		if(!recent){
			warning_printf(0, "Can't find recent list\n");
			return(list_remove(new));
		}
		new->id = recent->id + 1;
	}

	new->prev = NULL;
	new->next = first;

	return(new);
}


LIST *
list_insert_after (LIST *list, const char *tag, void *const data)
{
	LIST *recent;
	LIST *new = list_new(tag, data);

	if(!list) return(NULL);
	if(!new){
		error_printf(0, "Can't create new list\n");
		return(NULL);
	}

	recent = list_recent(list);
	if(!recent){
		warning_printf(0, "Can't find recent list\n");
		return(list_remove(new));
	}
	new->id = recent->id + 1;

	if(list->next){
		list->next->prev = new;
		new->next = list->next;
	}else{
		new->next = NULL;
	}
	new->prev = list;
	list->next = new;

	return(new);
}


LIST *
list_insert_before (LIST *list, const char *tag, void *const data)
{
	LIST *recent;
	LIST *new = list_new(tag, data);

	if(!list) return(NULL);
	if(!new){
		error_printf(0, "Can't create new list\n");
		return(NULL);
	}

	recent = list_recent(list);
	if(!recent){
		warning_printf(0, "Can't find recent list\n");
		return(list_remove(new));
	}

	new->id = recent->id + 1;
	if(list->prev){
		list->prev->next = new;
		new->prev = list->prev;
	}else{
		new->prev = NULL;
	}
	new->next = list;
	list->prev = new;

	return(new);
}


/* list : arbitrary position */
LIST *
list_first (LIST *list)
{
	if(list){
		void *p = (void *)list;
		while(list->prev){
			list = list->prev;
			if(p == (void *)list) return(NULL); /* ring list */
		}
	}

	return(list);
}


/* list : arbitrary position */
LIST *
list_last (LIST *list)
{
	if(list){
		void *p = (void *)list;
		while(list->next){
			list = list->next;
			if(p == (void *)list) return(NULL); /* ring list */
		}
	}

	return(list);
}


LIST *
list_next (LIST *list)
{
	if(list) return(list->next);
	return(NULL);
}


LIST *
list_prev (LIST *list)
{
	if(list) return(list->prev);
	return(NULL);
}


/* list : arbitrary position */
LIST *
list_nth (LIST *list, int n)
{
	if(!list) return(NULL);
	if(n <= 0) return(NULL);

	list = list_first(list);
	if(!list){
		error_printf(0, "Can't find head of the list\n");
		return(NULL);
	}

	for(n-- ; n>0; n--){
		if(!(list->next)){
			warning_printf(0, "Can't find N-th of the list\n");
			return(NULL);
		}
		list = list->next;
	}

	return(list);
}


/* list : arbitrary position */
int
list_length (LIST *list)
{
	LIST *first;
	int num = 0;

	if(!list) return(0);

	first = list_first(list);
	if(first){ /* link list */
		list = first;
		num = 1;
		while(list->next){
			num++;
			list = list->next;
		}
	}else{ /* ring list */
		void *p = (void *)list;
		do{
			if(!(list->next)){
				warning_printf(0, "This is a broken ring-list\n");
				return(0);
			}
			list = list->next;
			num++;
		}while(p != (void *)list);
	}

	return(num);
}


/* list : arbitrary position */
LIST *
list_recent (LIST *list)
{
	LIST *first, *tmp;

	if(!list) return(list);

	first = list_first(list);
	if(first){ /* link list */
		tmp = first;
		while(tmp->next){
			tmp = tmp->next;
			if(list->id < tmp->id) list = tmp;
		}
	}else{ /* ring list */
		void *p = (void *)list;
		tmp = list;
		do {
			if(!tmp->next){
				warning_printf(0, "This is a broken ring-list\n");
				return(NULL);
			}
			tmp = tmp->next;
			if(list->id < tmp->id) list = tmp;
		} while(p != (void *)tmp);
	}

	return(list);
}


/* list : arbitrary position */
LIST *
list_old (LIST *list)
{
	LIST *first, *tmp;

	if(!list) return(list);

	first = list_first(list);
	if(first){ /* link list */
		tmp = first;
		while(tmp->next){
			tmp = tmp->next;
			if(list->id > tmp->id) list = tmp;
		}
	}else{ /* ring list */
		void *p = (void *)list;
		tmp = list;
		do {
			if(!tmp->next){
				warning_printf(0, "This is a broken ring-list\n");
				return(NULL);
			}
			tmp = tmp->next;
			if(list->id > tmp->id) list = tmp;
		} while(p != (void *)tmp);
	}

	return(list);
}


/* remove corrent list.                           */
/* return neighborhood list or NULL(mean no list) */
LIST *
list_remove (LIST *list)
{
	LIST *ret = NULL;

	if(list){
		if(list->tag){
			; 
			if(mm_free(list->tag) != NORMAL_RC){
				error_printf(99999);
				return(NULL);
			}
		}
		if( (list->next) && (list->prev) ){
			list->next->prev = list->prev;
			list->prev->next = list->next;
			if(list->next == list) /* for ring_list */
				ret = NULL;
			else
				ret = list->next;
		}else if(list->prev){
			list->prev->next = NULL;
			ret = list->prev;
		}else if(list->next){
			list->next->prev = NULL;
			ret = list->next;
		}
		; 
		if(mm_free(list) != NORMAL_RC){
			error_printf(99999);
			return(NULL);
		}
	}

	return(ret);
}


/* list : arbitrary position           */
/* ATTENTION!!  return NULL if success */           
LIST *
list_remove_all (LIST *list)
{
	while(list)
		list = list_remove(list);

	return(list);
}


/* list : arbitrary position */
LIST *
list_remove_tag (LIST *list, const char *tag)
{
	if(!(list && tag )) return(NULL);

	list = list_search_tag(list, tag);
	return( list_remove(list) );
}


/* list : arbitrary position */
LIST *
list_search_tag (LIST *list, const char *tag)
{
	return( list_search(list, (void *)tag, tag_match) );
}


/* list : arbitrary position */
LIST *
list_search (LIST *list, void *data, int (*matchfunc)(LIST *, void *) )
{
	LIST *first;

	if(!(list && data)) return(NULL);

	first = list_first(list);
	if(first){ /* link list */
		list = first;
		while(list){
			if(matchfunc(list, data)) return(list);
			if(list->next)
				list = list->next;
			else
				return(NULL); /* END of List */
		}
	}else{ /* ring list */
		void *p = (void *)list;
		do{
			if(matchfunc(list, data)) return(list);
			if(list->next){
				list = list->next;
			}else{
				warning_printf(0, "This is a broken ring-listt\n");
				return(NULL);
			}
		}while(p != (void *)list);
	}

	return(NULL);
}


static int
tag_match (LIST *list, void *tag)
{
	if(tag)
		if(strcmp(list->tag, (char *)tag) == 0) return(1);

	return(0);
}


LIST *
list_sort_id (LIST *list)
{
	return( list_sort(list, compar_id) );
}


LIST *
list_sort (LIST *list, int (*compfunc)(void *, void *))
{
	return( sort_part(list_first(list), list_length(list), compfunc) );
}


static LIST *
sort_part (LIST *first, int num, int (*compfunc)(void *, void *))
{
	int ii1, small = 0, large = 0;
	LIST *list, *prev_first = NULL;

	if((!first) || num == 0) return(NULL);
	if(num == 1) return(first);

	list = list_next(first);
	for(ii1=0; (ii1<num-1) && list; ii1++){
		int ret = compfunc(first, list);
		if(ret == 1){
			if(small == 0) prev_first = list;
			list = move_before(first, list);
			small++;
		}else if(ret == 0){
			warning_printf(0, "This list includes same ID.\n");
		}else if(ret == -1){
			list = list_next(list);
			large++;
		}else{
			warning_printf(0, "Unknown return value of compare function.\n");
		}
	}

	if(prev_first && (small > 1))
		if( !sort_part(prev_first,  small, compfunc) ) return(NULL);

	if(first->next && (large > 1))
		if( !sort_part(first->next, large, compfunc) ) return(NULL);

	return(first);
}


static int
compar_id (void *id1, void *id2)
{
	if(id1 && id2){
		if((int *)id1 >  (int *)id2) return(1);
		if((int *)id1 == (int *)id2) return(0);
		if((int *)id1 <  (int *)id2) return(-1);
	}

	return(0);
}


static LIST *
move_before (LIST *list, LIST *target)
{
	LIST *remain = NULL;

	if(!(list && target)) return(NULL);

	/* remove target from corrent position */
	if( (target->next) && (target->prev) ){
		remain = target->next;
		target->next->prev = target->prev;
		target->prev->next = target->next;
	}else if(target->prev){
		target->prev->next = NULL;
	}else if(target->next){
		remain = target->next;
		target->next->prev = NULL;
	}

	/* insert target before "list" */
	if(list->prev){
		list->prev->next = target;
		target->prev = list->prev;
	}else{
		target->prev = NULL;
	}
	target->next = list;
	list->prev = target;

	return(remain);
}


/* list : arbitrary position */
LIST *
list_id_renum (LIST *list)
{
	int id = 0;
	LIST *first;

	if(!list) return(NULL);

	first = list_first(list);
	if(first){ /* link list */
		list = first;
		while(list){
			list->id = id;
			id++;
			if(!list->next) break;
			list = list->next;
		}
	}else{ /* ring list */
		void *p = (void *)list;
		do{
			list->id = id;
			id++;
			if(!list->next){
				warning_printf(0, "This is a broken ring-list\n");
				return(NULL);
			}
			list = list->next;
		}while(p != (void *)list);
	} 
	return(list);
}


/* list[1,2] : arbitrary position */
LIST *
list_cat (LIST *list1, LIST *list2)
{
	if(!(list1 &&list2)) return(NULL);

	list1 = list_last(list1);
	if(!list1){
		error_printf(0, "Can't find tail of the list\n");
		return(NULL);
	}

	list2 = list_first(list2);
	if(!list2){
		error_printf(0, "Can't find head of the list\n");
		return(NULL);
	}

	list1->next = list2;
	list2->prev = list1;

	list1 = list_id_renum(list1);

	return(list1);
}


LIST *
ring_list_new (int size)
{
	LIST *list, *first;

	list = list_new(NULL, NULL);
	if(!list){
		error_printf(0, "Can't create new list\n");
		return(NULL);
	}
	first = list;

	while(size){
		list = list_append(list, NULL, NULL);
		if(!list){
			error_printf(0, "Can't append new list\n");
			return(NULL);
		}
		size--;
	}

	list->next = first;
	first->prev = list;

	return(first);
}


void
list_print_info (LIST *list)
{
	if(list){
		fprintf(stderr, "==== ID: %d ====\n", list->id);
		fprintf(stderr, "list  = %p\n", list);
		fprintf(stderr, "data  = %p\n", list->data);
		fprintf(stderr, "prev  = %p\n", list->prev);
		fprintf(stderr, "next  = %p\n", list->next);
		fprintf(stderr, "tag   = %s\n", list->tag);
	}
}


void
list_print_info_all (LIST *list)
{
	LIST *first;

	if(list){
		first = list_first(list);
		if(first){ /* link list */
			list = first;
			while(list){
				list_print_info(list);
				list = list_next(list);
				if(!list->next) break;
			}
		}else{ /* ring list */
			void *p = (void *)list;
			do{
				list_print_info(list);
				if(!list->next){
					warning_printf(0, "This is a broken ring-list\n");
					break;
				}
				list = list->next;
			}while(p != (void *)list);
		} 
	}
}
