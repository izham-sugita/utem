/*********************************************************************
 * adv_component.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Yasuo ASAGA>
 *
 * Thanks to following contributors.
 *  <Taiki AOYAMA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: adv_component.c 410 2005-06-27 07:16:06Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "rc.h"
#include "fem_struct.h"
#include "Adv/AdvDocument.h"
#include "adv_component.h"

#define BUFF  (256)
#define ADVFILE_EXT   ".adv"
#define ADV_MODEL     "./model/advhddm_in"
#define ADV_RESULT    "./result/advhddm_out"
#define ADVFILE_TEMP  "adventure_shape_temp_file.adv"


static RC init_adv_model(const char *adv_model);
static RC init_adv_result(const char *adv_result);
static RC get_num_items_orig(const char *adv_model, const char *content,
                             const char *label, int *num_items_orig);
static RC get_num_parts(const char *adv_model, int *num_parts);
static RC make_node_index_p2g(const char *adv_model, int ipart, int  **ind_p2g);
static RC make_node_index_s2p(const char *adv_model, int ipart, int **ind_s2p);
static RC make_element_index_s2g(const char *adv_model, int ipart, int **ind_s2g);
static RC make_valid_node_index( FEM_NODE_ARRAY node, int **index);


RC adv_input_disp(const char *adv_model, const char *adv_result,
                  FEM_DISP_ARRAY *disp)
{
	int ipart, ii2, ii3;
	int ncount = 0, num_global = -1, num_items_orig = 0, num_parts = 0;
	int num_subdomains = 0, num_sub_items = 0;
	int *ind_s2p = NULL, *ind_p2g = NULL;
	char fname[BUFF];
	AdvDatabox *dbox = NULL;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	RC_TRY( get_num_items_orig(adv_model, "FEGenericAttribute",
	                           "NodeIndex_PartToGlobal", &num_items_orig) );
	RC_TRY( get_num_parts(adv_model, &num_parts) );
	RC_NEG_ZERO_CHK( num_items_orig );
	RC_NEG_ZERO_CHK( num_parts );
	RC_TRY( allocate_fem_disp_array(num_items_orig, disp) );
	disp->source = FEM_ADV;
	disp->size = num_items_orig;
	for(ipart=0; ipart<num_parts; ipart++){
		RC_TRY( make_node_index_p2g(adv_model, ipart, &ind_p2g) );
		RC_TRY( make_node_index_s2p(adv_model, ipart, &ind_s2p) );
		/***** read result file *****/
		sprintf(fname, "%s_%d%s", adv_result, ipart, ADVFILE_EXT);
		RC_NULL_CHK( dbox = adv_dbox_new() );
		if( !(adv_dbox_add(dbox, fname)) ){
			fprintf(stderr, "[%s]: error\n", fname);
			return(NULL_ERROR_RC);
		}
		doc = adv_dbox_find_by_property(dbox,NULL,"label","Displacement", NULL);
		RC_NULL_CHK( doc );
		adv_dio_get_property_int32(doc, "num_subdomains", &num_subdomains);
		RC_NEG_ZERO_CHK( num_subdomains );
		off = 0;
		ncount = 0;
		for(ii2=0; ii2<num_subdomains; ii2++){
			off += adv_dio_read_int32(doc, off, &num_sub_items);
			RC_NEG_ZERO_CHK( num_sub_items );
			for(ii3=0; ii3<num_sub_items; ii3++){
				RC_NEG_CHK( num_global = ind_p2g[ind_s2p[ncount]] );
				disp->array[num_global].node = num_global;
				off += adv_dio_read_float64(doc, off,
				                            &disp->array[num_global].v.x);
				off += adv_dio_read_float64(doc, off,
				                            &disp->array[num_global].v.y);
				off += adv_dio_read_float64(doc, off,
				                            &disp->array[num_global].v.z);
				ncount++;
			}
		}
		adv_dio_close( doc );
		adv_dbox_close( dbox );
		free( ind_s2p );
		free( ind_p2g );
	}

	RC_TRY( clean_fem_disp_array(disp) );

	return(NORMAL_RC);
}


static RC get_num_items_orig(const char *adv_model, const char *content,
                             const char *label, int *num_items_orig)
{
	char        fname[BUFF];
	AdvDatabox  *dbox = NULL;
	AdvDocument *doc = NULL;

	sprintf(fname, "%s_0%s", adv_model, ADVFILE_EXT);
	RC_NULL_CHK( dbox = adv_dbox_new() );
	if( !(adv_dbox_add(dbox, fname)) ){
		fprintf(stderr, "[%s]: error\n", fname);
		return(OPEN_ERROR_RC);
	}
	doc = adv_dbox_find_by_property(dbox, NULL, "content_type", content,
	                                "label", label, NULL);
	RC_NULL_CHK( doc );
	adv_dio_get_property_int32(doc, "num_items_orig", num_items_orig);
	RC_NEG_ZERO_CHK( *num_items_orig );
	adv_dio_close( doc );
	adv_dbox_close( dbox );

	return(NORMAL_RC);
}


static RC get_num_parts(const char *adv_model, int *num_parts)
{
	char        fname[BUFF];
	AdvDatabox  *dbox = NULL;
	AdvDocument *doc = NULL;
	
	sprintf(fname, "%s_0%s", adv_model, ADVFILE_EXT);
	RC_NULL_CHK( dbox = adv_dbox_new() );
	if( !(adv_dbox_add(dbox, fname)) ){
		fprintf(stderr, "[%s]: error\n", fname);
		return(OPEN_ERROR_RC);
	}
	doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                "DocumentList", "label", "HDDM_FEA_Model",
									NULL);
	RC_NULL_CHK( doc );
	adv_dio_get_property_int32(doc, "num_parts", num_parts);
	RC_NEG_ZERO_CHK( *num_parts );
	adv_dio_close( doc );
	adv_dbox_close( dbox );
	
	return(NORMAL_RC);
}


static RC make_node_index_p2g(const char *adv_model, int ipart, int **ind_p2g)
{
	char        fname[BUFF];
	int         num_items = 0;
	AdvDatabox  *dbox = NULL;
	AdvDocument *doc = NULL;

	/***** make model_file name *****/
	sprintf(fname, "%s_%d%s", adv_model, ipart, ADVFILE_EXT);
	RC_NULL_CHK( dbox = adv_dbox_new() );
	if( !(adv_dbox_add(dbox, fname)) ){
		fprintf(stderr, "[%s]: error\n", fname);
		return(OPEN_ERROR_RC);
	}
	doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                "FEGenericAttribute", "label",
	                                "NodeIndex_PartToGlobal", NULL);
	RC_NULL_CHK( doc );
	adv_dio_get_property_int32(doc, "num_items", &num_items);
	RC_NEG_ZERO_CHK( num_items );

	RC_NULL_CHK( (*ind_p2g) = (int *)malloc(num_items*sizeof(int)) );
	adv_dio_read_int32v(doc, 0, num_items, ind_p2g[0]);

	adv_dio_close( doc );
	adv_dbox_close( dbox );

	return(NORMAL_RC);
}


static RC make_node_index_s2p(const char *adv_model, int ipart, int **ind_s2p)
{
	char        fname[BUFF];
	int         ii1;
	int         ncount = 0, num_sub_items = 0, num_subdomains = 0;
	AdvDatabox  *dbox = NULL;
	AdvDocument *doc = NULL;
	adv_off_t   off = 0;

	/***** make model_file name *****/
	sprintf(fname, "%s_%d%s", adv_model, ipart, ADVFILE_EXT);
	RC_NULL_CHK( dbox = adv_dbox_new() );
	if( !(adv_dbox_add(dbox, fname)) ){
		fprintf(stderr, "[%s]: error\n", fname);
		return(OPEN_ERROR_RC);
	}
	doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                "HDDM_FEGenericAttribute", "label",
	                                "NodeIndex_SubdomainToPart", NULL);
	RC_NULL_CHK( doc );
	adv_dio_get_property_int32(doc, "num_subdomains", &num_subdomains);
	RC_NEG_ZERO_CHK( num_subdomains );
	off = 0;
	ncount = 0;
	off += adv_dio_read_int32(doc, off, &num_sub_items);
	RC_NULL_CHK( (*ind_s2p) = (int *)malloc( num_sub_items*sizeof(int)) );
	off += adv_dio_read_int32v(doc, off, num_sub_items, &((*ind_s2p)[ncount]));
	ncount += num_sub_items;
	for(ii1=0; ii1<num_subdomains-1; ii1++){
		off += adv_dio_read_int32(doc, off, &num_sub_items);
		RC_NEG_ZERO_CHK( num_sub_items );
		RC_NULL_CHK( (*ind_s2p) = (int *)realloc((*ind_s2p),
		                                         (num_sub_items+ncount)*
												 sizeof(int)) );
		off += adv_dio_read_int32v(doc, off, num_sub_items,
		                           &((*ind_s2p)[ncount]));
		ncount += num_sub_items;
	}
	adv_dio_close( doc );
	adv_dbox_close( dbox );

	return(NORMAL_RC);
}


RC adv_input_stress(const char *adv_model, const char *adv_result,
                    FEM_STRESS_ARRAY *stress)
{
	char fname[BUFF];
	int ii2, ii3, ipart;
	int ncount = 0, num_global = -1, num_items_orig = 0, num_parts = 0;
	int num_subdomains = 0, num_sub_items = 0;
	int *ind_s2g = NULL;
	AdvDatabox  *dbox = NULL;
	AdvDocument *doc = NULL;
	adv_off_t   off = 0;

	RC_TRY( get_num_items_orig(adv_model, "HDDM_FEGenericAttribute",
	                           "ElementIndex_SubdomainToGlobal",
	                           &num_items_orig) );
	RC_TRY( get_num_parts(adv_model, &num_parts) );
	RC_NEG_ZERO_CHK( num_items_orig );
	RC_NEG_ZERO_CHK( num_parts );
	RC_TRY( allocate_fem_stress_array(num_items_orig, stress) );
	stress->source = FEM_ADV;
	stress->size = num_items_orig;

	for(ipart=0; ipart<num_parts; ipart++){
		RC_TRY( make_element_index_s2g(adv_model, ipart, &ind_s2g) );
		/***** read result file *****/
		sprintf(fname, "%s_%d%s", adv_result, ipart, ADVFILE_EXT);
		RC_NULL_CHK( dbox = adv_dbox_new() );
		if( !(adv_dbox_add(dbox, fname)) ){
			fprintf(stderr, "[%s]: error\n", fname);
			return(NULL_ERROR_RC);
		}
		RC_NULL_CHK( doc = adv_dbox_find_by_property(dbox, NULL, "label",
		                                             "Stress", NULL) );
		adv_dio_get_property_int32(doc, "num_subdomains", &num_subdomains);
		RC_NEG_ZERO_CHK( num_subdomains );
		off = 0;
		ncount = 0;
		for(ii2=0; ii2<num_subdomains; ii2++){
			off += adv_dio_read_int32(doc, off, &num_sub_items);
			RC_NEG_ZERO_CHK( num_sub_items );
			for(ii3=0; ii3<num_sub_items; ii3++){
				RC_NEG_CHK( num_global = ind_s2g[ncount] );
				stress->array[num_global].element = num_global;
				off += adv_dio_read_float64(doc, off,
				                            &stress->array[num_global].v.x);
				off += adv_dio_read_float64(doc, off,
				                            &stress->array[num_global].v.y);
				off += adv_dio_read_float64(doc, off,
				                            &stress->array[num_global].v.z);
				off += adv_dio_read_float64(doc, off,
				                            &stress->array[num_global].v.xy);
				off += adv_dio_read_float64(doc, off,
				                            &stress->array[num_global].v.yz);
				off += adv_dio_read_float64(doc, off,
				                            &stress->array[num_global].v.zx);
				ncount++;
			}
		}
		adv_dio_close( doc );
		adv_dbox_close( dbox );
		free( ind_s2g );
	}

	RC_TRY( clean_fem_stress_array(stress) );

	return(NORMAL_RC);
}


static RC make_element_index_s2g(const char *adv_model, int ipart,
                                 int **ind_s2g)
{
	char        fname[BUFF];
	int         ii1;
	int         ncount = 0;
	int         num_subdomains = 0, num_sub_items = 0, sum_items = 0;
	AdvDatabox  *dbox = NULL;
	AdvDocument *doc = NULL;
	adv_off_t   off = 0;

	sprintf(fname, "%s_%d%s", adv_model, ipart, ADVFILE_EXT);
	/***** read model file *****/
	RC_NULL_CHK( dbox = adv_dbox_new() );
	if( !(adv_dbox_add(dbox, fname)) ){
		fprintf(stderr, "[%s]: error\n", fname);
		return(OPEN_ERROR_RC);
	}
	doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                "HDDM_FEGenericAttribute", "label",
	                                "ElementIndex_SubdomainToGlobal" , NULL);
	RC_NULL_CHK( doc );
	adv_dio_get_property_int32(doc, "num_subdomains", &num_subdomains);
	adv_dio_get_property_int32(doc, "sum_items", &sum_items);
	RC_NEG_ZERO_CHK( num_subdomains );
	RC_NEG_ZERO_CHK( sum_items );
	RC_NULL_CHK( (*ind_s2g) = (int *)malloc(sum_items*sizeof(int)) );
	off = 0;
	ncount = 0;
	for(ii1=0; ii1<num_subdomains; ii1++){
		off += adv_dio_read_int32(doc, off, &num_sub_items);
		RC_NEG_ZERO_CHK( num_sub_items );
		off += adv_dio_read_int32v(doc, off, num_sub_items,
		                           &((*ind_s2g)[ncount]));
		ncount += num_sub_items;
	}
	adv_dio_close( doc );
	adv_dbox_close( dbox );

	return(NORMAL_RC);
}


RC adv_input_strain(const char *adv_model, const char *adv_result,
                    FEM_STRAIN_ARRAY *strain)
{
	char fname[BUFF];
	int ii2, ii3, ipart;
	int ncount = 0, num_global = -1, num_items_orig = 0, num_parts = 0;
	int num_subdomains = 0, num_sub_items = 0;
	int *ind_s2g = NULL;
	AdvDatabox *dbox = NULL;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	RC_TRY( get_num_items_orig(adv_model, "HDDM_FEGenericAttribute",
	                           "ElementIndex_SubdomainToGlobal",
	                           &num_items_orig) );
	RC_TRY( get_num_parts(adv_model, &num_parts) );
	RC_NEG_ZERO_CHK( num_items_orig );
	RC_NEG_ZERO_CHK( num_parts );
	RC_TRY( allocate_fem_strain_array(num_items_orig, strain) );
	strain->source = FEM_ADV;
	strain->size = num_items_orig;

	for(ipart=0; ipart<num_parts; ipart++){
		RC_TRY( make_element_index_s2g(adv_model, ipart, &ind_s2g) );
		/***** read result file *****/
		sprintf(fname, "%s_%d%s", adv_result, ipart, ADVFILE_EXT);
		RC_NULL_CHK( dbox = adv_dbox_new() );
		if( !(adv_dbox_add(dbox, fname)) ){
			fprintf(stderr, "[%s]: error\n", fname);
			return(NULL_ERROR_RC);
		}
		RC_NULL_CHK( doc = adv_dbox_find_by_property(dbox, NULL, "label",
		                                             "Strain", NULL) );
		adv_dio_get_property_int32(doc, "num_subdomains", &num_subdomains);
		RC_NEG_ZERO_CHK( num_subdomains );
		off = 0;
		ncount = 0;
		for(ii2=0; ii2<num_subdomains; ii2++){
			off += adv_dio_read_int32(doc, off, &num_sub_items);
			RC_NEG_ZERO_CHK( num_sub_items );
			for(ii3=0; ii3<num_sub_items; ii3++){
				RC_NEG_CHK( num_global = ind_s2g[ncount] );
				strain->array[num_global].element = num_global;
				off += adv_dio_read_float64(doc, off,
				                            &strain->array[num_global].v.x);
				off += adv_dio_read_float64(doc, off,
				                            &strain->array[num_global].v.y);
				off += adv_dio_read_float64(doc, off,
				                            &strain->array[num_global].v.z);
				off += adv_dio_read_float64(doc, off,
				                            &strain->array[num_global].v.xy);
				off += adv_dio_read_float64(doc, off,
				                            &strain->array[num_global].v.yz);
				off += adv_dio_read_float64(doc, off,
				                            &strain->array[num_global].v.zx);
				ncount++;
			}
		}
		adv_dio_close( doc );
		adv_dbox_close( dbox );
		free( ind_s2g );
	}

	RC_TRY( clean_fem_strain_array(strain) );

	return(NORMAL_RC);
}


RC adv_output_dens_ratio_to_every_part(char *adv_model, char *adv_result,
                                       FEM_ELEMENT_ARRAY elem,
                                       FEM_MATERIAL_PROP_ARRAY dens)
{
	char fname[BUFF];
	char *content_type = "HDDM_FEGenericAttribute";
	char *label = "ElementIndex_SubdomainToGlobal";
	char doc_id[BUFF], doc_id_dens[BUFF];
	int ii1, ii2, ipart;
	int elem_index = -1, num_parts = 0, num_subdomains = 0;
	int num_sub_items = 0,sum_items = 0;
	int *ind_s2g = NULL;
	AdvDatabox *dbox = NULL;
	AdvDocFile *dfile_w = NULL;
	AdvDocument *doc_r = NULL, *doc_w = NULL;
	adv_off_t off_r = 0, off_w = 0;

	RC_TRY( get_num_parts(adv_model, &num_parts) );
	RC_NEG_ZERO_CHK( num_parts );

	for(ipart=0; ipart<num_parts; ipart++){
		/*RC_TRY( make_element_index_s2g(adv_model, ipart, &ind_s2g) );*/
		/***** read result file *****/
		sprintf(fname, "%s_%d%s", adv_model, ipart, ADVFILE_EXT);
		RC_NULL_CHK( dbox = adv_dbox_new() );
		if( !(adv_dbox_add(dbox, fname)) ){
			fprintf(stderr, "[%s]: error\n", fname);
			return(NULL_ERROR_RC);
		}
		RC_NULL_CHK( doc_r = adv_dbox_find_by_property(dbox, NULL,
		             "content_type", content_type, "label", label, NULL) );
		adv_dio_get_property_int32(doc_r, "num_subdomains", &num_subdomains);
		RC_NEG_ZERO_CHK( num_subdomains );
		adv_dio_get_property_int32(doc_r, "sum_items", &sum_items);
		RC_NEG_ZERO_CHK( sum_items );
		sprintf(fname, "%s_%d%s", adv_result, ipart, ADVFILE_EXT);
		RC_NULL_CHK( dfile_w = adv_dio_file_open(fname, "c") );
		sprintf(doc_id, "HDDM_FEGA@HDDM_Part[%d]",ipart);
		doc_w = adv_dio_create(dfile_w,
		                       adv_dio_make_documentid(doc_id));
		RC_NULL_CHK( doc_w );
		off_r = 0;
		off_w = 0;
		for(ii1=0; ii1<num_subdomains; ii1++){
			off_r += adv_dio_read_int32(doc_r, off_r, &num_sub_items);
			RC_NEG_ZERO_CHK( num_sub_items );
			ind_s2g = (int *)malloc(num_sub_items*sizeof(int));
			RC_NULL_CHK( ind_s2g );
			off_r += adv_dio_read_int32v(doc_r, off_r, num_sub_items, ind_s2g);
			adv_dio_set_property(doc_w, "content_type",
			                     "HDDM_FEGenericAttribute");
			adv_dio_set_property_int32(doc_w, "num_subdomains", num_subdomains);
			adv_dio_set_property(doc_w, "fega_type", "AllElementVariable");
			adv_dio_set_property(doc_w, "format", "f8");
			adv_dio_set_property(doc_w, "label", "DensityRatio");
			adv_dio_set_property_int32(doc_w, "index_byte", 4);
			adv_dio_set_property_int32(doc_w, "part_number", ipart);
			adv_dio_set_property(doc_w, "data_type", "s");
			adv_dio_set_property(doc_w, "subindex_type", "");
			adv_dio_set_property_int32(doc_w, "sum_items", sum_items);
			adv_dio_set_property_int32(doc_w, "num_items_orig", dens.size);
			off_w += adv_dio_write_int32(doc_w, off_w, num_sub_items);
			for(ii2=0; ii2<num_sub_items; ii2++){
				elem_index = search_fem_element_label(elem, ind_s2g[ii2]);
				RC_NEG_CHK( elem_index );
				off_w += adv_dio_write_float64(doc_w, off_w,
				                             dens.array[elem_index].rho);
			}
			free( ind_s2g );
		}
		sprintf(doc_id_dens, "%s", adv_dio_get_documentid(doc_w));
		adv_dio_close( doc_w );
		adv_dio_close( doc_r );
		adv_dbox_close( dbox );
		off_w = 0;
		sprintf(doc_id, "DocumentList@HDDM_Part[%d]",ipart);
		doc_w = adv_dio_create(dfile_w, adv_dio_make_documentid(doc_id));
		adv_dio_set_property(doc_w, "content_type", "DocumentList");
		adv_dio_set_property_int32(doc_w, "num_items", 1);
		adv_dio_set_property_int32(doc_w, "num_subdomains", num_subdomains);
		adv_dio_set_property_int32(doc_w, "num_parts", num_parts);
		adv_dio_set_property_int32(doc_w, "part_number", ipart);
		adv_dio_set_property(doc_w, "label", "HDDM_FEA_Result");
		adv_dio_set_property(doc_w, "result_type", "FinalResult");
		off_w += adv_dio_write_string(doc_w, off_w, doc_id_dens);
		adv_dio_close( doc_w );
		adv_dio_file_close( dfile_w );
	}

	return(NORMAL_RC);

}


RC adv_input_dens_ratio(AdvDatabox *dbox, FEM_MATERIAL_PROP_ARRAY *dens)
{
	int ii1;
	int num_items = 0;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	RC_NULL_CHK( doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                             "DensityRatio", NULL) );
	adv_dio_get_property_int32(doc, "num_items", &num_items);
	RC_NEG_ZERO_CHK( num_items );
	RC_TRY( allocate_fem_material_prop_array(num_items, dens) );
	dens->source = FEM_ADV;
	dens->size = num_items;

	for(ii1=0; ii1<(dens->size); ii1++){
		dens->array[ii1].label = ii1;
		off += adv_dio_read_float64(doc, off, &(dens->array[ii1].rho));
	}
	adv_dio_close( doc );
	RC_TRY( clean_fem_material_prop_array(dens) );

	return(NORMAL_RC);
}


RC adv_output_dens_ratio(AdvDocFile *dfile, FEM_MATERIAL_PROP_ARRAY dens)
{
	int ii1;
	int material_valid;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	if( dens.array == NULL ) return(ARG_ERROR_RC);
	RC_NEG_ZERO_CHK( dens.size );

	material_valid = count_valid_material_prop( dens );
	RC_NEG_ZERO_CHK( material_valid );

	doc = adv_dio_create(dfile, adv_dio_make_documentid("DensityRatio"));
	adv_dio_set_property(doc, "content_type", "DensityRatio");
	adv_dio_set_property_int32(doc, "num_items", material_valid);
	adv_dio_set_property(doc, "format", "f8" );

	off = 0;
	for(ii1=0; ii1<dens.size; ii1++){
		if( dens.array[ii1].label < 0 ) continue;
		off += adv_dio_write_float64(doc, off, dens.array[ii1].rho);
	}
	adv_dio_close( doc );

	return(NORMAL_RC);
}


RC adv_input_node(AdvDatabox *dbox, FEM_NODE_ARRAY *node)
{
	int ii1;
	int num_items = 0;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	RC_NULL_CHK( doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                             "Node", NULL) );
	adv_dio_get_property_int32(doc, "num_items", &num_items);
	RC_NEG_ZERO_CHK( num_items );
	RC_TRY( allocate_fem_node_array(num_items, node) );
	node->source = FEM_ADV;
	node->size = num_items;

	off = 0;
	for(ii1=0; ii1<(node->size); ii1++){
		off += adv_dio_read_float64(doc, off, &(node->array[ii1].p.x));
		off += adv_dio_read_float64(doc, off, &(node->array[ii1].p.y));
		off += adv_dio_read_float64(doc, off, &(node->array[ii1].p.z));
		node->array[ii1].label = ii1;
	}
	
	adv_dio_close( doc );
	RC_TRY( clean_fem_node_array(node) );

	return(NORMAL_RC);
}


RC adv_output_node(AdvDocFile *dfile, FEM_NODE_ARRAY node)
{
	int ii1;
	int node_valid = 0;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	if( node.array == NULL ) return(NULL_ERROR_RC);
	RC_NEG_ZERO_CHK( node_valid = count_valid_node(node) );
	RC_NEG_ZERO_CHK( node.size );

	RC_NULL_CHK( doc = adv_dio_create(dfile, adv_dio_make_documentid("Node")) );
	adv_dio_set_property(doc, "content_type", "Node");
	adv_dio_set_property_int32(doc, "num_items", node_valid);
	adv_dio_set_property(doc, "format", "f8f8f8" );
	adv_dio_set_property_int32(doc, "dimension", 3);

	off = 0;
	for(ii1=0; ii1<node.size; ii1++){
		if( node.array[ii1].label < 0 ) continue;
		off += adv_dio_write_float64(doc, off, node.array[ii1].p.x);
		off += adv_dio_write_float64(doc, off, node.array[ii1].p.y);
		off += adv_dio_write_float64(doc, off, node.array[ii1].p.z);
	}
	
	adv_dio_close( doc );

	return(NORMAL_RC);
}


RC adv_input_element(AdvDatabox *dbox, FEM_ELEMENT_ARRAY *elem)
{
	int ii1, ii2;
	int node_num = 0, num_items = 0;
	int i_dat[FEM_MAX_NODE+1];
	char element_type[BUFF];
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	RC_NULL_CHK( doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                             "Element", NULL) );
	adv_dio_get_property_int32(doc, "num_items", &num_items);
	RC_NEG_ZERO_CHK( num_items );
	RC_TRY( allocate_fem_element_array(num_items, elem) );
	elem->source = FEM_ADV;
	elem->size = num_items;

	off = 0;
	adv_dio_get_property_int32(doc, "num_nodes_per_element", &node_num);
	RC_NEG_ZERO_CHK( node_num );
	strcpy( element_type, adv_dio_get_property(doc, "element_type") );
	if( !strcmp(element_type,"3DLinearTetrahedron") ){
		for(ii1=0; ii1<num_items; ii1++){
			off += adv_dio_read_int32v(doc, off, 4, i_dat);
			elem->array[ii1].type     = ELEM_TETRA1;
			elem->array[ii1].node_num = node_num ;
			elem->array[ii1].label    = ii1;
			elem->array[ii1].node[0] = i_dat[0];
			elem->array[ii1].node[1] = i_dat[1];
			elem->array[ii1].node[2] = i_dat[3];
			elem->array[ii1].node[3] = i_dat[2];
			elem->array[ii1].material = 0;
			elem->array[ii1].physical = 0;
		}
	} else if( !strcmp( element_type, "3DQuadraticTetrahedron") ){
		for(ii1=0; ii1<num_items; ii1++){
			off += adv_dio_read_int32v(doc, off, 10, i_dat);
			elem->array[ii1].type     = ELEM_TETRA2;
			elem->array[ii1].node_num = node_num ;
			elem->array[ii1].label 	= ii1;
			elem->array[ii1].node[0]  = i_dat[0];
			elem->array[ii1].node[1]  = i_dat[1];
			elem->array[ii1].node[2]  = i_dat[3];
			elem->array[ii1].node[3]  = i_dat[2];
			elem->array[ii1].node[4]  = i_dat[9];
			elem->array[ii1].node[5]  = i_dat[6];
			elem->array[ii1].node[6]  = i_dat[4];
			elem->array[ii1].node[7]  = i_dat[5];
			elem->array[ii1].node[8]  = i_dat[7];
			elem->array[ii1].node[9]  = i_dat[8];
			elem->array[ii1].material = 0;
			elem->array[ii1].physical = 0;
		}
	} else if( !strcmp( element_type, "3DLinearHexahedron") ){
		for(ii1=0; ii1<num_items; ii1++){
			off += adv_dio_read_int32v(doc, off, 8, i_dat);
			elem->array[ii1].type     = ELEM_HEXA1;
			elem->array[ii1].node_num = node_num ;
			elem->array[ii1].label    = ii1;
			for(ii2=0; ii2<node_num; ii2++){
				elem->array[ii1].node[ii2] = i_dat[ii2];
			}
			elem->array[ii1].material = 0;
			elem->array[ii1].physical = 0;
		}
	} else if( !strcmp( element_type, "3DQuadraticHexahedron") ){
		for(ii1=0; ii1<num_items; ii1++){
			off += adv_dio_read_int32v(doc, off, 20, i_dat);
			elem->array[ii1].type     = ELEM_HEXA2;
			elem->array[ii1].node_num = node_num ;
			elem->array[ii1].label    = ii1;
			for(ii2=0; ii2<12; ii2++){
				elem->array[ii1].node[ii2] = i_dat[ii2];
			}
			elem->array[ii1].node[12] = i_dat[16];
			elem->array[ii1].node[13] = i_dat[17];
			elem->array[ii1].node[14] = i_dat[18];
			elem->array[ii1].node[15] = i_dat[19];
			elem->array[ii1].node[16] = i_dat[12];
			elem->array[ii1].node[17] = i_dat[13];
			elem->array[ii1].node[18] = i_dat[14];
			elem->array[ii1].node[19] = i_dat[15];
			elem->array[ii1].material = 0;
			elem->array[ii1].physical = 0;
		}
	} else{
	}
	adv_dio_close( doc );

	/***** MaterialID *****/
	off = 0;
	doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                "FEGenericAttribute", "label",
	                                "MaterialID", NULL);
	if( doc != NULL ){
		adv_dio_get_property_int32(doc, "num_items", &num_items);
		RC_NEG_ZERO_CHK( num_items );
		for(ii1=0; ii1<num_items; ii1++){
			off += adv_dio_read_int32(doc, off, &(elem->array[ii1].material));
		}
		adv_dio_close( doc );
	}

	RC_TRY( clean_fem_element_array(elem) );

	return(NORMAL_RC);
}


RC adv_output_element(AdvDocFile *dfile, FEM_ELEMENT_ARRAY elem,
                      FEM_NODE_ARRAY node)
{
	int ii1, ii2;
	int element_valid = 0, index = -1;
	int *node_index = NULL;
	AdvDocument *doc = NULL;
	adv_off_t   off = 0;

	if( elem.array == NULL ) return(NULL_ERROR_RC);
	RC_NEG_ZERO_CHK( element_valid = count_valid_element(elem) );
	RC_NEG_ZERO_CHK( elem.size );
	doc = adv_dio_create(dfile, adv_dio_make_documentid("Element"));
	RC_NULL_CHK( doc );
	adv_dio_set_property(doc, "content_type", "Element");
	adv_dio_set_property_int32(doc, "num_items", element_valid);

	off = 0;
	RC_TRY( make_valid_node_index(node, &node_index) );
	for(ii1=0; ii1<elem.size; ii1++){
		if( elem.array[ii1].label < 0 ) continue;
		switch ( elem.array[ii1].type ) {
		case ELEM_TETRA1 :
			adv_dio_set_property(doc, "format", "i4i4i4i4");
			adv_dio_set_property_int32(doc, "num_nodes_per_element", 4);
			adv_dio_set_property(doc, "element_type", "3DLinearTetrahedron");
			adv_dio_set_property_int32(doc, "dimension", 3);
			adv_dio_set_property_int32(doc, "index_byte", 4);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                               elem.array[ii1].node[0])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                               elem.array[ii1].node[1])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                               elem.array[ii1].node[3])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                               elem.array[ii1].node[2])] );
			off += adv_dio_write_int32(doc, off, index);
			break;

		case ELEM_TETRA2 :
			adv_dio_set_property(doc, "format", "i4i4i4i4i4i4i4i4i4i4");
			adv_dio_set_property_int32(doc, "num_nodes_per_element", 10);
			adv_dio_set_property(doc, "element_type", "3DQuadraticTetrahedron");
			adv_dio_set_property_int32(doc, "dimension", 3);
			adv_dio_set_property_int32(doc, "index_byte", 4);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[0])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[1])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[3])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[2])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[6])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[7])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[5])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[8])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[9])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[4])] );
			off += adv_dio_write_int32(doc, off, index);
			break;

		case ELEM_HEXA1 :
			adv_dio_set_property(doc, "format", "i4i4i4i4i4i4i4i4" );
			adv_dio_set_property_int32(doc, "num_nodes_per_element", 8);
			adv_dio_set_property(doc, "element_type", "3DLinearHexahedron");
			adv_dio_set_property_int32(doc, "dimension", 3);
			adv_dio_set_property_int32(doc, "index_byte", 4);
			for(ii2=0; ii2<8; ii2++){
				RC_NEG_CHK( index = node_index[search_fem_node_label(node,
				                    elem.array[ii1].node[ii2])] );
				off += adv_dio_write_int32(doc, off, index);
			}
			break;
		
		case ELEM_HEXA2 :
			adv_dio_set_property(doc, "format", "i4i4i4i4i4i4i4i4i4i4"
			                                     "i4i4i4i4i4i4i4i4i4i4" );
			adv_dio_set_property_int32(doc, "num_nodes_per_element", 20);
			adv_dio_set_property(doc, "element_type", "3DQuadraticHexahedron");
			adv_dio_set_property_int32(doc, "dimension", 3);
			adv_dio_set_property_int32(doc, "index_byte", 4);
			for(ii2=0; ii2<12; ii2++){
				RC_NEG_CHK( index = node_index[search_fem_node_label(node,
				                    elem.array[ii1].node[ii2])] );
				off += adv_dio_write_int32(doc, off, index);
			}
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[16])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[17])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[18])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[19])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[12])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[13])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[14])] );
			off += adv_dio_write_int32(doc, off, index);
			RC_NEG_CHK( index = node_index[search_fem_node_label(node,
			                    elem.array[ii1].node[15])] );
			off += adv_dio_write_int32(doc, off, index);
			break;

		default :
			return(IMPLEMENT_ERROR_RC);
		}
	}
	adv_dio_close( doc );

	free( node_index );
	
	return(NORMAL_RC);
}


RC adv_input_restraint(AdvDatabox *dbox, FEM_BC_ARRAY *rest)
{
	int nflag, num_items = 0, ncount, pre_value;
	int i_dat[2];
	double d_dat;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                                "FEGenericAttribute", "label",
	                                "ForceDisplacement", NULL);
	if( doc == NULL ){
		doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
		                                "FEGenericAttribute", "label",
		                                "ForcedDisplacement", NULL);
	}
	if( doc == NULL ){
		RC_TRY( allocate_fem_bc_array(0, rest) );
		rest->source = FEM_ADV;
		return(NORMAL_RC);
	}
	adv_dio_get_property_int32(doc, "num_items", &num_items);
	RC_NEG_CHK( num_items );
	RC_TRY( allocate_fem_bc_array(num_items, rest) );
	rest->source = FEM_ADV;
	if( num_items == 0 ) return(NORMAL_RC);

	off = 0;
	nflag = 0;
	ncount = 0;
	off += adv_dio_read_int32v(doc, off, 2, i_dat);
	off += adv_dio_read_float64(doc, off, &d_dat);
	while( nflag < num_items ){
		pre_value = i_dat[0];
		rest->array[ncount].node = i_dat[0];  /*** i_dat[0] : Node Number ***/
		do {
			switch ( i_dat[1] ){
			case 0 :       /***** Fixed Node Of X-Direction *****/
				rest->array[ncount].v_type.x = BC_FIX;
				rest->array[ncount].v.x      = d_dat;
				break;
			case 1 :       /***** Fixed Node Of Y-Direction *****/
				rest->array[ncount].v_type.y = BC_FIX;
				rest->array[ncount].v.y      = d_dat;
				break;
			case 2 :       /***** Fixed Node Of Z-Direction *****/
				rest->array[ncount].v_type.z = BC_FIX;
				rest->array[ncount].v.z      = d_dat;
				break;
			default :
				return(READ_ERROR_RC);
			}
			nflag++;
			off += adv_dio_read_int32v(doc, off, 2, i_dat);
			off += adv_dio_read_float64(doc, off, &d_dat);
		} while( pre_value == i_dat[0] );
		ncount ++;
	}
	rest->size = ncount ;

	adv_dio_close( doc );
	RC_TRY( clean_fem_bc_array(rest) );

	return(NORMAL_RC);
}


RC adv_output_restraint(AdvDocFile *dfile, FEM_BC_ARRAY rest,
                        FEM_NODE_ARRAY node)
{
	int ii1;
	int num_rest = 0, index = -1;
	int *node_index = NULL;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	if( rest.size <= 0 ) return(NORMAL_RC);

	num_rest = 0;
	for(ii1=0; ii1<rest.size; ii1++){
		if( rest.array[ii1].node < 0 ) continue;
		if( rest.array[ii1].v_type.x == BC_FIX ) num_rest++;
		if( rest.array[ii1].v_type.y == BC_FIX ) num_rest++;
		if( rest.array[ii1].v_type.z == BC_FIX ) num_rest++;
	}
	RC_NEG_ZERO_CHK( num_rest );

	RC_NULL_CHK( doc = adv_dio_create(dfile,
	             adv_dio_make_documentid("ForcedDisplacement")) );
	adv_dio_set_property(doc, "content_type", "FEGenericAttribute");
	adv_dio_set_property_int32(doc, "num_items", num_rest);
	adv_dio_set_property(doc, "fega_type", "NodeVariable");
	adv_dio_set_property(doc, "label", "ForcedDisplacement");
	adv_dio_set_property(doc, "format", "i4f8");
	adv_dio_set_property_int32(doc, "index_byte", 4);

	off = 0;
	RC_TRY( make_valid_node_index(node, &node_index) );
	for(ii1=0; ii1<rest.size; ii1++){
		if( rest.array[ii1].node < 0 ) continue;
		RC_NEG_CHK( index = search_fem_node_label(node, rest.array[ii1].node) );
		RC_NEG_CHK( index = node_index[index] );
		if( rest.array[ii1].v_type.x == BC_FIX ){
			off += adv_dio_write_int32(doc, off, index);
			off += adv_dio_write_int32(doc, off, 0);
			off += adv_dio_write_float64(doc, off, rest.array[ii1].v.x);
		}
		if( rest.array[ii1].v_type.y == BC_FIX ){
			off += adv_dio_write_int32(doc, off, index);
			off += adv_dio_write_int32(doc, off, 1);
			off += adv_dio_write_float64(doc, off, rest.array[ii1].v.y);
		}
		if( rest.array[ii1].v_type.z == BC_FIX ){
			off += adv_dio_write_int32(doc, off, index);
			off += adv_dio_write_int32(doc, off, 2);
			off += adv_dio_write_float64(doc, off, rest.array[ii1].v.z);
		}
	}
	adv_dio_close( doc );

	free( node_index );

	return(NORMAL_RC);
}


RC adv_input_force(AdvDatabox *dbox, FEM_BC_ARRAY *force)
{
	int nflag, num_items = 0, ncount, pre_value;
	int i_dat[2];
	double       d_dat;
	AdvDocument  *doc = NULL;
	adv_off_t    off = 0;

	doc = adv_dbox_find_by_property(dbox, NULL, "content_type",
	                               "FEGenericAttribute", "label", "Load", NULL);
	if( doc == NULL ){
		RC_TRY( allocate_fem_bc_array(0, force) );
		force->source = FEM_ADV;
		return(NORMAL_RC);
	}

	adv_dio_get_property_int32(doc, "num_items", &num_items);
	RC_NEG_CHK( num_items );
	RC_TRY( allocate_fem_bc_array(num_items, force) );
	force->size = num_items;
	force->source = FEM_ADV;
	if( num_items == 0 ) return(NORMAL_RC);

	off = 0;
	ncount = 0;
	nflag = 0;
	off += adv_dio_read_int32v(doc, off, 2, i_dat);
	off += adv_dio_read_float64(doc, off, &d_dat);
	while( nflag < num_items ){
		pre_value = i_dat[0];
		force->array[ncount].node = i_dat[0];
		do {
			switch ( i_dat[1] ) {
			case 0 :        /***** Load Of X-Direction *****/
				force->array[ncount].v_type.x = BC_FORCE;
				force->array[ncount].v.x      = d_dat;
				break;
			case 1 :        /***** Load Of Y-Direction *****/
				force->array[ncount].v_type.y = BC_FORCE;
				force->array[ncount].v.y      = d_dat;
				break;
			case 2 :        /***** Load Of Z-Direction *****/
				force->array[ncount].v_type.z = BC_FORCE;
				force->array[ncount].v.z      = d_dat;
				break;
			default :
				return(READ_ERROR_RC);
			}
			nflag++;
		off += adv_dio_read_int32v(doc, off, 2, i_dat);
		off += adv_dio_read_float64(doc, off, &d_dat);
		} while( pre_value == i_dat[0] );
		ncount ++;
	}
	adv_dio_close( doc );
	RC_TRY( clean_fem_bc_array(force) );

	return(NORMAL_RC);
}


RC adv_output_force(AdvDocFile *dfile, FEM_BC_ARRAY force, FEM_NODE_ARRAY node)
{
	int ii1;
	int num_force = 0, index = -1;
	int *node_index = NULL;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	if( force.size <= 0 ) return(NORMAL_RC);

	num_force = 0;
	for(ii1=0; ii1<force.size; ii1++){
		if( force.array[ii1].node < 0 ) continue;
		if( force.array[ii1].v_type.x == BC_FORCE ) num_force++;
		if( force.array[ii1].v_type.y == BC_FORCE ) num_force++;
		if( force.array[ii1].v_type.z == BC_FORCE ) num_force++;
	}
	RC_NEG_ZERO_CHK( num_force );

	RC_NULL_CHK( doc = adv_dio_create(dfile, adv_dio_make_documentid("Load")) );
	adv_dio_set_property(doc, "content_type", "FEGenericAttribute");
	adv_dio_set_property_int32(doc, "num_items", num_force);
	adv_dio_set_property(doc, "fega_type", "NodeVariable");
	adv_dio_set_property(doc, "label", "Load");
	adv_dio_set_property(doc, "format", "i4f8" );
	adv_dio_set_property_int32(doc, "index_byte", 4);

	off = 0;
	RC_TRY( make_valid_node_index(node, &node_index) );
	for(ii1=0; ii1<force.size; ii1++){
		if( force.array[ii1].node < 0 ) continue;
		RC_NEG_CHK( index = search_fem_node_label(node, force.array[ii1].node));
		RC_NEG_CHK( index = node_index[index] );
		if( force.array[ii1].v_type.x == BC_FORCE ){
			off += adv_dio_write_int32(doc, off, index);
			off += adv_dio_write_int32(doc, off, 0);
			off += adv_dio_write_float64(doc, off, force.array[ii1].v.x);
		}
		if( force.array[ii1].v_type.y == BC_FORCE ){
			off += adv_dio_write_int32(doc, off, index);
			off += adv_dio_write_int32(doc, off, 1);
			off += adv_dio_write_float64(doc, off, force.array[ii1].v.y);
		}
		if( force.array[ii1].v_type.z == BC_FORCE ){
			off += adv_dio_write_int32(doc, off, index);
			off += adv_dio_write_int32(doc, off, 2);
			off += adv_dio_write_float64(doc, off, force.array[ii1].v.z);
		}
	}
	adv_dio_close( doc );

	free( node_index );

	return(NORMAL_RC);
}


RC adv_input_material(AdvDatabox *dbox, FEM_MATERIAL_PROP_ARRAY *material)
{
	char *content_type = "FEGenericAttribute";
	int nmulti = 0;
	int ncount = 0, num_items = 0;
	double gravity[3];
	AdvDocument *doc = NULL, *prev = NULL;
	adv_off_t off = 0;

	doc = adv_dbox_find_by_property(dbox, NULL, "content_type", content_type,
	                                "label", "MaterialID", NULL);
	if( doc != NULL ){
		nmulti = 1;
		adv_dio_get_property_int32(doc, "num_items", &num_items);
		RC_NEG_CHK( num_items );
		RC_TRY( allocate_fem_material_prop_array(num_items, material) );
		material->size = num_items;
		adv_dio_close( doc );
	}else{
		RC_TRY( allocate_fem_material_prop_array(1, material) );
		material->size = 1;
	}
	material->source = FEM_ADV;

	/***** Density *****/
	ncount = 0;
	prev = NULL;
	while( (doc = adv_dbox_find_by_property(dbox, prev, "content_type",
		    content_type, "label", "Density", NULL)) != NULL ){
		if( nmulti == 1 ){
			adv_dio_get_property_int32(doc, "material_id",
			                           &(material->array[ncount].label) );
		}else{
			material->array[ncount].label = 0;
		}
		adv_dio_read_float64(doc, 0 , &(material->array[ncount].rho) );
		prev = doc;
		ncount++;
		adv_dio_close(doc);
	}

	/***** PoissonRatio *****/
	ncount = 0;
	prev = NULL;
	while( (doc = adv_dbox_find_by_property(dbox, prev, "content_type",
		    content_type, "label", "PoissonRatio", NULL )) != NULL ){
		if( nmulti == 1 ){
			adv_dio_get_property_int32(doc, "material_id",
			                           &(material->array[ncount].label) );
		}else{
			material->array[ncount].label = 0;
		}
		adv_dio_read_float64(doc, 0 , &(material->array[ncount].nu) );
		prev = doc;
		ncount++;
		adv_dio_close(doc);
	}

	/***** YoungModule *****/
	ncount = 0;
	prev = NULL;
	while( (doc = adv_dbox_find_by_property(dbox, prev, "content_type",
		    content_type, "label", "YoungModule", NULL )) != NULL ){
		adv_dio_read_float64(doc, 0 , &(material->array[ncount].E) );
		material->array[ncount].label = 0;
		material->array[ncount].mat_type = MAT_ISOTROPIC;
		material->array[ncount].ana_type = ANA_3D;
		prev = doc;
		ncount++;
		adv_dio_close(doc);
	}
	prev = NULL;
	while( (doc = adv_dbox_find_by_property(dbox, prev, "content_type",
		    content_type, "label", "YoungModulus", NULL )) != NULL ){
		adv_dio_read_float64(doc, 0 , &(material->array[ncount].E) );
		if( nmulti == 1 ){
			adv_dio_get_property_int32(doc, "material_id",
			                           &(material->array[ncount].label) );
		}else{
			material->array[ncount].label = 0;
		}
		material->array[ncount].mat_type = MAT_ISOTROPIC;
		material->array[ncount].ana_type = ANA_3D;
		prev = doc;
		ncount++;
		adv_dio_close(doc);
	}
	
	/***** GravityAcceleration *****/
	ncount = 0;
	doc = adv_dbox_find_by_property(dbox, NULL, "content_type", content_type,
	                               "label", "GravityAcceleration", NULL);
	if( doc != NULL ){
		off = 0;
		off += adv_dio_read_float64v(doc, off, 3, gravity);
		material->array[0].d_info[0] = gravity[0];
		material->array[0].d_info[1] = gravity[1];
		material->array[0].d_info[2] = gravity[2];
		material->array[0].i_info[0] = 1;
		adv_dio_close( doc );
	}else{
		material->array[0].d_info[0] = 0.0;
		material->array[0].d_info[1] = 0.0;
		material->array[0].d_info[2] = 0.0;
		material->array[0].i_info[0] = 0;
	}

	RC_TRY( clean_fem_material_prop_array(material) );

	return(NORMAL_RC);
}


RC adv_output_material(AdvDocFile *dfile, FEM_MATERIAL_PROP_ARRAY material,
                       FEM_ELEMENT_ARRAY elem)
{
	char *content_type = "FEGenericAttribute";
	char *fega_type = "AllElementConstant";
	char *format = "f8";
	int ii1;
	int element_valid = 0;
	AdvDocument *doc = NULL;
	adv_off_t off = 0;

	if( material.size <= 0 ) return(NORMAL_RC);

	/***** Density *****/
	if( material.size > 1 ){
		for(ii1=0; ii1<material.size; ii1++){
			if( material.array[ii1].label < 0 ) continue;
			doc = adv_dio_create(dfile, adv_dio_make_documentid("Density"));
			RC_NULL_CHK( doc );
			adv_dio_set_property(doc, "content_type", content_type);
			adv_dio_set_property_int32(doc, "num_items", 1);
			adv_dio_set_property(doc, "fega_type", fega_type);
			adv_dio_set_property(doc, "label", "Density");
			adv_dio_set_property(doc, "format", format);
			adv_dio_set_property_int32(doc, "material_id", 
			                           material.array[ii1].label);
			adv_dio_set_property_int32(doc, "index_byte", 4);
			adv_dio_write_float64(doc, off, material.array[ii1].rho);
			adv_dio_close( doc );
		}
	}else{
		doc = adv_dio_create(dfile, adv_dio_make_documentid("Density"));
		RC_NULL_CHK( doc );
		adv_dio_set_property(doc, "content_type", content_type);
		adv_dio_set_property_int32(doc, "num_items", 1);
		adv_dio_set_property(doc, "fega_type", fega_type);
		adv_dio_set_property(doc, "label", "Density");
		adv_dio_set_property(doc, "format", format);
		adv_dio_set_property_int32(doc, "index_byte", 4);
		adv_dio_write_float64(doc, off, material.array[0].rho);
		adv_dio_close( doc );
	}

	/***** PoissonRatio *****/
	if( material.size > 1 ){
		for(ii1=0; ii1<material.size; ii1++){
			if( material.array[ii1].label < 0 ) continue;
			doc = adv_dio_create(dfile,adv_dio_make_documentid("PoissonRatio"));
			RC_NULL_CHK( doc );
			adv_dio_set_property(doc, "content_type", content_type);
			adv_dio_set_property_int32(doc, "num_items", 1);
			adv_dio_set_property(doc, "fega_type", fega_type);
			adv_dio_set_property(doc, "label", "PoissonRatio");
			adv_dio_set_property(doc, "format", format);
			adv_dio_set_property_int32(doc, "material_id",
			                           material.array[ii1].label);
			adv_dio_set_property_int32(doc, "index_byte", 4);
			adv_dio_write_float64(doc, 0, material.array[ii1].nu);
			adv_dio_close( doc );
		}
	}else{
		doc = adv_dio_create(dfile,adv_dio_make_documentid("PoissonRatio"));
		RC_NULL_CHK( doc );
		adv_dio_set_property(doc, "content_type", content_type);
		adv_dio_set_property_int32(doc, "num_items", 1);
		adv_dio_set_property(doc, "fega_type", fega_type);
		adv_dio_set_property(doc, "label", "PoissonRatio");
		adv_dio_set_property(doc, "format", format);
		adv_dio_set_property_int32(doc, "index_byte", 4);
		adv_dio_write_float64(doc, 0, material.array[0].nu);
		adv_dio_close( doc );
	}

	/***** YoungModulus *****/
	if( material.size > 1 ){
		for(ii1=0; ii1<material.size; ii1++){
			if( material.array[ii1].label < 0 ) continue;
			doc = adv_dio_create(dfile,adv_dio_make_documentid("YoungModulus"));
			RC_NULL_CHK( doc );
			adv_dio_set_property(doc, "content_type", "FEGenericAttribute");
			adv_dio_set_property_int32(doc, "num_items", 1);
			adv_dio_set_property(doc, "fega_type", fega_type);
			adv_dio_set_property(doc, "label", "YoungModulus");
			adv_dio_set_property(doc, "format", "f8" );
			adv_dio_set_property_int32(doc, "material_id",
			                           material.array[ii1].label);
			adv_dio_set_property_int32(doc, "index_byte", 4);
			adv_dio_write_float64(doc, 0, material.array[ii1].E);
			adv_dio_close( doc );
		}
	}else{
		doc = adv_dio_create(dfile,adv_dio_make_documentid("YoungModulus"));
		RC_NULL_CHK( doc );
		adv_dio_set_property(doc, "content_type", "FEGenericAttribute");
		adv_dio_set_property_int32(doc, "num_items", 1);
		adv_dio_set_property(doc, "fega_type", fega_type);
		adv_dio_set_property(doc, "label", "YoungModulus");
		adv_dio_set_property(doc, "format", "f8" );
		adv_dio_set_property_int32(doc, "index_byte", 4);
		adv_dio_write_float64(doc, 0, material.array[0].E);
		adv_dio_close( doc );
	}

	/***** GravityAcceleration *****/
	if( material.array[0].i_info[0] == 1 ){
		format = "f8f8f8";
		RC_NULL_CHK( doc = adv_dio_create(dfile,
		             adv_dio_make_documentid("GravityAcceleration")) );
		adv_dio_set_property( doc, "content_type", content_type);
		adv_dio_set_property_int32(doc, "num_items", 1);
		adv_dio_set_property(doc, "fega_type", fega_type);
		adv_dio_set_property(doc, "label", "GravityAcceleration");
		adv_dio_set_property(doc, "format", format);
		adv_dio_set_property_int32(doc, "index_byte", 4);
		off = 0;
		off += adv_dio_write_float64(doc, off, material.array[0].d_info[0]);
		off += adv_dio_write_float64(doc, off, material.array[0].d_info[1]);
		off += adv_dio_write_float64(doc, off, material.array[0].d_info[2]);
		adv_dio_close( doc );
	}
	
	/***** MaterialID *****/
	if( material.size > 1 ){
		off = 0;
		element_valid = count_valid_element( elem );
		RC_NEG_ZERO_CHK( element_valid );
		doc = adv_dio_create(dfile, adv_dio_make_documentid("MaterialID"));
		RC_NULL_CHK( doc );
		adv_dio_set_property(doc, "content_type", content_type);
		adv_dio_set_property_int32(doc, "num_items", element_valid);
		adv_dio_set_property(doc, "fega_type", "AllElementVariable");
		adv_dio_set_property(doc, "label", "MaterialID");
		adv_dio_set_property(doc, "format", "i4");
		adv_dio_set_property_int32(doc, "index_byte", 4);
		for(ii1=0; ii1<elem.size; ii1++){
			if( elem.array[ii1].label < 0 ) continue;
			off += adv_dio_write_int32(doc, off, elem.array[ii1].material);
		}
		adv_dio_close( doc );
	}

	return(NORMAL_RC);
}


RC adv_input_physical(AdvDatabox *dbox, FEM_PHYSICAL_PROP_ARRAY *physical)
{
	RC_TRY( allocate_fem_physical_prop_array(1, physical) );
	physical->source = FEM_ADV;
	physical->size = 1;
	physical->array[0].label = 1;
	RC_TRY( clean_fem_physical_prop_array(physical) );
	return(NORMAL_RC);
}


static RC make_valid_node_index(FEM_NODE_ARRAY node, int **node_index)
{
	int ii1, ncount;

	RC_NEG_ZERO_CHK( node.size );
	RC_NULL_CHK( (*node_index) = (int *)malloc(node.size*sizeof(int)) );

	ncount = 0;
	for(ii1=0; ii1<node.size; ii1++){
		if( node.array[ii1].label >= 0 ){
			(*node_index)[ii1] = ncount;
			ncount++;
		} else{
			(*node_index)[ii1] = -1;
		}
	}

	return(NORMAL_RC);
}


static RC init_adv_model(const char *adv_model)
{
	adv_model = ADV_MODEL;
	return(NORMAL_RC);
}


static RC init_adv_result(const char *adv_result)
{
	adv_result = ADV_RESULT;
	return(NORMAL_RC);
}


RC adv_fem(const char *adv_model, const char *adv_result,
           const char *adv_metis, const char *adv_solid,
           FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
           FEM_MATERIAL_PROP_ARRAY material, FEM_BC_ARRAY rest,
           FEM_BC_ARRAY force, FEM_DISP_ARRAY *disp, FEM_STRESS_ARRAY *stress,
           FEM_STRAIN_ARRAY *strain)
{
	char       command[BUFF];
	AdvDocFile *dfile = NULL;

	if( adv_model == NULL ) RC_TRY( init_adv_model(adv_model) );
	if( adv_result == NULL ) RC_TRY( init_adv_result(adv_result) );

	RC_NULL_CHK( dfile = adv_dio_file_open(ADVFILE_TEMP,"c") );
	RC_TRY( adv_output_element(dfile, elem, node) );
	RC_TRY( adv_output_node(dfile, node) );
	RC_TRY( adv_output_restraint(dfile, rest, node) );
	RC_TRY( adv_output_force(dfile, force, node) );
	RC_TRY( adv_output_material(dfile, material, elem) );
	adv_dio_file_close( dfile );

	sprintf(command,"/bin/sh %s", adv_metis);
	system(command);

	sprintf(command,"advsolid -conf %s", adv_solid);
	system(command);

	if( disp != NULL )
		RC_TRY( adv_input_disp(adv_model, adv_result, disp) );
	if( stress != NULL )
		RC_TRY( adv_input_stress(adv_model, adv_result, stress) );
	if( strain != NULL )
		RC_TRY( adv_input_strain(adv_model, adv_result, strain) );

	remove(ADVFILE_TEMP);

	return(NORMAL_RC);
}


RC adv_speed_analysis(const char *adv_model, const char *adv_result,
                      const char *adv_metis, const char *adv_solid,
                      int num_case, FEM_NODE_ARRAY node,
                      FEM_ELEMENT_ARRAY elem, FEM_MATERIAL_PROP_ARRAY material,
                      FEM_BC_ARRAY rest, const FEM_BC_ARRAY *trac,
                      FEM_DISP_ARRAY *disp, FEM_STRAIN_ARRAY *strain)
{
	char       command[BUFF];
	int        ii1; 
	AdvDocFile *dfile = NULL;

	if( adv_model == NULL ) RC_TRY( init_adv_model(adv_model) );
	if( adv_result == NULL ) RC_TRY( init_adv_result(adv_result) );

	for(ii1=0; ii1<num_case; ii1++){
		RC_NULL_CHK( dfile = adv_dio_file_open(ADVFILE_TEMP,"c") );
		RC_TRY( adv_output_element(dfile, elem, node) );
		RC_TRY( adv_output_node(dfile, node) );
		RC_TRY( adv_output_restraint(dfile, rest, node) );
		RC_TRY( adv_output_force(dfile, trac[ii1], node) );
		RC_TRY( adv_output_material(dfile, material, elem) );
		adv_dio_file_close( dfile );

		sprintf(command,"/bin/sh %s", adv_metis);
		system(command);

		sprintf(command,"advsolid -conf %s", adv_solid);
		system(command);

		if( disp != NULL )
			RC_TRY( adv_input_disp(adv_model, adv_result, &(disp[ii1])) );
		if( strain != NULL )
			RC_TRY( adv_input_strain(adv_model, adv_result, &(strain[ii1])) );
	}

	remove(ADVFILE_TEMP);

	return(NORMAL_RC);
}


RC adv_input_model(const char *fname, FEM_NODE_ARRAY *node,
                   FEM_ELEMENT_ARRAY *elem, FEM_BC_ARRAY *rest,
                   FEM_BC_ARRAY *force, FEM_MATERIAL_PROP_ARRAY *material,
                   FEM_PHYSICAL_PROP_ARRAY *physical)
{
	AdvDatabox *dbox;

	RC_NULL_CHK( dbox = adv_dbox_new() );
	if( !(adv_dbox_add(dbox, fname)) ){
		fprintf(stderr, "[%s]: error\n", fname);
		return(OPEN_ERROR_RC);
	}
	RC_TRY( adv_input_node(dbox, node) );
	RC_TRY( adv_input_element(dbox, elem) );
	RC_TRY( adv_input_restraint(dbox, rest) );
	RC_TRY( adv_input_force(dbox, force) );
	RC_TRY( adv_input_material(dbox, material) );
	RC_TRY( adv_input_physical(dbox, physical) );
	adv_dbox_close(dbox);

	return(NORMAL_RC);
}


RC adv_output_model(char *fname, FEM_NODE_ARRAY node, FEM_ELEMENT_ARRAY elem,
                    FEM_BC_ARRAY rest, FEM_BC_ARRAY force,
                    FEM_MATERIAL_PROP_ARRAY material)
{
	AdvDocFile *dfile;

	RC_NULL_CHK( dfile = adv_dio_file_open(fname, "c") );
	RC_TRY( adv_output_element(dfile, elem, node) );
	RC_TRY( adv_output_node(dfile, node) );
	RC_TRY( adv_output_restraint(dfile, rest, node) );
	RC_TRY( adv_output_force(dfile, force, node) );
	RC_TRY( adv_output_material(dfile, material, elem) );
	adv_dio_file_close( dfile );

	return(NORMAL_RC);
}






