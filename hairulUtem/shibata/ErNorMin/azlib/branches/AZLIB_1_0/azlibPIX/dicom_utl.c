/*********************************************************************
 * dicom_utl.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA> <Yoshihito KURIZUKA> <Kazuma ADACHI>
 *
 * Thanks to following contributors.
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: dicom_utl.c,v 1.10 2004/01/13 10:20:41 sasaoka Exp $ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "rc.h"
#include "bin_io.h"
#include "dicom_utl.h"
#include "graphic_utl.h"

#define BUFFER_SIZE (512)

static RC read_tag(FILE *fp, DICOM_HEADER *dicom, unsigned long *ul);
static RC read_tag_elem(FILE *fp, char tmp[], unsigned long tag,
                        int vr_flag);
static RC skip_tag_elemOB(FILE *fp);
static RC skip_tag_elemAE(FILE *fp, int vr_flag);
static RC skip_tag_elemNONE(FILE *fp);
static int search_dicom_tag(const int size, const unsigned long target);
static RC check_vr(FILE *fp, int *flag);

static DICOM_ELEM Dicom_elem[] = {
	{0x00000000, UL, "Group Length 0000"},/* $B%0%k!<%WD9$5(B 0000 */

	{0x00020000, UL, "Group_0002_Length"},            /* $B%0%k!<%WD9$5(B 0002 */
	{0x00020001, OB, "File_Meta_Information_Version"},/* $B%U%!%$%k%a%?>pJsHG(B */
	{0x00020002, UI, "Media Stored SOP Class UID "},  /* $BG^BNJ]B8(BSOP$B%/%i%9(BUID */
	{0x00020003, UI, "Media Stored SOP Instance UID "},
	                                          /* $BG^BNJ]B8(BSOP$B%$%s%9%?%s%9(BUID */
	{0x00020010, UI, "Transfer_Syntax_UID "},        /* $BE>Aw9=J8(BUID */
	{0x00020012, UI, "Implementation_Class_UID"},    /* $B<BAu%/%i%9(BUID */
	{0x00020013, SH, " Implementation_Version_Name"},/* $B<BAuHGL>(B */
	
	{0x00080000, UL, "Group_0008_Length"},           /* $B%0%k!<%WD9$5(B 0008 */ 
	{0x00080001, NONE, "Length to End"},
	{0x00080005, CS, "Specific_Character_Set"},      /* $BFCDjJ8;z=89g(B */
	{0x00080008, CS, "Image_Type" },                 /* $B2hA|%?%$%W(B */
	{0x00080010, NONE, "Recognition Code"},
	{0x00080012, DA, "Instance_Creation_Date"},      /* $B%$%s%9%?%s%9:n@.F|IU(B */
	{0x00080013, TM, "Instance_Creation_Time"},      /* $B%$%s%9%?%s%9:n@.;~9o(B */
	{0x00080014, UI, "Instance_Creator_UID"},        /* $B%$%s%9%?%s%9:n@.<T(BUID */
	{0x00080016, UI, "SOP_Class_UID"},               /* SOP$B%/%i%9(BUID */
	{0x00080018, UI, "SOP_Instance_UID"},            /* SOP$B%$%s%9%?%s%9(BUID */
	{0x00080020, DA, "Study_Date"},                  /* $B8!::F|IU(B */
	{0x00080021, DA, "Series_Date"},                 /* $B%7%j!<%:F|IU(B */
	{0x00080022, DA, "Acquisition_Date"},            /* $B<}=8F|IU(B */
	{0x00080023, DA, "Image_Date"},                  /* $B2hA|F|IU(B */
	{0x00080024, DA, "Overlay_Date"},                /* $B%*!<%P!<%l%$F|IU(B */
	{0x00080025, DA, "Curve_Date"},                  /* $B%+!<%VF|IU(B */
	{0x00080030, TM, "Study_Time"},                  /* $B8!::;~9o(B */
	{0x00080031, TM, "Series_Time"},                 /* $B%7%j!<%:;~9o(B */
	{0x00080032, TM, "Acquisition_Time"},            /* $B<}=8;~9o(B */
	{0x00080033, TM, "Image_Time"},                  /* $B2hA|;~9o(B */
	{0x00080034, TM, "Overlay_Time"},                /* $B%*!<%P!<%l%$;~9o(B */
	{0x00080035, TM, "Curve_Time"},                  /* $B%+!<%V;~9o(B */
	{0x00080040, NONE, "Data Set Type"},
	{0x00080041, NONE, "Data Set Subtype"},
	{0x00080042, CS, "Nuclear_Medicine_Series_Type"},/* $B3K0e3X%7%j!<%:%?%$%W(B */
	{0x00080050, SH, "Accession_Number"},            /* $B<uIUHV9f(B */
	{0x00080052, CS, "Query_Retrieve_Level"},        /* $B4V9g$;(B/$B<hF@%l%Y%k(B */
	{0x00080054, AE, "Retrieve_AE_Title"},           /* $B<hF@(BAE$BL>>N(B */
	{0x00080058, UI, "Failed_SOP_Instance_UID_List"},
	                                          /* $B<:GT(BSOP$B%$%s%9%?%s%9(BUID$B%j%9%H(B */
	{0x00080060, CS, "Modality"},                     /* $B%b%@%j%F%#(B */
	{0x00080061, CS, "Modality in Study"},            /* $B8!::$NCf$N%b%@%j%F%#(B */
	{0x00080064, CS, "Conversion_Type"},              /* $BJQ497A<0(B */
	{0x00080068, CS, "Presentation Intent Type"},     /* $BDs<(0U?^%?%$%W(B */
	{0x00080070, LO, "Manufacturer"},                 /* $B@=B$<T(B */
	{0x00080080, LO, "Institution_Name"},             /*  $B;\@_L>(B */
	{0x00080081, ST, "Institution_Address"},          /* $B;\@_=;=j(B */
	{0x00080082, SQ, "Institution_Code_Sequence"},    /* $B;\@_%3!<%I%7!<%1%s%9(B */
	{0x00080090, PN, "Referring_Physicians_Name"},      /* $B>H2q0e;U$NL>A0(B */
	{0x00080092, ST, "Referring_Physicians_Address"},   /* $B>R2p0e;U$N=;=j(B */
	{0x00081010, SH, "Station_Name"},                   /* $B%9%F!<%7%g%sL>(B */
	{0x00081030, LO, "Study_Description"},              /* $B8!::5-=R(B */
	{0x00081040, LO, "Institutional_Department_Name"},  /* $B;\@_ItLgL>(B */
	{0x00081050, PN, "Attending_Physicians_Name"},      /* $B<B;\0e;U$NL>A0(B */
	{0x00081060, PN, "Name_of_Physician_Reading_Study"},/* $B8!::FI1F0e;U$NL>A0(B */
	{0x00081070, PN, "Operators_Name"},                 /* $BA`:n<T$NL>A0(B */
	{0x00081080, LO, "Admitting_Diagnoses_Description"},/* $B<u?G;~?GCG5-=R(B */
	{0x00081090, LO, "Manufacturers_Model_Name"},       /* $B@=B$<T$N%b%G%kL>(B */
	{0x00082120, SH, "Stage_Name"},                     /* $B%9%F!<%8L>(B */
	{0x00082122, IS, "Stage_Number"},                   /* $B%9%F!<%8HV9f(B */
	{0x00082124, IS, "Number_of_Stages"},               /* $B%9%F!<%8$N?t(B */
	{0x00082128, IS, "View_Number"},                    /* $B%S%e!<HV9f(B */
	{0x0008212a, IS, "Number_of_Views_in_Stage"},/* $B%9%F!<%8$NCf$N%S%e!<$N?t(B */

	{0x00100000, UL, "Group_0010_Length"},     /* $B%0%k!<%WD9$5(B 0010 */
	{0x00100010, PN, "Patients_Name"},         /* $B45<T$NL>A0(B */
	{0x00100020, LO, "Patient_ID"},            /* $B45<T(BID */
	{0x00100030, DA, "Patients_Birth_Date"},   /* $B45<T$NCB@8F|(B */
	{0x00100032, DA, "Patients_Birth_Time"},   /* $B45<T$NCB@8;~9o(B */
	{0x00100040, CS, "Patients_Sex"},          /* $B45<T$N@-JL(B */
	{0x00101000, LO, "Other_Patient_IDs"},     /* $B45<T$NB>$N(BID */
	{0x00101010, AS, "Patients_Age"},          /* $B45<T$NG/Na(B */
	{0x00101020, DS, "Patients_Size"},         /* $B45<T$N?HD9(B */
	{0x00101030, DS, "Patients_Weight"},       /* $B45<T$NBN=E(B */
	{0x00104000, LT, "Patient_Comments"},      /* $B45<T%3%a%s%H(B */

	{0x00180000, UL, "Group_0018_Length"},       /* $B%0%k!<%WD9$5(B 0018 */
	{0x00180010, LO, "Contrast_Bolus_Agent"},    /* $BB$1F:^(B/$B%\!<%i%9:^(B */
	{0x00180015, CS, "Body_Part_Examined"},      /* $B8!::$5$l$kIt0L(B */
	{0x00180020, CS, "Scanning_Sequence"},       /* $B%9%-%c%K%s%0%7!<%1%s%9(B */
	{0x00180021, CS, "Sequence_Variant"},        /* $B%7!<%1%s%9JQ7A(B */
	{0x00180022, CS, "Scan_Options"},            /* $B%9%-%c%s%*%W%7%g%s(B */
	{0x00180023, CS, "MR_Acquisition_Type"},     /* MR$B<}=8%?%$%W(B */
	{0x00180024, SH, "Sequence_Name"},           /* $B%7!<%1%s%9L>(B */
	{0x00180050, DS, "Slice_Thickness"},         /* $B%9%i%$%98|$5(B */
	{0x00180060, DS, "KVP"},
	{0x00180081, DS, "Echo_Time"},               /* $B%(%3!<;~4V(B */
	{0x00180083, DS, "Number_of_Averages"},      /* $BJ?6Q2=$N2s?t(B */
	{0x00180086, IS, "Echo_Numbers"},            /* $B%(%3!<HV9f(B */
	{0x00180091, IS, "Echo_Train_Length"},       /* $B%(%3!<%H%l%sD9$5(B */
	{0x00181000, LO, "Device_Serial_Number"},    /* $BAuCV%7%j%"%kHV9f(B */
	{0x00181020, LO, "Software_Versions "},      /* $B%=%U%H%&%'%"HG(B */
	{0x00181030, LO, "Protocol_Name "},          /* $B%W%m%H%3%kL>(B */
	{0x00181100, DS, "Reconstruction_Diameter "},/* $B:F9=@.D>7B(B */
	{0x00181120, DS, "Gantry_Detector_Tilt "},   /* $B2MBf(B $B8!=P4o79(B */
	{0x00181130, DS, "Table_Height"},            /* $B%F!<%V%k9b$5(B */
	{0x00181140, CS, "Rotation_Direction"},      /* $B2sE>J}8~(B */
	{0x00181150, DS, "Exposure_Time"},           /* $BGx<M;~4V(B */
	{0x00181151, IS, "X_ray_Tube_Current "},     /* X$B@~4IEEN.(B */
	{0x00181152, IS, "Exposure "},               /* $BGx<MNL(B */
	{0x00181210, DS, "Convolution_Kernel"},      /* $B%3%s%\%j!<%7%g%s%+!<%M%k(B */
	{0x00181250, SH, "Receiving_Coil"},          /* $B<u?.%3%$%k(B */
	{0x00181251, SH, "Transmitting_Coil"},       /* $BAw?.%3%$%k(B */
	{0x00181314, DS, "Flip_Angle"},              /* $B%U%j%C%W3Q(B */
	{0x00185100, CS, "Patient_Position"},        /* $B45<T0LCV(B */
	{0x00185101, CS, "View_Position"},           /* $B;kLn0LCV(B */

	{0x00200000, UL, "Group_0020_Length"},          /* $B%0%k!<%WD9$5(B 2000 */
	{0x0020000d, UI, "Study_Instance_UID"},         /* $B8!::%$%s%9%?%s%9(BID */
	{0x0020000e, UI, "Series_Instance_UID"},        /* $B%7%j!<%:%$%s%9%?%s%9(BID */
	{0x00200010, SH, "Study_ID"},                   /* $B8!::(BID */
	{0x00200011, IS, "Series_Number"},              /* $B%7%j!<%:HV9f(B */
	{0x00200012, IS, "Acquisition_Number"},         /* $B<}=8HV9f(B */
	{0x00200013, IS, "Image_Number"},               /* $B%$%s%9%?%s%9HV9f(B */
	{0x00200020, CS, "Patient_Orientation"},        /* $B45<TJ}8~(B */
	{0x00200032, DS, "Image_Position_Patient"},     /* $B2hA|0LCV(B($B45<T(B) */
	{0x00200037, DS, "Image_Orientation_Patient"},  /* $B2hA|J}8~(B($B45<T(B) */
	{0x00200052, UI, "Frame_of_Reference_UID"},     /* $B4p=`:BI87O(BUID */
	{0x00200060, CS, "Laterality"},                 /* $B:81&(B */
	{0x00201002, IS, "Images_in_Acquisition"},      /* $B<}=8$NCf$N2hA|(B */
	{0x00201040, LO, "Position_Reference_Indicator"},/* $B0LCV4p=`I8<1(B */
	{0x00201041, DS, "Slice_Location"},              /* $B%9%i%$%90LCV(B */
	{0x00204000, LT, "Image_Comments"},              /* $B2hA|%3%a%s%H(B */

	{0x00280000, UL, "Group_0028_Length"},         /* $B%0%k!<%WD9$5(B 0028 */
	{0x00280002, US, "Samples_per_Pixel"},         /* $B2hAG$"$?$j%5%s%W%k(B */
	{0x00280004, CS, "Photometric_Interpretation"},/* $B8wEYB,Dj2r<a(B */
	{0x00280006, US, "Planar_Configuration"},      /* $BLL9=@.(B */
	{0x00280010, US, "Rows"},                      /* $B2#9T(B */
	{0x00280011, US, "Columns"},                   /* $B=D9T(B */
	{0x00280030, DS, "Pixel_Spacing"},             /* $B2hAG4V3V(B */
	{0x00280034, IS, "Pixel_Aspect_Ratio"},        /* $B%"%9%Z%/%HHf(B */
	{0x00280100, US, "Bits_Allocated"},            /* $B3d$jEv$F%S%C%H(B */
	{0x00280101, US, "Bits_Stored"},               /* $B3JG<%S%C%H(B */
	{0x00280102, US, "High_Bit"},                  /* $B9b0L%S%C%H(B */
	{0x00280103, US, "Pixel_Representation"},      /* $B2hAGI=8=(B */
	{0x00280106, US, "Smallest_Image_Pixel_Value"},/* $B:G>.2hA|2hAGCM(B */
	{0x00280107, US, "Largest_Image_Pixel_Value"}, /* $B:GBg2hA|2hAGCM(B */
	{0x00281050, DS, "Window_Center"},             /* $B%&%#%s%I%&Cf?4(B */
	{0x00281051, DS, "Window_Width"},              /* $B%&%#%s%I%&I}(B */
	{0x00281052, DS, "Rescale_Intercept"},         /* $B%j%9%1!<%k@ZJR(B */
	{0x00281053, DS, "Rescale_Slope"},             /* $B%j%9%1!<%k79<P(B */
                                                
	{0x7FE00000, UL, "Group_7FE0_Length"},         /* $B%0%k!<%WD9$5(B 7FE0 */
	{0x7FE00010, OW, "Pixel_Data"}                 /* $B2hAG%G!<%?(B */
};

static int dicom_elem_size = (sizeof(Dicom_elem)/sizeof(DICOM_ELEM));

RC read_dicom(FILE *fp, DICOM_HEADER *dicom, GRAPHIC *gr)
{
	int pimage;

	init_dicom(dicom);
	init_graphic(gr);
	RC_TRY( read_dicom_header(fp, dicom, &pimage) );
	RC_TRY( set_dicom_header(*dicom, gr) );
	RC_TRY( read_dicom_image(fp, *dicom, gr, pimage) );
	/*
	print_dicom_header(stderr, *dicom);
	RC_TRY( print_dicom_tag(stderr, *dicom) );
	*/
	return(NORMAL_RC);
}

RC read_dicom_header(FILE *fp, DICOM_HEADER *dicom, int *pimage)
{
	RC rc;
	unsigned long ul;
	unsigned int ui;
	unsigned int vr;
	int index;
	char ctmp[DICOM_STR_SIZE];
	char c;
	float ftmp, ftmp1, ftmp2;
	int vr_flag ;

	rc = read_tag(fp, dicom, &ul);
	if(rc == END_RC){
		return(NORMAL_RC);
	}else if(rc != NORMAL_RC){
		return(rc);
	}

	if(ul == Group_0000_Length){
		if(fseek(fp, 128L, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	}else{
		rewind(fp);
	}

	while(1){
		rc = read_tag(fp, dicom, &ul);
		if(rc == END_RC){
			return(NORMAL_RC);
		}else if(rc != NORMAL_RC){
			return(rc);
		}
		RC_TRY( check_vr(fp, &vr_flag) );

		switch(ul){
		case Group_0008_Length_to_End:
			RC_TRY( lread_4bytes(fp, &ul) );
			RC_TRY( lread_4bytes(fp, &ul) );
			dicom->length_to_end = ul;
	        break;

		case Study_Date: 
			RC_TRY( read_tag_elem(fp, dicom->study_date, ul, vr_flag) );
	        break;

		case Study_Time:
			RC_TRY( read_tag_elem(fp, dicom->study_time, ul, vr_flag) );
			break;

		case Modality:
			RC_TRY( read_tag_elem(fp, dicom->modality, ul, vr_flag) );
			break;

		case Manufacturer:
			RC_TRY( read_tag_elem(fp, dicom->manufacturer, ul, vr_flag) );
			break;

		case Station_Name:
			RC_TRY( read_tag_elem(fp, dicom->station_name, ul, vr_flag) );
			break;

		case Patients_Sex:
			RC_TRY( read_tag_elem(fp, dicom->patients_sex, ul, vr_flag) );
			break;

		case Patients_Age:
			RC_TRY( read_tag_elem(fp, dicom->patients_age, ul, vr_flag) );
			break;
			
		case Slice_Thickness:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f", &ftmp) != 1) return(READ_ERROR_RC);
		    dicom->slice_thickness = (double)ftmp;
			break;

		case Gantry_Detector_Tilt:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f", &ftmp) != 1) return(READ_ERROR_RC);
		    dicom->gantry_detector = (double)ftmp;
			break;

		case Table_Height:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f", &ftmp) != 1) return(READ_ERROR_RC);
		    dicom->table_height = (double)ftmp;
			break;

		case Rotation_Direction:
			RC_TRY( read_tag_elem(fp, dicom->rotation_direction, ul, vr_flag) );
			break;

		case Patient_Position:
			RC_TRY( read_tag_elem(fp, dicom->pation_position, ul, vr_flag) );
			break;

		case Patient_Orientation:
			RC_TRY( read_tag_elem(fp, dicom->patient_orientation, ul,vr_flag) );
			break;

		case Smallest_Image_Pixel_Value:
				RC_TRY( lread_2bytes(fp, &ui) );
				RC_TRY( lread_2bytes(fp, &ui) );
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->smallest_image_pixel_value = (unsigned long)ui;
			break;

		case Largest_Image_Pixel_Value:
				RC_TRY( lread_2bytes(fp, &ui) );
				RC_TRY( lread_2bytes(fp, &ui) );
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->largest_image_pixel_value = (unsigned long)ui;
			break;

		case Slice_Location:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f", &ftmp) != 1) return(READ_ERROR_RC);
		    dicom->slice_location = (double)ftmp;
			break;

		case Samples_per_Pixel:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else if(vr_flag == 0){
				RC_TRY( lread_4bytes(fp, &ul) );
			}else{
				return(UNKNOWN_ERROR_RC);
			}
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->samples_per_pixel = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->samples_per_pixel = ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Photometric_Interpretation:
			RC_TRY( read_tag_elem(fp, dicom->photometric_interpretation,
			                      ul, vr_flag) );
			break;

		case Planar_Configuration:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else if(vr_flag == 0){
				RC_TRY( lread_4bytes(fp, &ul) );
			}else{
				return(UNKNOWN_ERROR_RC);
			}
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->planar_configuration = ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->planar_configuration = (unsigned int)ul;
			}else{
				return(READ_ERROR_RC);
			}
			break;

		case Rows:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}
			
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->rows = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->rows = ul;
			}else{
				fprintf(stderr, "unknown size : case [Rows] <");
				return(READ_ERROR_RC);
			}
			break;

		case Columns:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}
			
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->columns = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->columns = ul;
			}else{
				fprintf(stderr, "unknown size : case [Columns] <");
				return(READ_ERROR_RC);
			}
			break;

		case Pixel_Spacing:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			if(sscanf(ctmp, "%f%c%f", &ftmp1,&c,&ftmp2) != 3){
				return(READ_ERROR_RC);
			}
		    dicom->pixcel_spacing_x = (double)ftmp1;
		    dicom->pixcel_spacing_y = (double)ftmp2;
			break;

		case Pixel_Aspect_Ratio:
			RC_TRY( read_tag_elem(fp, ctmp, ul, vr_flag) );
			break;

		case Bits_Allocated:
			if(vr_flag == 1){
			RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}
			
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->bits_allocated = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->bits_allocated = ul;
			}else{
				fprintf(stderr, "Unknown size : case [bits_allocated] <");
				return(READ_ERROR_RC);
			}
			break;

		case Bits_Stored:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}
			
			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->bits_stored = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->bits_stored = ul;
			}else{
				fprintf(stderr, "Unknown size : case [bits_stored <");
				return(READ_ERROR_RC);
			}
			break;

		case Pixel_Representation:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}

			if(ul == 2){
				RC_TRY( lread_2bytes(fp, &ui) );
				dicom->pixel_representation = (unsigned long)ui;
			}else if(ul == 4){
				RC_TRY( lread_4bytes(fp, &ul) );
				dicom->pixel_representation = ul;
			}else{
				fprintf(stderr, "Unknown size : case"
				                " [pixel_representation] <");
				return(READ_ERROR_RC);
			}
			break;

		case Pixel_Data:
			if(vr_flag == 1){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				RC_TRY( lread_4bytes(fp, &ul) );
			}else{
				RC_TRY( lread_4bytes(fp, &ul) );
			}
			
			dicom->pixel_data = ul;
			*pimage = ftell(fp);
			if(fseek(fp, ul, SEEK_CUR) != 0){
				return(SEEK_ERROR_RC);
			}
			break;

		default:
			index = search_dicom_tag(dicom_elem_size, ul);
			if(index >= 0){
				if( (Dicom_elem[index].vr == OB) ||
					(Dicom_elem[index].vr == OW) ||
					(Dicom_elem[index].vr == SQ) ||
					(Dicom_elem[index].vr == UN)){
					RC_TRY( skip_tag_elemOB(fp) );
				}else if(Dicom_elem[index].vr == NONE){
					RC_TRY( skip_tag_elemNONE(fp) );
				}else{
					RC_TRY( skip_tag_elemAE(fp, vr_flag) );
				}
			}else{
				/* $B4q?tHV9f$r;}$D;dE*%G!<%?MWAG$KBP$9$k;CDjA<CV(B */
				if(((ul>>16) % 2) != 0){
					RC_TRY( skip_tag_elemAE(fp, vr_flag) );
				/* $BITL@HV9f%?%0$rI=<($7$F=*N;(B */
				}else{
					fprintf(stderr, "tag = %0#10lx\n", ul);
					return(UNKNOWN_ERROR_RC);
				}
			}
			break;
		}
	}

	return(NORMAL_RC);
}


RC set_dicom_header(DICOM_HEADER dicom, GRAPHIC *gr)
{
	if((dicom.rows < 1) || (dicom.columns < 1)){
		fprintf(stderr,"NO DATA\n");
		return(READ_ERROR_RC);
	}
	
	gr->source = GRAPHIC_DICOM;
	gr->size.x = dicom.columns;
	gr->size.y = dicom.rows;
	gr->size.z = 1;
	gr->resolution.x = dicom.pixcel_spacing_x;
	gr->resolution.y = dicom.pixcel_spacing_y;
	gr->resolution.z = 0.0;
	
	if(strncmp(dicom.photometric_interpretation, "MONOCHROME1", 11) == 0){
		if(dicom.bits_stored == 16){
			gr->type = GRAPHIC_MONO16;
			gr->i_info[0] = 0;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"MONOCHROME2", 11) == 0){
		if(dicom.bits_stored == 16){
			gr->type = GRAPHIC_MONO16;
			gr->i_info[0] = 1;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation,"PALETTE COLOR",13) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 2;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"RGB", 3) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 3;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"HSV", 3) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 4;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"ARGB", 4) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 5;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"CMYK", 4) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 6;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"YBR_FULL", 8) == 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 7;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation ,"YBR_FULL_422", 12)== 0){
		if(dicom.bits_stored == 8){
			gr->type = GRAPHIC_RGB24;
			gr->i_info[0] = 8;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else if(strncmp(dicom.photometric_interpretation , "\0", 1)== 0){
		if(dicom.bits_stored == 16){
			gr->type = GRAPHIC_MONO16;
			gr->i_info[0] = 9;
		}else{
			fprintf(stderr,"Can't read this Graphic Type!\n");
			return(UNKNOWN_ERROR_RC);
		}
	}else{
		fprintf(stderr,"Unknown Graphic Type!\n");
		return(UNKNOWN_ERROR_RC);
	}
	return(NORMAL_RC);
}


RC read_dicom_image(FILE *fp, DICOM_HEADER dicom, GRAPHIC *gr, int pimage)
{
	unsigned int ui;
	int cl;
	
	if(pimage <= 0){
		fprintf(stderr,"Bad Image File Position <");
		return(READ_ERROR_RC);
	}
	if(fseek(fp, pimage, SEEK_SET) != 0) return(SEEK_ERROR_RC);

	if((gr->type == GRAPHIC_RGB24) || (gr->type == GRAPHIC_MONO16)){
		if(gr != NULL) RC_TRY( allocate_graphic_data(gr) );
	}else{
		fprintf(stderr,"Bad File Type ! <");
		return(UNKNOWN_ERROR_RC);
	}
		
	switch(gr->i_info[0]){
	case 1:
	case 9:
		if((dicom.bits_stored == 8) && (dicom.bits_allocated == 8)){
			GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
				if((cl = fgetc(fp)) == EOF) return(READ_ERROR_RC);
				gr->data.mono16[0][ii1][ii2] = (long)cl;
			}GRAPHIC_SIZE_LOOP_2D_END;
		}else if((dicom.bits_stored == 16) && (dicom.bits_allocated == 16)){
			GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
				RC_TRY(lread_2bytes(fp, &ui) );
				gr->data.mono16[0][ii1][ii2] = (long)ui;
			}GRAPHIC_SIZE_LOOP_2D_END;
		}else{
			fprintf(stderr,"Unknown bits_stored & bits_allocated\n");
			return(READ_ERROR_RC);
		}
		break;
	case 3:
		if((dicom.bits_stored == 8) && (dicom.bits_allocated == 8)){
			if(dicom.planar_configuration == 0){ 
				GRAPHIC_SIZE_LOOP_2D(gr->size.y, gr->size.x){
				  if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].r = (unsigned char)cl;
				  if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].g = (unsigned char)cl;
				  if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
				  gr->data.rgb24[0][ii1][ii2].b = (unsigned char)cl;
				}GRAPHIC_SIZE_LOOP_2D_END;
			}else if(dicom.planar_configuration > 0){
				int ii1, ii2;
				for(ii1=0; (unsigned long)ii1<gr->size.y; ii1++){
					for(ii2=0; (unsigned long)ii2<gr->size.x; ii2++){
						if( (cl=fgetc(fp)) == EOF ) return(READ_ERROR_RC);
						gr->data.rgb24[0][ii1][ii2].r = (unsigned char)cl;
					}
					for(ii2=0; (unsigned long)ii2<gr->size.x; ii2++){
						if( (cl=fgetc(fp)) == EOF) return(READ_ERROR_RC);
						gr->data.rgb24[0][ii1][ii2].g = (unsigned char)cl;
					}
					for(ii2=0; (unsigned long)ii2<gr->size.x; ii2++){
						if( (cl=fgetc(fp)) == EOF) return(READ_ERROR_RC);
						gr->data.rgb24[0][ii1][ii2].b = (unsigned char)cl;
					}
				}
			}else{
				fprintf(stderr, "Unknown planar_configuration type($BLL9=@.(B)\n");
				return(UNKNOWN_ERROR_RC);
			}
		}else{
			fprintf(stderr,"Unknown bits_stored & bits_allocated");
			return(READ_ERROR_RC);
		}
		break;
		
	default: 
		fprintf(stderr,"Unknown Error\n");
		return(UNKNOWN_ERROR_RC);
	}

	return(NORMAL_RC);
}


static RC read_tag(FILE *fp, DICOM_HEADER *dicom, unsigned long *ul)
{
	unsigned int ui1, ui2;
	int current_size;

	if(fgetc(fp) == EOF) return(END_RC);
	if(fseek(fp, -1, SEEK_CUR) != 0) return(SEEK_ERROR_RC);

	RC_TRY( lread_2bytes(fp, &ui1) );
	RC_TRY( lread_2bytes(fp, &ui2) );

	*ul = (unsigned long)ui1;
	*ul <<= 16;
	*ul += ui2;
	
	current_size = dicom->realloc_tag_size;
	if(dicom->tag_size == current_size){
		dicom->realloc_tag_size += DICOM_TAG_SIZE;
		dicom->tag = (unsigned long *)realloc((void *)(dicom->tag),
		              sizeof(unsigned long) * (dicom->realloc_tag_size));
		if((dicom->tag) == (unsigned long *)NULL){
			fprintf(stderr,"realloc < ");
			return(ALLOC_ERROR_RC);
	    }
		current_size = dicom->realloc_tag_size;
	}

	dicom->tag[dicom->tag_size] = *ul;
	dicom->tag_size++;

	return(NORMAL_RC);
}


/* OB,OW,SQ,$B$^$?$O(BUN$B$NL@<(E*(BVR$B$r$b$D%G!<%?MWAG(B */
static RC skip_tag_elemOB(FILE *fp)
{
	unsigned long ul;
	unsigned int ui;
	unsigned int vr;
	
	RC_TRY( lread_2bytes(fp, &ui) );
	RC_TRY( lread_2bytes(fp, &vr) );
	RC_TRY( lread_4bytes(fp, &ul) );
	if (ul == 0) return(NORMAL_RC);
	if(fseek(fp, ul, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	return(NORMAL_RC);
}


/* OB,OW,SQ,$B$^$?$O(BUN$B0J30$NL@<(E*(BVR$B$r$b$D%G!<%?MWAG(B */
static RC skip_tag_elemAE(FILE *fp, int vr_flag)
{
	unsigned int ui;
	unsigned int vr;
	unsigned long ul;

	if(vr_flag == 1){
		RC_TRY( lread_2bytes(fp, &vr) );
		RC_TRY( lread_2bytes(fp, &ui) );
		if(fseek(fp, ui, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	}else if(vr_flag == 0){
		RC_TRY( lread_4bytes(fp, &ul) );
		if(fseek(fp, ul, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	}else{
		return(UNKNOWN_ERROR_RC);
	}
	return(NORMAL_RC);
}


/* $B0EL[E*(BVR$B$r$b$D%G!<%?MWAG(B */
static RC skip_tag_elemNONE (FILE * fp)
{
	unsigned long ul;

	RC_TRY( lread_4bytes(fp, &ul) );
	if (ul == 0) return(NORMAL_RC);
	if(fseek(fp, ul, SEEK_CUR) != 0) return(SEEK_ERROR_RC);
	return(NORMAL_RC);
}


static RC read_tag_elem(FILE *fp, char tmp[], unsigned long tag, int vr_flag)
{
	unsigned long ul;
	unsigned int  ui;
	unsigned int  vr;
	int index;

	if (vr_flag == 1){
		index = search_dicom_tag(dicom_elem_size, tag);
		if(index >= 0){
			if( (Dicom_elem[index].vr == OB) ||
				(Dicom_elem[index].vr == OW) ||
				(Dicom_elem[index].vr == SQ) ||
				(Dicom_elem[index].vr == UN)){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				RC_TRY( lread_4bytes(fp, &ul) );
			}else if(Dicom_elem[index].vr == NONE){
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}else{
				RC_TRY( lread_2bytes(fp, &vr) );
				RC_TRY( lread_2bytes(fp, &ui) );
				ul = (unsigned long)ui;
			}
		}
	}else if(vr_flag == 0){
		RC_TRY( lread_4bytes(fp, &ul) );
	}else{
		return(UNKNOWN_ERROR_RC);
	}

	if(ul == 0){
		tmp[0] = (char)'\0';
		return(NORMAL_RC);
	}else if(ul >= DICOM_STR_SIZE){
		if( fread((void *)tmp, sizeof(char), DICOM_STR_SIZE-1, fp)
		    != DICOM_STR_SIZE-1){
			return(READ_ERROR_RC);
		}	
		tmp[DICOM_STR_SIZE-1] = (char)'\0';
		if(fseek(fp, (long)(ul - (DICOM_STR_SIZE-1)), SEEK_CUR) != 0){
			return(SEEK_ERROR_RC);
		}
	}else{
		if(fread((void *)tmp,sizeof(char),ul,fp) != ul) return(READ_ERROR_RC);
		tmp[ul] = (char)'\0';
	}

	return(NORMAL_RC);
}


static int search_dicom_tag(const int size, const unsigned long target)
{
	int found_index = -1;
	int low = 0;
	int high = size;
	int middle;

	while(low <= high){
		middle = (low + high) / 2;
		if(target == Dicom_elem[middle].tag){
			found_index = middle;
			break;
		}else if(target < Dicom_elem[middle].tag){
			high = middle - 1;
		}else{
			low = middle + 1;
		}
	}
	return(found_index);
}


void init_dicom(DICOM_HEADER *dicom)
{
	dicom->tag = NULL;
	dicom->tag_size = 0;
	dicom->realloc_tag_size = 0;
	dicom->length_to_end = 0;
	dicom->study_date[0] = (char)'\0';
	dicom->study_time[0] = (char)'\0';
	dicom->modality[0] = (char)'\0';
	dicom->manufacturer[0] = (char)'\0';
	dicom->station_name[0] = (char)'\0';
	dicom->patients_sex[0] = (char)'\0';
	dicom->patients_age[0] = (char)'\0';
	dicom->slice_thickness = 0.0; 
	dicom->gantry_detector = 0.0;
	dicom->table_height = 0.0;
	dicom->rotation_direction[0] = (char)'\0';
	dicom->pation_position[0] = (char)'\0';
	dicom->patient_orientation[0] = (char)'\0';
	dicom->slice_location = 0.0;
	dicom->samples_per_pixel = 0;
	dicom->photometric_interpretation[0] = (char)'\0';
	dicom->planar_configuration = 0;
	dicom->rows = 0;
	dicom->columns = 0;
	dicom->pixcel_spacing_x = 0.0;
	dicom->pixcel_spacing_y = 0.0;
	dicom->bits_allocated = 0;
	dicom->bits_stored = 0;
	dicom->pixel_representation = 0;
	dicom->smallest_image_pixel_value = 0;
	dicom->largest_image_pixel_value = 0;
	dicom->pixel_data = 0;
}

void print_dicom_header(FILE *fp, DICOM_HEADER dicom)
{
/*
	fprintf(fp, "Length to End = %ld\n", dicom.length_to_end);
	fprintf(fp, "$B8!::F|IU(B = %s\n", dicom.study_date );
	fprintf(fp, "$B8!::;~4V(B = %s\n", dicom.study_time );
	fprintf(fp, "$B%b%@%j%F%#(B = %s\n", dicom.modality);
	fprintf(fp, "$B@=B$<TL>(B = %s\n", dicom.manufacturer);
	fprintf(fp, "$B%9%F!<%7%g%sL>(B = %s\n", dicom.station_name);
	fprintf(fp, "$B45<T@-JL(B = %s\n", dicom.patients_sex);
	fprintf(fp, "$B45<TG/Na(B = %s\n", dicom.patients_age);
	fprintf(fp, "$B%9%i%$%98|$5(B = %8.3f [mm]\n", dicom.slice_thickness);
	fprintf(fp, "$B2MBf(B/$B8!=P4o79$-(B = %8.3f [deg]\n", dicom.gantry_detector);
	fprintf(fp, "$B%F!<%V%k9b$5(B = %8.3f [mm]\n", dicom.table_height);
	fprintf(fp, "$B2sE>J}8~(B = %s: CW:$B;~7W2s$j(B  CC:$BH?;~7W2s$j!@(Bn",
	             dicom.rotation_direction);
	fprintf(fp, "$B45<T0LCV(B = %s: HFP:Head First-Prone($B$&$D$`$-(B)\n"
	            "                 HFS:Head First-Supinei($B6D8~$1(B)\n"
	            "                 HFDR:Head First-Decubitus Right($B2i>2;Q@*(B)\n"
	            "                 HFDL:Head First-Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Right\n"
	            "                 FFP:Feet First_Prone\n"
	            "                 FFS:Feet First_Supine\n",
	             dicom.pation_position);
	fprintf(fp, "$B45<TJ}8~(B = %s: A(anterior:$BA0LL(B), P(posterior:$BGXLL(B)\n"
	            "                 R(right:$B1&(B), L(left:$B:8(B)\n"
	            "                 H(head:$BF,(B), F(foot:$BB-(B)\n",
	             dicom.patient_orientation);
	fprintf(fp, "$B%9%i%$%90LCV(B = %8.3f\n", dicom.slice_location);
	fprintf(fp, "$B2hAG$"$?$j$N%5%s%W%k(B = %ld\n", dicom.samples_per_pixel);
	fprintf(fp, "$B8wEYB,Dj2r<a(B = %s\n", dicom.photometric_interpretation);
	fprintf(fp, "$BLL9=@.(B = %d\n", dicom.planar_configuration);
	fprintf(fp, "$B2#9T(BX size = %ld\n", dicom.rows);
	fprintf(fp, "$B=D9T(BY size = %ld\n", dicom.columns);
	fprintf(fp, "$B2hAG4V3V(BX = %6.4g [mm]\n", dicom.pixcel_spacing_x);
	fprintf(fp, "$B2hAG4V3V(BY = %6.4g [mm]\n", dicom.pixcel_spacing_y);
	fprintf(fp, "$B3d$jEv$F%S%C%H(B = %ld [bit]\n", dicom.bits_allocated);
	fprintf(fp, "$B3JG<%S%C%H(B = %ld [bit]\n", dicom.bits_stored);
	fprintf(fp, "$B2hAGI=8=(B = %ld\n", dicom.pixel_representation);
	fprintf(fp, "$B:G>.2hAG?t(B = %ld\n", dicom.smallest_image_pixel_value);
	fprintf(fp, "$B:GBg2hAG?t(B = %ld\n", dicom.largest_image_pixel_value);
	fprintf(fp, "pixel size[X*Y*ColorBit] = %ld\n", dicom.pixel_data);
*/
	fprintf(fp, "Length to End = %ld\n", dicom.length_to_end);
	fprintf(fp, "Study Date = %s\n", dicom.study_date );
	fprintf(fp, "Study Time = %s\n", dicom.study_time );
	fprintf(fp, "Modality = %s\n", dicom.modality);
	fprintf(fp, "Manufacturer = %s\n", dicom.manufacturer);
	fprintf(fp, "Station Name = %s\n", dicom.station_name);
	fprintf(fp, "Patients Sex = %s\n", dicom.patients_sex);
	fprintf(fp, "Patients Age = %s\n", dicom.patients_age);
	fprintf(fp, "Slice Thickness = %8.3f [mm]\n", dicom.slice_thickness);
	fprintf(fp, "Gantry/Detector Tilt = %8.3f [deg]\n", dicom.gantry_detector);
	fprintf(fp, "Table Height = %8.3f [mm]\n", dicom.table_height);
	fprintf(fp, " Rotation Direction = %s: CW  CC$B!@(Bn",
	             dicom.rotation_direction);
	fprintf(fp, "Patient Position = %s: HFP:Head First-Prone\n"
	            "                 HFS:Head First-Supinei\n"
	            "                 HFDR:Head First-Decubitus Right\n"
	            "                 HFDL:Head First-Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Right\n"
	            "                 FFP:Feet First_Prone\n"
	            "                 FFS:Feet First_Supine\n",
	             dicom.pation_position);
	fprintf(fp, "Patient Orientation = %s: A(anterior), P(posterior)\n"
	            "                 R(right), L(left)\n"
	            "                 H(head), F(foot)\n",
	             dicom.patient_orientation);
	fprintf(fp, "Slice Location = %8.3f\n", dicom.slice_location);
	fprintf(fp, "Samples Per Pixel = %ld\n", dicom.samples_per_pixel);
	fprintf(fp, "Photometric Interpretation = %s\n",
	             dicom.photometric_interpretation);
	fprintf(fp, "Planar Configuration = %d\n", dicom.planar_configuration);
	fprintf(fp, "RowsX Size = %ld\n", dicom.rows);
	fprintf(fp, "ColumnsY Size = %ld\n", dicom.columns);
	fprintf(fp, "Pixel SpacingX = %6.4g [mm]\n", dicom.pixcel_spacing_x);
	fprintf(fp, "Pixel SpacingY = %6.4g [mm]\n", dicom.pixcel_spacing_y);
	fprintf(fp, "Bits Allocated = %ld [bit]\n", dicom.bits_allocated);
	fprintf(fp, "Bits Stored = %ld [bit]\n", dicom.bits_stored);
	fprintf(fp, "Pixel Representation = %ld\n", dicom.pixel_representation);
	fprintf(fp, "Smallest Image Pixel Value = %ld\n",
	             dicom.smallest_image_pixel_value);
	fprintf(fp, "Largest Image Pixel Value = %ld\n",
	             dicom.largest_image_pixel_value);
	fprintf(fp, "pixel size[X*Y*ColorBit] = %ld\n", dicom.pixel_data);
}


RC print_dicom_tag(FILE *fp, DICOM_HEADER dicom)
{
	int ii1, index;
	unsigned long current_num;
	int dicom_elem_size = (sizeof(Dicom_elem)/sizeof(DICOM_ELEM));

	if(dicom.tag == NULL) return(UNKNOWN_ERROR_RC);
	
	current_num = dicom.tag[0]>>16;
	for(ii1=0;ii1<dicom.tag_size;ii1++){
		index = search_dicom_tag(dicom_elem_size, dicom.tag[ii1]);
		if(index >= 0){
			if(current_num != (dicom.tag[ii1]>>16)){
				fprintf(fp,"\n");
				current_num = dicom.tag[ii1]>>16;
			}
			fprintf(fp, "%0#10lx [%s]\n",dicom.tag[ii1],Dicom_elem[index].name);
		}else{
			if((dicom.tag[ii1]>>16 % 2) != 0){
				if(current_num != (dicom.tag[ii1]>>16)){
					fprintf(fp,"\n");
					current_num = dicom.tag[ii1]>>16;
				}
				fprintf(fp, "%0#10lx [$B;dE*MWAG(B]\n", dicom.tag[ii1]);
			}else{
				fprintf(stderr, "tag = %0#10lx\n", dicom.tag[ii1]);
				return(UNKNOWN_ERROR_RC);
			}
		}
	}
	return(NORMAL_RC);
}


static RC check_vr(FILE *fp, int *vr_flag)
{
	int ii1;
	char tmp[3];
	char *vr[] = {"OB","OW","SQ","UN","AE","AS",
	              "AT","CS","DA","DS","DT","FL",
	              "IS","LO","LT","PN","SH","SL",
	              "SS","ST","TM","UI","UL","US","UT"};

	if((tmp[0] = (char)fgetc(fp)) == EOF) return(READ_ERROR_RC);
	if((tmp[1] = (char)fgetc(fp)) == EOF) return(READ_ERROR_RC);

	tmp[2] = (char)'\0';
	if(fseek(fp, -2L, SEEK_CUR) != 0) return(SEEK_ERROR_RC);

	for(ii1=0;ii1<25;ii1++){
		if(strncmp(vr[ii1], tmp, 2) == 0){
			*vr_flag = 1;
			return(NORMAL_RC);
		}
	}
	*vr_flag = 0;
	return(NORMAL_RC);
}


/*
 * Pixel_Representation($B2hAGI=8=(B)$B$,%;%C%H$5$l$F$$$k$H!$$=$N%b%N%/%m(B
 * $B2hA|$O(B 2$B$NJd?tI=8=$r$7$F$$$k!%$=$l$r(BazlibPIX$B8~$1$N(Bmono16$B$KJQ49$9$k(B
 * Pixel_Representation$B$,%;%C%H$5$l$F$$$F$b!$(B2$B$NJd?t$G$J$$%T%/%;%k$,(B
 * $B4^$^$l$F$$$k$3$H$,$"$k$N$GCm0U!%3F%T%/%;%k$4$H$KD4$Y$k$7$+$J$$!%(B
 */ 
RC conv_2complement_to_mono16(GRAPHIC src, GRAPHIC *dest, DICOM_HEADER dicom)
{
	RC_NULL_CHK(src.data.mono16)
	RC_TRY( copy_graphic_info(src, dest) );
	RC_TRY( allocate_graphic_data(dest) );

	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x)
	{
		/* src 2$B$NJd?tI=8=$+$I$&$+(B2$B=E$N%A%'%C%/(B */
		if( (dicom.pixel_representation != 0) &&
			((src.data.mono16[ii1][ii2][ii3] & (1 << 15)) != 0) ){
			dest->data.mono16[ii1][ii2][ii3] =
			  -( (~(src.data.mono16[ii1][ii2][ii3]-1))&(0xFFFF) );
		}else
			dest->data.mono16[ii1][ii2][ii3] = src.data.mono16[ii1][ii2][ii3];
	}GRAPHIC_SIZE_LOOP_3D_END;
	return(NORMAL_RC);
}

