/*********************************************************************
 * dicom_utl.h
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

/* $Id: dicom_utl.h,v 1.6 2003/12/02 03:24:37 sasaoka Exp $ */

#ifndef DICOM_UTL_H
#define DICOM_UTL_H

#include <stdio.h>
#include "rc.h"
#include "graphic_utl.h"

#define DICOM_STR_SIZE (32)
#define DICOM_TAG_SIZE (20)

typedef struct {
	unsigned long *tag;
	int tag_size;
	int realloc_tag_size;
	unsigned long length_to_end;
	char study_date[DICOM_STR_SIZE];
	char study_time[DICOM_STR_SIZE];
	char modality[DICOM_STR_SIZE];
	char manufacturer[DICOM_STR_SIZE];
	char station_name[DICOM_STR_SIZE];
	char patients_sex[DICOM_STR_SIZE];
	char patients_age[DICOM_STR_SIZE];
	double slice_thickness;
	double gantry_detector;
	double table_height;
	char rotation_direction[DICOM_STR_SIZE];
	char pation_position[DICOM_STR_SIZE];
	char patient_orientation[DICOM_STR_SIZE];
	double slice_location;
	unsigned long samples_per_pixel;	
	char photometric_interpretation[DICOM_STR_SIZE];	
	unsigned int planar_configuration;	
	unsigned long rows;
	unsigned long columns;
	double pixcel_spacing_x;
	double pixcel_spacing_y;
	unsigned long bits_allocated;
	unsigned long bits_stored;
	unsigned long pixel_representation;
	unsigned long smallest_image_pixel_value;
	unsigned long largest_image_pixel_value;
	unsigned long pixel_data;
} DICOM_HEADER;

typedef enum{
/* OB,OW,SQ$B$^$?$O(BUN */
	OB, /* $B$=$NB>$N%P%$%HNs(B */
	OW, /* $B$=$NB>$N%o!<%INs(B */
	SQ, /* $B9`L\$N%7!<%1%s%9(B */
	UN, /* $BL$CN(B */
/* OB,OW,SQ$B$^$?$O(BUN$B0J30$NL@<(E*(BVR */
/* (VR$BL>(B), (VR$BL>$N0UL#(B)  ($BCM$ND9$5(B)			*/
	AE, /* $B1~MQ%(%s%G%#%F%#(B  (max 16 bytes) */
	AS, /* $BG/NaNs(B            (4bytes char) */
	AT, /* $BB0@-%?%0(B          (4bytes) */
	CS, /* $B%3!<%INs(B          (max 16 bytes) */
	DA, /* $BF|IU(B              (8bytes) */
	DS, /* 10$B?J?tNs(B          (max 16 bytes) */
	DT, /* $BF|;~(B              (max 16 bytes) */
	FL, /* $BC1@:EYIbF0>.?tE@(B  (4bytes) */
	IS, /* $BG\@:EYIbF0>.?tE@(B  (8bytes) */
	LO, /* $BD9Ns(B              (max 16 bytes) */
	LT, /* $BD9%F%-%9%H(B        (max 10240 char) */
	PN, /* $B?ML>(B              (max 64 char) */
	SH, /* $BC;Ns(B              (max 16 char) */
	SL, /* $BId9fIUD9@0?t(B      (4bytes) */
	SS, /* $BId9fIUC;@0?t(B      (2bytes) */
	ST, /* $BC;%F%-%9%H(B        (max 1024 bytes) */
	TM, /* $B;~4V(B              (max 16 bytes) */
	UI, /* $B8DM-<1JL;R(B        (max 64 bytes) */
	UL, /* $BId9fL5$7D9@0?t(B    (4bytes) */
	US, /* $BId9fL5$7C;@0?t(B    (2bytes) */
	UT, /* $BL5@)8B%F%-%9%H(B    (2^32-2 ?)*/
/* NONE $B$O(B VR$BL5$7(B($B0E<(E*(B) */
	NONE /* VR$BL5$7(B */
} VR;

typedef struct {
	unsigned long tag;
	VR vr;
	char name[DICOM_STR_SIZE];
} DICOM_ELEM;

/* 
*  $B;29MJ88%(B: *  DICOM($B0eNE$K$*$1$k%G%#%8%?%k2hA|$HDL?.(B)
*    $B4,(B5: $B%G!<%?9=B$$HId9f2=(B  $BI=(B6.2-1 DICOM$BCMI=8=(B(p.p.13-17)
*    $B4,(B6: $B%G!<%?<-=q(B    6. DICOM$B%G!<%?MWAG$NEPO?(B(p.p.5-21)
*/
/*****************************************************************************
 * DICOM $B%G!<%?MWAG(B                                                          *
 *****************************************************************************/
typedef enum {
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0000_Length = 0x00000000,
 
 /* Group Length to End $B=*$o$j$^$G$ND9$5(B */
 Group_0000_Length_to_End = 0x00000001,
 
 /* Affected SOP Class UID $B1F6A$5$l$k(BSOP$B%/%i%9(BUID */
 Affected_SOP_Class_UID = 0x00000002,
 
 /* Requested SOP Class UID $BMW5a$5$l$?(BSOP$B%/%i%9(BUID */
 Requested_SOP_Class_UID = 0x00000003,
 
 /* Recognition Code $BG'<1%3!<%I(B */
 Recognition_Code = 0x00000010,
 
 /* Command Field $B%3%^%s%INN0h(B */
 Command_Field = 0x00000100,
 
 /* Message ID $B%a%C%;!<%8(BID */
 Message_ID = 0x00000110,
 
 /* Message Id being Responded to $B1~Ez$5$l$F$$$k%a%C%;!<%8(BID */
 Message_Id_being_Responded_to = 0x00000120,
 
 /* Initiator $B5/F0B&(B */
 Initiator = 0x00000200,
 
 /* Receiver $B<u?.B&(B */
 Receiver = 0x00000300,
 
 /* Find Location FIND$B0LCV(B */
 Find_Location = 0x00000400,
 
 /* Move Destination MOVE$B08@h(B */
 Move_Destination = 0x00000600,
 
 /* Priority $BM%@hEY(B */
 Priority = 0x00000700,
 
 /* Data Set Type $B%G!<%?=89g%?%$%W(B */
 Data_Set_Type = 0x00000800,
 
 /* Number of Matches $BHf3S$N?t(B */
 Number_of_Matches = 0x00000850, 
 
 /* Response Sequence Number $B1~Ez%7!<%1%s%9HV9f(B */
 Response_Sequence_Number = 0x00000860,
 
 /* Status $B>uBV(B */
 Status = 0x00000900,
 
 /* Offending Element $B0cH?MWAG(B */
 Offending_Element = 0x00000901,
 
 /* Error Comment $B%(%i!<%3%a%s%H(B */
 Error_Comment = 0x00000902, 
 
 /* Error ID $B%(%i!<(BID */
 Error_ID = 0x00000903,
 
 /* Affected SOP Instance UID $B1F6A$5$l$k(BSOP$B%$%s%9%?%s%9(BUID */
 Affected_SOP_Instance_UID = 0x00001000,
 
 /* Requested SOP Instance UID $BMW5a$5$l$?(Bsop$B%$%s%9%?%s%9(BUID */
 Requested_SOP_Instance_UID = 0x00001001,
 
 /* Event Type ID $B%$%Y%s%H%?%$%W(BID */
 Event_Type_ID = 0x00001002,
 
 /* Attribute Identifier List $BB0@-<1JL;R%j%9%H(B */
 Attribute_Identifier_List = 0x00001005,
 
 /* Action Type ID $BF0:n%?%$%W(BID */
 Action_Type_ID = 0x00001008,
 
 /* Requested SOP Instance UID List  */
 Requested_SOP_Instance_UID_List = 0x00001012,
 
 /* Number of Remaining Sub operations $B;DB8I{A`:n$N?t(B */
 Number_of_Remaining_Sub_operations = 0x00001020,
 
 /* Number of Completed Sub operations $B40N;I{A`:n$N?t(B */
 Number_of_Completed_Sub_operations = 0x00001021,
 
 /* Number of Failed Sub operations $B<:GTI{A`:n$N?t(B */
 Number_of_Failed_Sub_operations = 0x00001022,
 
 /* Number of Warning Sub operations $B7Y9pI{A`:n$N?t(B */
 Number_of_Warning_Sub_operations = 0x00001023,
 
 /* Move Originator Application Entity Title MOVE$BH/9T851~MQ%(%s%F%#%F%#L>>N(B */
 Move_Originator_Application_Entity_Title = 0x00001030,
 
 /* Move Originator Message ID MOVE$BH/9T85%a%C%;!<%8(BID */
 Move_Originator_Message_ID = 0x00001031,
 
 /* Message Set ID $B%a%C%;!<%8=89g(BID */
 Message_Set_ID = 0x00005010,
 
 /* End Message Set ID $B:G8e$N%a%C%;!<%8=89g(BID */
 End_Message_Set_ID = 0x00005020,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0002_Length = 0x00020000,
 
 /* File Meta Information Version $B%U%!%$%k%a%?>pJsHG(B */
 File_Meta_Information_Version = 0x00020001,
 
 /* Media Stored SOP Class UID $BG^BNJ]B8(Bsop$B%/%i%9(BUID */
 Media_Stored_SOP_Class_UID = 0x00020002,
 
 /* Media Stored SOP Instance UID $BG^BNJ]B8(Bsop$B%$%s%9%?%s%9(BUID */
 Media_Stored_SOP_Instance_UID = 0x00020003,
 
 /* Transfer Syntax UID $BE>Aw9=J8(BUID */
 Transfer_Syntax_UID = 0x00020010,
 
 /* Implementation Class UID $B<BAu%/%i%9(BUID */
 Implementation_Class_UID = 0x00020012,
 
 /* Implementation Version Name $B<BAuHGL>(B */
 Implementation_Version_Name = 0x00020013,
 
 /* Source Application Entity Title $BH/@8851~MQ%(%s%F%#%F%#L>>N(B */
 Source_Application_Entity_Title = 0x00020016,
 
 /* Private Information Creator UID $B;dE*>pJs@8@.<T(BUID */
 Private_Information_Creator_UID = 0x00020100,
 
 /* Private Information $B;dE*>pJs(B */
 Private_Information = 0x00020102,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0004_Length = 0x00040000,
 
 /* File set ID $B%U%!%$%k=89g(BID */
 File_set_ID = 0x00041130, 
 
 /* File set Desctiptor File File ID $B%U%!%$%k=89g5-=R%U%!%$%k(BID */
 File_set_Descriptor_File_File_ID = 0x00041141,
 
 /* File set Descriptor File Format $B%U%!%$%k=89g5-=R%U%!%$%k$NFCDjJ8;z=89g(B */
 File_set_Descriptor_File_Format = 0x00041142,
 
 /* Root Directory Entitys First Directory Record Offset */
 /* $B%k!<%H%G%#%l%/%H%j%(%s%F%#%F%#$N:G=i$N%G%#%l%/%H%j%l%3!<%I$N%*%U%;%C%H(B */
 Root_Directory_Entitys_First_Directory_Record_Offset = 0x00041200,
 
 /* Root Directory Entitys Last Directory Record Offset */
 /* $B%k!<%H%G%#%l%/%H%j%(%s%F%#%F%#$N:G8e$N%G%#%l%/%H%j%l%3!<%I$N%*%U%;%C%H(B */
 Root_Directory_Entitys_Last_Directory_Record_Offset = 0x00041202,
 
 /* File set Consistence Flag $B%U%!%$%k=89g0l47@-%U%i%0(B */
 File_set_Consistence_Flag = 0x00041212,
 
 /* Directory Record Sequence $B%G%#%l%/%H%j%l%3!<%I%7!<%1%s%9(B */
 Directory_Record_Sequence = 0x00041220,
 
 /* Next Directory Record Offset $B<!$N%G%#%l%/%H%j%l%3!<%I$N%*%U%;%C%H(B */
 Next_Directory_Record_Offset = 0x00041400,
 
 /* Record In use Flag $B%l%3!<%I;HMQCf%U%i%0(B */
 Record_In_use_Flag = 0x00041410,
 
 /* Referenced Lower level Directory Entity Offset */
 /* $B;2>H2<0L%G%#%l%/%H%j%(%s%F%#%F%#$N%*%U%;%C%H(B */
 Referenced_Lower_level_Directory_Entity_Offset = 0x00041420,
 
 /* Directory Record Type $B%G%#%l%/%H%j%l%3!<%I%?%$%W(B */
 Directory_Record_Type = 0x00041430,
 
 /* Private Record UID $B;dE*%l%3!<%I(BUID */
 Private_Record_UID = 0x00041432,
 
 /* Referenced File ID $B;2>H%U%!%$%k(BID */
 Referenced_File_ID = 0x00041500,
 
 /* Referenced SOP Class UID in File $B%U%!%$%k$NCf$N;2>H(Bsop$B%/%i%9(BUID */
 Referenced_SOP_Class_UID_in_File = 0x00041510,
 
 /* Referenced SOP Instance UID in File */
 /* $B%U%!%$%k$NCf$N;2>H(Bsop$B%$%s%9%?%s%9(BUID */
 Referenced_SOP_Instance_UID_in_File = 0x00041511,
 
 /* Number of References $B;2>H$N?t(B */
 Number_of_References = 0x00041600,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0008_Length = 0x00080000,
 
 /* Length to End  */
 Group_0008_Length_to_End = 0x00080001,
 
 /* Specific Character Set $BFCDjJ8;z=89g(B */
 Specific_Character_Set = 0x00080005,
 
 /* Image Type $B2hA|%?%$%W(B */
 Image_Type = 0x00080008,
 
 /* Recognition Code  */
 Recognition_Code_0008 = 0x00080010,
 
 /* Instance Creation Date $B%$%s%9%?%s%9@8@.F|IU(B */
 Instance_Creation_Date = 0x00080012,
 
 /* Instance Creation Time $B%$%s%9%?%s%9@8@.;~9o(B */
 Instance_Creation_Time = 0x00080013,
 
 /* Instance Creator UID $B%$%s%9%?%s%9@8@.<T(BUID */
 Instance_Creator_UID = 0x00080014,
 
 /* SOP Class UID sop$B%/%i%9(BUID */
 SOP_Class_UID = 0x00080016,
 
 /* SOP Instance UID sop$B%$%s%9%?%s%9(BUID */
 SOP_Instance_UID = 0x00080018,
 
 /* Study Date $B8!::F|IU(B */
 Study_Date = 0x00080020,
 
 /* Series Date $B%7%j!<%:F|IU(B */
 Series_Date = 0x00080021,
 
 /* Acquisition Date $B<}=8F|IU(B */
 Acquisition_Date = 0x00080022,
 
 /* Image Date $B2hA|F|IU(B */
 Image_Date = 0x00080023,
 
 /* Overlay Date $B%*!<%P%l%$F|IU(B */
 Overlay_Date = 0x00080024,
 
 /* Curve Date $B%+!<%VF|IU(B */
 Curve_Date = 0x00080025,
 
 /* Study Time $B8!::;~9o(B */
 Study_Time = 0x00080030,
 
 /* Series Time $B%7%j!<%:;~9o(B */
 Series_Time = 0x00080031,
 
 /* Acquisition Time $B<}=8;~9o(B */
 Acquisition_Time = 0x00080032,
 
 /* Image Time $B2hA|;~9o(B */
 Image_Time = 0x00080033,
 
 /* Overlay Time $B%*!<%P%l%$;~9o(B */
 Overlay_Time = 0x00080034,
 
 /* Curve Time $B%+!<%V;~9o(B */
 Curve_Time = 0x00080035,
 
 /* Data Set Type */
 Data_Set_Type_0008 = 0x00080040,
 
 /* Data Set Subtype */
 Data_Set_Subtype = 0x00080041,
 
 /* Nuclear Medicine Series Type $B3K0e3X%7%j!<%:%?%$%W(B */
 Nuclear_Medicine_Series_Type = 0x00080042,
 
 /* Accession Number $B<uIUHV9f(B */
 Accession_Number = 0x00080050,
 
 /* Query/Retrieve Level $BLd9g$;!?<hF@%l%Y%k(B */
 Query_Retrieve_Level = 0x00080052,
 
 /* Retrieve AE Title $B<hF@(BAE$BL>>N(B */
 Retrieve_AE_Title = 0x00080054,
 
 /* Failed SOP Instance UID List $B<:GT(Bsop$B%$%s%9%?%s%9(BUID$B%j%9%H(B */
 Failed_SOP_Instance_UID_List = 0x00080058,
 
 /* Modality $B%b%@%j%F%#(B */
 Modality = 0x00080060,
 
 /* Conversion Type $BJQ497A<0(B */
 Conversion_Type = 0x00080064,
 
 /* Manufacturer $B@=B$<T(B */
 Manufacturer = 0x00080070,
 
 /* Institution Name $B;\@_L>(B */
 Institution_Name = 0x00080080,
 
 /* Institution Address $B;\@_$N=;=j(B */
 Institution_Address = 0x00080081,
 
 /* Institution Code Sequence $B;\@_%3!<%I%7!<%1%s%9(B */
 Institution_Code_Sequence = 0x00080082,
 
 /* Referring Physician's Name $B>H2q0e;UL>(B */
 Referring_Physicians_Name = 0x00080090,
 
 /* Referring Physician's Address $B>H2q0e;U=;=j(B */
 Referring_Physicians_Address = 0x00080092,
 
 /* Referring Physician's Telephone Numbers $B>H2q0e;UEEOCHV9f(B */
 Referring_Physicians_Telephone_Numbers = 0x00080094,
 
 /* Code Value $B%3!<%ICM(B */
 Code_Value = 0x00080100,
 
 /* Coding Scheme Designator $BId9f2=BN7O;XDj;R(B */
 Coding_Scheme_Designator = 0x00080102,
 
 /* Code Meaning $B%3!<%I0UL#(B */
 Code_Meaning = 0x00080104,
 
 /* Network ID  */
 Network_ID = 0x00081000,
 
 /* Station Name $B%9%F!<%7%g%sL>(B */
 Station_Name = 0x00081010,
 
 /* Study Description $B8!::5-=R(B */
 Study_Description = 0x00081030,
 
 /* Procedure Code Sequence $B=hCV%3!<%I%7!<%1%s%9(B */
 Procedure_Code_Sequence = 0x00081032,
 
 /* Series Description $B%7%j!<%:5-=R(B */
 Series_Description = 0x0008103E,
 
 /* Institutional Department Name $B;\@_ItLgL>(B */
 Institutional_Department_Name = 0x00081040,
 
 /* Performing Physician's Name $B<B;\0e;U$NL>A0(B */
 Attending_Physicians_Name = 0x00081050,
 
 /* Name of  Physician(s) Reading Study $B8!::FI1F0e;UL>(B */
 Name_of_Physician_Reading_Study = 0x00081060,
 
 /* Operators' Name $BA`:n<T$NL>A0(B */
 Operators_Name = 0x00081070,
 
 /* Admitting Diagnoses Description $B<u?G;~?GCG5-=R(B */
 Admitting_Diagnoses_Description = 0x00081080,
 
 /* Admitting Diagnosis Code Sequence $B<u?G;~?GCG%3!<%I%7!<%1%s%9(B */
 Admitting_Diagnosis_Code_Sequence = 0x00081084,
 
 /* Manufacturer's Model Name $B@=B$<T%b%G%kL>(B */
 Manufacturers_Model_Name = 0x00081090,
 
 /* Referenced Results Sequence $B;2>H7k2L%7!<%1%s%9(B */
 Referenced_Results_Sequence = 0x00081100,
 
 /* Referenced Study Sequence $B;2>H8!::%7!<%1%s%9(B */
 Referenced_Study_Sequence = 0x00081110,
 
 /* Referenced Study Component Sequence $B;2>H8!::9=@.MWAG%7!<%1%s%9(B */
 Referenced_Study_Component_Sequence = 0x00081111,
 
 /* Referenced Series Sequence $B;2>H%7%j!<%:%7!<%1%s%9(B */
 Referenced_Series_Sequence = 0x00081115,
 
 /* Referenced Patient Sequence $B;2>H45<T%7!<%1%s%9(B */
 Referenced_Patient_Sequence = 0x00081120,
 
 /* Referenced Visit Sequence $B;2>HMh1!%7!<%1%s%9(B */
 Referenced_Visit_Sequence = 0x00081125,
 
 /* Referenced Overlay Sequence $B;2>H%*!<%P%l%$%7!<%1%s%9(B */
 Referenced_Overlay_Sequence = 0x00081130,
 
 /* Referenced Image Sequence $B;2>H2hA|%7!<%1%s%9(B */
 Referenced_Image_Sequence = 0x00081140,
 
 /* Referenced Curve Sequence $B;2>H%+!<%V%7!<%1%s%9(B */
 Referenced_Curve_Sequence = 0x00081145,
 
 /* Referenced SOP Class UID $B;2>H(Bsop$B%/%i%9(BUID */
 Referenced_SOP_Class_UID = 0x00081150,
 
 /* Referenced SOP Instance UID $B;2>H(Bsop$B%$%s%9%?%s%9(BUID */
 Referenced_SOP_Instance_UID = 0x00081155,
 
 /* Derivation Description $BF3=P5-=R(B */
 Derivation_Description = 0x00082111,
 
 /* Source Image Sequence $BH/@8852hA|%7!<%1%s%9(B */
 Source_Image_Sequence = 0x00082112,
 
 /* Stage Name $B%9%F!<%8L>(B */
 Stage_Name = 0x00082120,
 
 /* Stage Number $B%9%F!<%8HV9f(B */
 Stage_Number = 0x00082122,
 
 /* Number of Stages $B%9%F!<%8?t(B */
 Number_of_Stages = 0x00082124,
 
 /* Number of Event Timers $B%$%Y%s%H%?%$%^$N?t(B */
 Number_of_Event_Timers = 0x00082129,
 
 /* View Number $B%S%e!<HV9f(B */
 View_Number = 0x00082128,
 
 /* Number of Views in Stage $B%9%F!<%8$NCf$N%S%e!<?t(B */
 Number_of_Views_in_Stage = 0x0008212A,
 
 /* Event Elapsed Time(s) $B%$%Y%s%H7P2a;~4V(B */
 Event_Elapsed_Time = 0x00082130,
 
 /* Event Timer Name(s) $B%$%Y%s%H%?%$%^L>(B */
 Event_Timer_Name = 0x00082132,
 
 /* Start Trim $B3+;O%H%j%`(B */
 Start_Trim = 0x00082142,
 
 /* Stop Trim $BDd;_%H%j%`(B */
 Stop_Trim = 0x00082143,
 
 /* Recommended Display Frame Rate $B?d>)I=<(%U%l!<%`N((B */
 Recommended_Display_Frame_Rate = 0x00082144,
 
 /* Transducer Position $BC5?(;R$N0LCV(B */
 Transducer_Position = 0x00082200,
 
 /* Transducer Orientation $BC5?(;R$NJ}8~(B */
 Transducer_Orientation = 0x00082204,
 
 /* Anatomic Structure $B2rK63XE*9=B$(B */
 Anatomic_Structure = 0x00082208, 
 
 /* Comments  */
 Group_0008_Comments = 0x00084000,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0010_Length = 0x00100000,
 
 /* Patient's Name $B45<T$NL>A0(B */
 Patients_Name = 0x00100010,
 
 /* Patient ID $B45<T(BID */
 Patient_ID = 0x00100020,
 
 /* Issuer of Patient ID $B45<T(BID$B$NH/9T<T(B */
 Issuer_of_Patient_ID = 0x00100021,
 
 /* Patient's Birth Date $B45<T$NCB@8F|(B */
 Patients_Birth_Date = 0x00100030,
 
 /* Patient's Birth Time $B45<T$NCB@8;~9o(B */
 Patients_Birth_Time = 0x00100032,
 
 /* Patient's Sex $B45<T$N@-JL(B */
 Patients_Sex = 0x00100040,
 
 /*  Patient's Social Security Number $B45<T$N<R2qJ]>cHV9f(B */
 Patients_Social_Security_Number = 0x00100042,
 
 /* Patient's Insurance Plan Code Sequence */
 /* $B45<T$NJ]817W2h%3!<%I%7!<%1%s%9(B */
 Patients_Insurance_Plan_Code_Sequence = 0x00100050,
 
 /* Other Patient IDs $B45<T$NB>$N(BID */
 Other_Patient_IDs = 0x00101000,
 
 /* Other Patient Names $B45<T$NB>$NL>A0(B */
 Other_Patient_Names = 0x00101001,
 
 /* Patient's Birth Name $B45<T$NCB@8L>(B */
 Patients_Maiden_Name = 0x00101005,
 
 /* Patient's Age $B45<T$NG/Np(B */
 Patients_Age = 0x00101010,
 
 /* Patient's Size $B45<T$N?HD9(B */
 Patients_Size = 0x00101020,
 
 /* Patient's Weight $B45<T$NBN=E(B */
 Patients_Weight = 0x00101030,
 
 /* Patient's Address $B45<T$N=;=j(B */
 Patients_Address = 0x00101040,
 
 /* Insurance Plan Identification  */
 Insurance_Plan_Identification = 0x00101050,
 
 /* Patient's Mother's Birth Name $B45<T$NJl$NCB@8L>(B */
 Patients_Mothers_Maiden_Name = 0x00101060,
 
 /* Military Rank $B73$N3,5i(B */
 Military_Rank = 0x00101080,
 
 /* Branch of Service $B%5!<%S%9ItLg(B */
 Branch_of_Service = 0x00101081,
 
 /* Medical Record Locator $B0eNE5-O?=j:_<1JL;R(B */
 Medical_Record_Locator = 0x00101090,
 
 /* Medical Alerts $B0eNE7Y9p(B */
 Medical_Alerts = 0x00102000,
 
 /* Contrast Allergies $BB$1F:^%"%l%k%.!<(B */
 Contrast_Allergies = 0x00102110,
 
 /* Country of Residence $B5o=;$N9q(B */
 Country_of_Residence = 0x00102150,
 
 /* Region of Residence $B5o=;$NCO0h(B */
 Region_of_Residence = 0x00102152, 
 
 /* Patient's Telephone Numbers $B45<T$NEEOCHV9f(B */
 Patients_Telephone_Numbers = 0x00102154,
 
 /* Ethnic Group $BL1B2%0%k!<%W(B */
 Ethnic_Group = 0x00102160,
 
 /* Occupation $B?&6H(B */
 Occupation = 0x00102180,
 
 /* Smoking Status $B5J1l$N>uBV(B */
 Smoking_Status = 0x001021A0,
 
 /* Additional Patient History $B45<T$NDI2CIBNr(B */
 Additional_Patient_History = 0x001021B0, 
 
 /* Pregnancy Status $BG%?1$N>uBV(B */
 Pregnancy_Status = 0x001021C0, 
 
 /* Last Menstrual Date $B:G=*7n7PF|(B */
 Last_Menstrual_Date = 0x001021D0,
 
 /* Patient's Religious Preference $B45<T$N?.$8$k=!65(B */
 Patients_Religious_Preference = 0x001021F0,
 
 /* Patient Comments $B45<T%3%a%s%H(B */
 Patient_Comments = 0x00104000,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0018_Length = 0x00180000,
 
 /* Contrast/Bolus Agent $BB$1F:^!?%\!<%i%9Lt:^(B */
 Contrast_Bolus_Agent = 0x00180010,
 
 /* Body Part Examined $B8!::$5$l$kIt0L(B */
 Body_Part_Examined = 0x00180015,
 
 /* Scanning Sequence $B%9%-%c%K%s%0%7!<%1%s%9(B */
 Scanning_Sequence = 0x00180020,
 
 /* Sequence Variant $B%7!<%1%s%9JQ7A(B */
 Sequence_Variant = 0x00180021,
 
 /* Scan Options $B%9%-%c%s%*%W%7%g%s(B */
 Scan_Options = 0x00180022,
 
 /* MR Acquisition Type $B#M#R<}=8%?%$%W(B */
 MR_Acquisition_Type = 0x00180023, 
 
 /* Sequence Name $B%7!<%1%s%9L>(B */
 Sequence_Name = 0x00180024,
 
 /* Angio Flag $B%"%s%8%*%U%i%0(B */
 Angio_Flag = 0x00180025,
 
 /* Radionuclide $BJ|<M@~3K<o(B */
 Radionuclide = 0x00180030,
 
 /* Radio-pharmaceutical $BJ|<M@-0eLtIJ(B */
 Radiopharmaceutical = 0x00180031,
 
 /* Energy Window Centerline $B%(%M%k%.%&%$%s%I%&Cf?4@~(B */
 Energy_Window_Centerline = 0x00180032,
 
 /* Energy Window Total Width $B%(%M%k%.%&%$%s%I%&A4I}(B */
 Energy_Window_Total_Width = 0x00180033,
 
 /* Intervention Drug Name $B%$%s%?%Y%s%7%g%sLtIJL>(B */
 Intervention_Drug_Name = 0x00180034,
 
 /* Intervention Drug Start Time $B%$%s%?%Y%s%7%g%sLtIJ3+;O;~9o(B */
 Intervention_Drug_Start_Time = 0x00180035,
 
 /* Cine Rate $B%7%M$N%U%l!<%`N((B */
 Cine_Rate = 0x00180040, 
 
 /* Slice Thickness $B%9%i%$%98|$5(B */
 Slice_Thickness = 0x00180050,
 
 /* KVP $B#K#V#P(B */
 KVP = 0x00180060,
 
 /* Counts Accumulated $B@Q;;%+%&%s%H(B */
 Counts_Accumulated = 0x00180070, 
 
 /* Acquisition Termination Condition $B<}=8=*N;>r7o(B */
 Acquisition_Termination_Condition = 0x00180071,
 
 /* Effective Series Duration $B<B8z%7%j!<%:;}B3;~4V(B */
 Effective_Series_Duration = 0x00180072,
 
 /* Echo Time $B%(%3!<;~4V(B */
 Echo_Time = 0x00180081, 
 
 /* Inversion Time $BH?E>;~4V(B */
 Inversion_Time = 0x00180082,
 
 /* Number of Averages $BJ?6Q2=2s?t(B */
 Number_of_Averages = 0x00180083,
 
 /* Imaging Frequency $B2hA|<~GH?t(B */
 Imaging_Frequency = 0x00180084,
 
 /* Imaged Nucleus $B2hA|86;R3K(B */
 Imaged_Nucleus = 0x00180085, 
 
 /* Echo Number(s) $B%(%3!<HV9f(B */
 Echo_Numbers = 0x00180086, 
 
 /* Magnetic Field Strength $B<'>l6/EY(B */
 Magnetic_Field_Strength = 0x00180087,
 
 /* Spacing Between Slices $B%9%i%$%94V4V3V(B */
 Spacing_Between_Slices = 0x00180088,
 
 /* Number of Phase Encoding Steps $B0LAjId9f2=%9%F%C%W$N?t(B */
 Number_of_Phase_Encoding_Steps = 0x00180089,
 
 /* Data Collection Diameter $B%G!<%?:N=8D>7B(B */
 Data_Collection_Diameter = 0x00180090,
 
 /* Echo Train Length $B%(%3!<%H%l%sD9$5(B */
 Echo_Train_Length = 0x00180091,
 
 /* Percent Sampling $B%5%s%W%j%s%0I4J,N((B */
 Percent_Sampling = 0x00180093, 
 
 /* Percent Phase Field of View $B0LAj;kLnI4J,N((B */
 Percent_Phase_Field_of_View = 0x00180094,
 
 /* Pixel Bandwidth $B2hAG%P%s%II}(B */
 Pixel_Bandwidth = 0x00180095,
 
 /* Device Serial Number $BAuCV%7%j%"%kHV9f(B */
 Device_Serial_Number = 0x00181000,
 
 /* Plate ID $B%W%l!<%H(BID */
 Plate_ID = 0x00181004, 
 
 /* Secondary Capture Device ID $BFs<!<hF@AuCV(BID */
 Secondary_Capture_Device_ID = 0x00181010,
 
 /* Date of Secondary Capture $BFs<!<hF@$NF|IU(B */
 Date_of_Secondary_Capture = 0x00181012,
 
 /* Time of Secondary Capture $BFs<!<hF@$N;~9o(B */
 Time_of_Secondary_Capture = 0x00181014,
 
 /* Secondary Capture Device Manufacturer $BFs<!<hF@AuCV@=B$<T(B */
 Secondary_Capture_Device_Manufacturer = 0x00181016,
 
 /* Secondary Capture Device Manufacturer's Model Name */
 /* $BFs<!<hF@AuCV@=B$<T7?<0L>(B */
 Secondary_Capture_Device_Manufacturers_Model_Name = 0x00181018,
 
 /* Secondary Capture Device Software Version(s)*/
 /* $BFs<!<hF@AuCV%=%U%H%&%'%"HG(B */
 Secondary_Capture_Device_Software_Version = 0x00181019,
 
 /* Software Version(s) $B%=%U%H%&%'%"HG(B */
 Software_Versions = 0x00181020,
 
 /* Video Image Format Acquired $B<hF@$7$?%S%G%*2hA|$N7A<0(B */
 Video_Image_Format_Acquired = 0x00181022,
 
 /* Digital Image Format Acquired $B<hF@$7$?%G%8%?%k2hA|$N7A<0(B */
 Digital_Image_Format_Acquired = 0x00181023,
 
 /* Protocol Name $B%W%m%H%3%kL>(B */
 Protocol_Name = 0x00181030,
 
 /* Contrast/Bolus Route $BB$1F:^!?%\!<%i%97PO)(B */
 Contrast_Bolus_Route = 0x00181040,
 
 /* Contrast/Bolus Volume $BB$1F:^!?%\!<%i%9MFNL(B */
 Contrast_Bolus_Volume = 0x00181041, 
 
 /* Contrast/Bolus Start Time $BB$1F:^!?%\!<%i%93+;O;~9o(B */
 Contrast_Bolus_Start_Time = 0x00181042,
 
 /* Contrast/Bolus Stop Time $BB$1F:^!?%\!<%i%9Dd;_;~9o(B */
 Contrast_Bolus_Stop_Time = 0x00181043,
 
 /* Contrast/Bolus Total Dose $BB$1F:^!?%\!<%i%9A4EjM?NL(B */
 Contrast_Bolus_Total_Dose = 0x00181044,
 
 /* Syringe counts $BCmF~%+%&%s%H(B */
 Syringe_Counts = 0x00181045,
 
 /* Spatial Resolution $B6u4VJ,2rG=(B */
 Spatial_Resolution = 0x00181050, 
 
 /* Trigger Time $B%H%j%,;~4V(B */
 Trigger_Time = 0x00181060,
 
 /* Trigger Source or Type $B%H%j%,8;$^$?$O%?%$%W(B */
 Trigger_Source_or_Type = 0x00181061,
 
 /* Nominal Interval $B8x>N#R!]#R4V3V(B */
 Nominal_Interval = 0x00181062,
 
 /* Frame Time $B%U%l!<%`;~4V(B */
 Frame_Time = 0x00181063,
 
 /* Framing Type $B%U%l!<%_%s%0%?%$%W(B */
 Framing_Type = 0x00181064,
 
 /* Frame Time Vector $B%U%l!<%`;~4V%Y%/%H%k(B */
 Frame_Time_Vector = 0x00181065,
 
 /* Frame Delay $B%U%l!<%`CY$l(B */
 Frame_Delay = 0x00181066, 
 
 /* Radionuclide Route $BJ|<M@-3K<o7PO)(B */
 Radionuclide_Route = 0x00181070, 
 
 /* Radionuclide Volume $BJ|<M@-3K<oMFNL(B */
 Radionuclide_Volume = 0x00181071,
 
 /* Radionuclide Start Time $BJ|<M@-3K<o3+;O;~4V(B */
 Radionuclide_Start_Time = 0x00181072,
 
 /* Radionuclide Stop Time $BJ|<M@-3K<oDd;_;~4V(B */
 Radionuclide_Stop_Time = 0x00181073,
 
 /* Radionuclide Total Dose $BJ|<M@-3K<oAmEjM?NL(B */
 Radionuclide_Total_Dose = 0x00181074, 
 
 /* Beat Rejection Flag $BGoF0=|5n%U%i%0(B */
 Beat_Rejection_Flag = 0x00181080, 
 
 /* Low R-R Value $B2<8B#R!]#RCM(B */
 Low_R_R_Value = 0x00181081,
 
 /* High R-R Value $B>e8B#R!]#RCM(B */
 High_R_R_Value = 0x00181082, 
 
 /* Intervals Acquired $B<hF@4V3V(B */
 Intervals_Acquired = 0x00181083,
 
 /* Intervals Rejected $B=|5n4V3V(B */
 Intervals_Rejected = 0x00181084, 
 
 /* PVC Rejection $B#P#V#C$N=|5n(B */
 PVC_Rejection = 0x00181085, 
 
 /* Skip Beats $B%9%-%C%WGoF0(B */
 Skip_Beats = 0x00181086,
 
 /* Heart Rate $B?4Go?t(B */
 Heart_Rate = 0x00181088,
 
 /* Cardiac Number of Images $B?42hA|$N?t(B */
 Cardiac_Number_of_Images = 0x00181090,
 
 /* Trigger Window $B%H%j%,%&%$%s%I%&(B */
 Trigger_Window = 0x00181094,
 
 /* Reconstruction Diameter $B:F9=@.D>7B(B */
 Reconstruction_Diameter = 0x00181100,
 
 /* Distance Source to Detector $B@~8;8!=P4o4V5wN%(B */
 Distance_Source_to_Detector = 0x00181110, 
 
 /* Distance Source to Patient $B@~8;45<T4V5wN%(B */
 Distance_Source_to_Patient = 0x00181111,
 
 /* Gantry/Detector Tilt $B2MBf!?8!=P4o79$-(B */
 Gantry_Detector_Tilt = 0x00181120,
 
 /* Table Height $B%F!<%V%k9b$5(B */
 Table_Height = 0x00181130,
 
 /* Table Traverse $B%F!<%V%k%H%i%P!<%9(B */
 Table_Traverse = 0x00181131,
 
 /* Rotation Direction $B2sE>J}8~(B */
 Rotation_Direction = 0x00181140,
 
 /* Angular Position $B3QEY0LCV(B */
 Angular_Position = 0x00181141,
 
 /* Radial Position $B7BJ}8~0LCV(B */
 Radial_Position = 0x00181142,
 
 /* Scan Arc $B%9%-%c%s%"!<%/(B */
 Scan_Arc = 0x00181143,
 
 /* Angular Step $B3QEY%9%F%C%W(B */
 Angular_Step = 0x00181144,
 
 /* Center of Rotation Offset $B2sE>%*%U%;%C%H$NCf?4(B */
 Center_of_Rotation_Offset = 0x00181145, 
 
 /* Rotation Offset $B2sE>%*%U%;%C%H(B */
 Rotation_Offset = 0x00181146,
 
 /* Field of View Shape $B;kLn$N7A>u(B */
 Field_of_View_Shape = 0x00181147,
 
 /* Field of View Dimension(s) $B;kLn$N@#K!(B */
 Field_of_View_Dimensions = 0x00181149, 
 
 /* Exposure Time $BGx<M;~4V(B */
 Exposure_Time = 0x00181150,
 
 /* X-ray Tube Current $B#X@~4IEEN.(B */
 X_ray_Tube_Current = 0x00181151,
 
 /* Exposure $BGx<MNL(B */
 Exposure = 0x00181152, 
 
 /* Filter Type $B%U%#%k%?%?%$%W(B */
 Filter_Type = 0x00181160,
 
 /* Generator Power $BH/@8AuCV=PNO(B */
 Generator_Power = 0x00181170,
 
 /* Collimator/grid Name $B%3%j%a!<%?!?%0%j%C%IL>(B */
 Collimator_grid_Name = 0x00181180, 
 
 /* Collimator Type $B%3%j%a!<%?%?%$%W(B */
 Collimator_Type = 0x00181181,
 
 /* Focal Distance $B>GE@5wN%(B */
 Focal_Distance = 0x00181182,
 
 /* X Focus Center X $B>GE@Cf?4(B */
 X_Focus_Center = 0x00181183,
 
 /* Y Focus Center Y $B>GE@Cf?4(B */
 Y_Focus_Center = 0x00181184,
 
 /* Focal Spot(s) $B>GE@(B */
 Focal_Spot = 0x00181190, 
 
 /* Date of Last Calibration $B:G=*3S@5F|IU(B */
 Date_of_Last_Calibration = 0x00181200,
 
 /* Time of Last Calibration $B:G=*3S@5;~9o(B */
 Time_of_Last_Calibration = 0x00181201, 
 
 /* Convolution Kernel $B%3%s%\%j%e!<%7%g%s%+!<%M%k(B */
 Convolution_Kernel = 0x00181210, 
 
 /* Upper/Lower Pixel Values  */
 Upper_Lower_Pixel_Values = 0x00181240,
 
 /* Actual Frame Duration $B<B%U%l!<%`;}B3;~4V(B */
 Actual_Frame_Duration = 0x00181242, 
 
 /* Count Rate $B%+%&%s%HN((B */
 Count_Rate = 0x00181243,
 
 /* Receiving Coil $B<u?.%3%$%k(B */
 Receiving_Coil = 0x00181250,
 
 /* X-ray Tube Current $B#X@~4IEEN.(B */
 Transmitting_Coil = 0x00181151,
 
 /* Filter Type $B%U%#%k%?%?%$%W(B */
 Screen_Type = 0x00181160, 
 
 /* Phosphor Type $B7V8wBN%?%$%W(B */
 Phosphor_Type = 0x00181261,
 
 /* Scan Velocity $B%9%-%c%sB.EY(B */
 Scan_Velocity = 0x00181300, 
 
 /* Whole Body Technique $BA4?H5;=Q(B */
 Whole_Body_Technique = 0x00181301,
 
 /* Scan Length $B%9%-%c%sD9$5(B */
 Scan_Length = 0x00181302,
 
 /* Acquisition Matrix $B<}=8%^%H%j%C%/%9(B */
 Acquisition_Matrix = 0x00181310,
 
 /* Phase Encoding Direction $B0LAjId9f2=J}8~(B */
 Phase_Encoding_Direction = 0x00181312,
 
 /* Flip Angle $B%U%j%C%W3Q(B */
 Flip_Angle = 0x00181314,
 
 /* Variable Flip Angle Flag $B2DJQ%U%j%C%W3Q%U%i%0(B */
 Variable_Flip_Angle_Flag = 0x00181315,
 
 /* SAR */
 SAR = 0x00181316, 
 
 /* dB/dt dB/dt */
 dB_dt = 0x00181318,
 
 /* Acquisition Device Processing Description $B<}=8AuCV=hM}5-=R(B */
 Acquisition_Device_Processing_Description = 0x00181400, 
 
 /* Acquisition Device Processing Code $B<}=8AuCV=hM}%3!<%I(B */
 Acquisition_Device_Processing_Code = 0x00181401, 
 
 /* Cassette Orientation $B%+%;%C%FJ}0L(B */
 Cassette_Orientation = 0x00181402, 
 
 /* Cassette Size $B%+%;%C%F%5%$%:(B */
 Cassette_Size = 0x00181403,
 
 /* Exposures on Plate $B%W%l!<%H>e$NGx<MNL(B */
 Exposures_on_Plate = 0x00181404, 
 
 /* Relative X-ray Exposure $BAjBP#X@~Gx<MNL(B */
 Relative_X_ray_Exposure = 0x00181405,
 
 /* Comments  */
 Group_0018_Comments = 0x00184000,
 
 /* Output Power $B=PNO(B */
 Output_Power = 0x00185000,
 
 /* Transducer Data $BC5?(;R%G!<%?(B */
 Transducer_Data = 0x00185010,
 
 /* Focus Depth $B>GE@?<$5(B */
 Focus_Depth = 0x00185012,
 
 /* Preprocessing Function $BA0=hM}4X?t(B */
 Preprocessing_Function = 0x00185020,
 
 /* Postprocessing Function $B8e=hM}4X?t(B */
 Postprocessing_Function = 0x00185021,
 
 /* Mechanical Index MI */
 Mechanical_Index = 0x00185022,
 
 /* Thermal Index $B9|(BTI */
 Thermal_Index = 0x00185024, 
 
 /* Cranial Thermal Index $BF,38(BTI */
 Cranial_Thermal_Index = 0x00185026, 
 
 /* Soft Tissue Thermal Index $BFpAH?%(BTI */
 Soft_Tissue_Thermal_Index = 0x00185027, 
 
 /* Soft Tissue-focus Thermal Index $BFpAH?%>GE@(BTI */
 Soft_Tissue_focus_Thermal_Index = 0x00185028,
 
 /* Soft Tissue-surface Thermal Index $BFpAH?%I=LL(BTI */
 Soft_Tissue_surface_Thermal_Index = 0x00185029,
 
 /* Dynamic Range  */
 Dynamic_Range = 0x00185030, 
 
 /* Total Gain  */
 Total_Gain = 0x00185040, 
 
 /* Depth of Scan Field $BI=<($N?<$5(B */
 Depth_of_Scan_Field = 0x00185050,
 
 /* Patient Position $B45<T0LCV(B */
 Patient_Position = 0x00185100,
 
 /* View Position $B;kLn0LCV(B */
 View_Position = 0x00185101,
 
 /* Image Transformation Matrix $B2hA|JQ49%^%H%j%C%/%9(B */
 Image_Transformation_Matrix = 0x00185210,
 
 /* Image Translation Vector $B2hA|JQ49%Y%/%H%k(B */
 Image_Translation_Vector = 0x00185212, 
 
 /* Sensitivity $B46EY(B */
 Sensitivity = 0x00186000, 
 
 /* Sequence of Ultrasound Regions $BD62;GHNN0h%7!<%1%s%9(B */
 Sequence_of_Ultrasound_Regions = 0x00186011, 
 
 /* Region Spatial Format $BNN0h6u4V%U%)!<%^%C%H(B */
 Region_Spatial_Format = 0x00186012, 
 
 /* Region Data Type $BNN0h$N%G!<%?%?%$%W(B */
 Region_Data_Type = 0x00186014,
 
 /* Region Flags $BNN0h%U%i%0(B */
 Region_Flags = 0x00186016, 
 
 /* Region Location Min X0 $BNN0h0LCV(B Min x0 */
 Region_Location_Min_X0 = 0x00186018,
 
 /* Region Location Min Y0 $BNN0h0LCV(B Min y0 */
 Region_Location_Min_Y0 = 0x0018601A, 
 
 /* Region Location Max X1 $BNN0h0LCV(B Max x1 */
 Region_Location_Max_X1 = 0x0018601C,
 
 /* Region Location Max Y1 $BNN0h0LCV(B Max y1 */
 Region_Location_Max_Y1 = 0x0018601E,
 
 /* Reference Pixel X0 $B;2>H2hAG(B x0 */
 Reference_Pixel_X0 = 0x00186020, 
 
 /* Reference Pixel Y0 $B;2>H2hAG(B y0 */
 Reference_Pixel_Y0 = 0x00186022, 
 
 /* Physical Units X Direction $BJ*M}C10L(B X $BJ}8~(B */
 Physical_Units_X_Direction = 0x00186024,
 
 /* Physical Units Y Direction $BJ*M}C10L(B Y $BJ}8~(B */
 Physical_Units_Y_Direction = 0x00186026,
 
 /* Reference Pixel Physical Value X $B;2>H2hAGJ*M}CM(B X */
 Reference_Pixel_Physical_Value_X = 0x00186028, 
 
 /* Reference Pixel Physical Value Y $B;2>H2hAGJ*M}CM(B Y */
 Reference_Pixel_Physical_Value_Y = 0x0018602A,
 
 /* Physical Delta X $BJ*M}JQ2=NL(B X */
 Physical_Delta_X = 0x0018602C,
 
 /* Physical Delta Y $BJ*M}JQ2=NL(B Y */
 Physical_Delta_Y = 0x0018602E,
 
 /* Transducer Frequency $BC5?(;R<~GH?t(B */
 Transducer_Frequency = 0x00186030,
 
 /* Transducer Type $BC5?(;R%?%$%W(B */
 Transducer_Type = 0x00186031,
 
 /* Pulse Repetition Frequency $B%Q%k%97+JV$7<~GH?t(B */
 Pulse_Repetition_Frequency = 0x00186032, 
 
 /* Doppler Correction Angle $B%I%W%i3QEYJd@5(B */
 Doppler_Correction_Angle = 0x00186034, 
 
 /* Steering Angle $B%9%F%"%j%s%03QEY(B */
 Sterring_Angle = 0x00186036,
 
 /* Doppler Sample Volume X Position $B%I%W%i%5%s%W%kMF@Q(B X $B0LCV(B */
 Doppler_Sample_Volume_X_Position = 0x00186038, 
 
 /* Doppler Sample Volume Y Position $B%I%W%i%5%s%W%kMF@Q(B Y $B0LCV(B */
 Doppler_Sample_Volume_Y_Position = 0x0018603A,
 
 /* TM-Line Position X0 TM$B@~$N0LCV(B x0 */
 TM_Line_Position_X0 = 0x0018603C,
 
 /* TM-Line Position Y0 TM$B@~$N0LCV(B y0 */
 TM_Line_Position_Y0 = 0x0018603E,
 
 /* TM-Line Position X1 TM$B@~$N0LCV(B x1 */
 TM_Line_Position_X1 = 0x00186040, 
 
 /* TM-Line Position Y1 TM$B@~$N0LCV(B y1 */
 TM_Line_Position_Y1 = 0x00186042,
 
 /* Pixel Component Organization $B2hAG9=@.MWAG$NJ}<0(B */
 Pixel_Component_Organization = 0x00186044,
 
 /* Pixel Component Mask $B2hAG9=@.MWAG%^%9%/(B */
 Pixel_Component_Mask = 0x00186046,
 
 /* Pixel Component Range Start $B2hAG9=@.MWAGHO0O$N;OE@(B */
 Pixel_Component_Range_Start = 0x00186048,
 
 /* Pixel Component Range Stop $B2hAG9=@.MWAGHO0O$N=*E@(B */
 Pixel_Component_Range_Stop = 0x0018604A, 
 
 /* Pixel Component Physical Units $B2hAG9=@.MWAGJ*M}C10L(B */
 Pixel_Component_Physical_Units = 0x0018604C,
 
 /* Pixel Component Data Type $B2hAG9=@.MWAG%G!<%?%?%$%W(B */
 Pixel_Component_Data_Type = 0x0018604E, 
 
 /* Number of Table Break Points $BI=%V%l%$%/%]%$%s%H$N?t(B */
 Number_of_Table_Break_Points = 0x00186050,
 
 /* Table of X Break Points X $B%V%l%$%/%]%$%s%H$NI=(B */
 Table_of_X_Break_Points = 0x00186052,
 
 /* Table of Y Break Points Y $B%V%l%$%/%]%$%s%H$NI=(B */
 Table_of_Y_Break_Points = 0x00186054,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0020_Length = 0x00200000, 
 
 /* Study Instance UID $B8!::%$%s%9%?%s%9(BUID */
 Study_Instance_UID = 0x0020000D, 
 
 /* Series Instance UID $B%7%j!<%:%$%s%9%?%s%9(BUID */
 Series_Instance_UID = 0x0020000E,
 
 /* Study ID $B8!::(BID */
 Study_ID = 0x00200010, 
 
 /* Series Number $B%7%j!<%:HV9f(B */
 Series_Number = 0x00200011, 
 
 /* Acquisition Number $B<}=8HV9f(B */
 Acquisition_Number = 0x00200012,
 
 /* Image Number $B2hA|HV9f(B */
 Image_Number = 0x00200013,
 
 /* Isotope Number $BF10L85AGHV9f(B */
 Isotope_Number = 0x00200014, 
 
 /* Phase Number $B0LAjHV9f(B */
 Phase_Number = 0x00200015,  
 
 /* Interval Number $B4V3VHV9f(B */
 Interval_Number = 0x00200016,
 
 /* Time Slot Number $B;~4V%9%m%C%HHV9f(B */
 Time_Slot_Number = 0x00200017,  
 
 /* Angle Number $B3QEYHV9f(B */
 Angle_Number = 0x00200018,  
 
 /* Patient Orientation $B45<TJ}8~(B */
 Patient_Orientation = 0x00200020,  
 
 /* Overlay Number $B%*!<%P%l%$HV9f(B */
 Overlay_Number = 0x00200022,
 
 /* Curve Number $B%+!<%VHV9f(B */
 Curve_Number = 0x00200024,
 
 /* Image Position  */
 Image_Position_0020 = 0x00200030,
 
 /* Image Position (Patient) $B2hA|0LCV!J45<T!K(B */
 Image_Position_Patient = 0x00200032,
 
 /* Image Orientation  */
 Image_Orientation = 0x00200035,
 
 /* Image Orientation (Patient) $B2hA|J}8~!J45<T!K(B */
 Image_Orientation_Patient = 0x00200037,
 
 /* Location  */
 Location = 0x00200050,
 
 /* Frame of Reference UID $B4p=`:BI87O(BUID */
 Frame_of_Reference_UID = 0x00200052,
 
 /* Laterality $B:81&(B */
 Laterality = 0x00200060,
 
 /* Image Geometry Type  */
 Image_Geometry_Type = 0x00200070,
 
 /* Masking Image  */
 Masking_Image_UID = 0x00200080,
 
 /* Temporal Position Identifier $B;~4V0LCV<1JL;R(B */
 Temporal_Position_Identifier = 0x00200100, 
 
 /* Number of Temporal Positions $B;~4VE*0LCV$N?t(B */
 Number_of_Temporal_Positions = 0x00200105,
 
 /* Temporal Resolution $B;~4VJ,2rG=(B */
 Temporal_Resolution = 0x00200110,
 
 /* Series in Study $B8!::$NCf$N%7%j!<%:(B */
 Series_in_Study = 0x00201000,
 
 /* Acquisitions in Series  */
 Acquisitions_in_Series = 0x00201001,
 
 /* Images in Acquisition $B<}=8$NCf$N2hA|(B */
 Images_in_Acquisition = 0x00201002,
 
 /* Acquisition in Study $B8!::$NCf$N<}=8(B */
 Acquisition_in_Study = 0x00201004,
 
 /* Reference  */
 Reference = 0x00201020, 
 
 /* Position Reference Indicator $B0LCV;2>H%$%s%8%1!<%?(B */
 Position_Reference_Indicator = 0x00201040, 
 
 /* Slice Location $B%9%i%$%90LCV(B */
 Slice_Location = 0x00201041, 
 
 /* Other Study Numbers $BB>$N8!::HV9f(B */
 Other_Study_Numbers = 0x00201070,
 
 /* Number of Patient Related Studies $B45<T$K4X78$7$?8!::$N?t(B */
 Number_of_Patient_Related_Studies = 0x00201200,
 
 /* Number of Patient Related Series $B45<T$K4X78$7$?%7%j!<%:$N?t(B */
 Number_of_Patient_Related_Series = 0x00201202,
 
 /* Number of Patient Related Images $B45<T$K4X78$7$?2hA|$N?t(B */
 Number_of_Patient_Related_Images = 0x00201204,
 
 /* Number of Study Related Series $B8!::$K4X78$7$?%7%j!<%:$N?t(B */
 Number_of_Study_Related_Series = 0x00201206, 
 
 /* Number of Study Related Images $B8!::$K4X78$7$?2hA|$N?t(B */
 Number_of_Study_Related_Images = 0x00201208,
 
 /* Source Image IDs  */
 Source_Image_IDs = 0x00203100,

 /* Modifying Device ID  */
 Modifying_Device_ID = 0x00203401, 
 
 /* Modified Image ID  */
 Modified_Image_ID = 0x00203402, 
 
 /* Modified Image Date  */
 Modified_Image_Date = 0x00203403,
 
 /* Modifying Device Manufacturer  */
 Modifying_Device_Manufacturer = 0x00203404,
 
 /* Modified Image Time  */
 Modified_Image_Time = 0x00203405,
 
 /* Modified Image Description  */
 Modified_Image_Description = 0x00203406,
 
 /* Image Comments $B2hA|%3%a%s%H(B */
 Image_Comments = 0x00204000,
 
 /* Original Image Identification  */
 Original_Image_Identification = 0x00205000,
 
 /* Original Image Identification Nomenclature  */
 Original_Image_Identification_Nomenclature = 0x00205002, 
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0028_Length = 0x00280000,
 
 /* Samples per Pixel $B2hAG$"$?$j%5%s%W%k(B */
 Samples_per_Pixel = 0x00280002, 
 
 /* Photometric Interpretation $B8wEYB,Dj2r<a(B */
 Photometric_Interpretation = 0x00280004, 
 
 /* Image Dimensions  */
 Image_Dimensions = 0x00280005,
 
 /* Planar Configuration $BLL9=@.(B */
 Planar_Configuration = 0x00280006, 
 
 /* Number of Frames $B%U%l!<%`$N?t(B */
 Number_of_Frames = 0x00280008,
 
 /* Frame Increment Pointer $B%U%l!<%`A}J,%]%$%s%?(B */
 Frame_Increment_Pointer = 0x00280009,
 
 /* Rows $B9T(B */
 Rows = 0x00280010,
 
 /* Columns $BNs(B */
 Columns = 0x00280011,
 
 /* Pixel Spacing $B2hAG4V3V(B */
 Pixel_Spacing = 0x00280030,
 
 /* Zoom Factor $B%:!<%`N((B */
 Zoom_Factor = 0x00280031,
 
 /* Zoom Center $B%:!<%`Cf?4(B */
 Zoom_Center = 0x00280032,
 
 /* Pixel Aspect Ratio $B2hAG%"%9%Z%/%HHf(B */
 Pixel_Aspect_Ratio = 0x00280034, 
 
 /* Image Format  */
 Image_Format = 0x00280040,
 
 /* Manipulated Image  */
 Manipulated_Image = 0x00280050,
 
 /* Corrected Image $BJd@52hA|(B */
 Corrected_Image = 0x00280051,
 
 /* Compression Code  */
 Compression_Code_0028 = 0x00280060,
 
 /* Bits Allocated $B3d$jEv$F%S%C%H(B */
 Bits_Allocated = 0x00280100,
 
 /* Bits Stored $B3JG<%S%C%H(B */
 Bits_Stored = 0x00280101,
 
 /* High Bit $B9b0L%S%C%H(B */
 High_Bit = 0x00280102,
 
 /* Pixel Representation $B2hAGI=8=(B */
 Pixel_Representation = 0x00280103,
 
 /* Smallest Valid Pixel Value  */
 Smallest_Valid_Pixel_Value = 0x00280104, 
 
 /* Largest Valid Pixel Value  */
 Largest_Valid_Pixel_Value = 0x00280105,
 
 /* Smallest Image Pixel Value $B:G>.2hA|2hAGCM(B */
 Smallest_Image_Pixel_Value = 0x00280106,
 
 /* Largest Image Pixel Value $B:GBg2hA|2hAGCM(B */
 Largest_Image_Pixel_Value = 0x00280107, 
 
 /* Smallest Pixel Value in Series $B%7%j!<%:$NCf$N:G>.2hAGCM(B */
 Smallest_Pixel_Value_in_Series = 0x00280108, 
 
 /* Largest Pixel Value in Series $B%7%j!<%:$NCf$N:GBg2hAGCM(B */
 Largest_Pixel_Value_in_Series = 0x00280109, 
 
 /* Pixel Padding Value $B2hAG%Q%G%#%s%0CM(B */
 Pixel_Padding_Value = 0x00280120,
 
 /* Image Location  */
 Image_Location = 0x00280200,
 
 /* Window Center $B%&%#%s%I%&Cf?4(B */
 Window_Center = 0x00281050, 
 
 /* Window Width $B%&%#%s%I%&I}(B */
 Window_Width = 0x00281051,
 
 /* Rescale Intercept $B%j%9%1!<%k@ZJR(B */
 Rescale_Intercept = 0x00281052, 
 
 /* Rescale Slope $B%j%9%1!<%k79<P(B */
 Rescale_Slope = 0x00281053,
 
 /* Rescale type $B%j%9%1!<%k%?%$%W(B */
 Rescale_Type = 0x00281054,
 
 /*  */
 Window_Center_Width_Explanation = 0x00281055, 
 
 /* Gray Scale  */
 Gray_Scale = 0x00281080,
 
 /* Gray Lookup Table Descriptor  */
 Gray_Lookup_Table_Descriptor = 0x00281100,
 
 /* Red Palette Color Lookup Table Descriptor */
 /* $B%l%C%I%Q%l%C%H%+%i!<%k%C%/%"%C%W%F!<%V%k5-=R;R(B */
 Red_Palette_Color_Lookup_Table_Descriptor = 0x00281101, 
 
 /* Green Palette Color Lookup Table Descriptor */
 /* $B%0%j!<%s%Q%l%C%H%+%i!<%k%C%/%"%C%W%F!<%V%k5-=R;R(B */
 Green_Palette_Color_Lookup_Table_Descriptor = 0x00281102,
 
 /* Blue Palette Color Lookup Table Descriptor */
 /*$B%V%k!<%Q%l%C%H%+%i!<%k%C%/%"%C%W%F!<%V%k5-=R;R(B */
 Blue_Palette_Color_Lookup_Table_Descriptor = 0x00281103,
 
 /* Gray Lookup Table Data */
 Gray_Lookup_Table_Data = 0x00281200,
 
 /* Red Palette Color Lookup Table Data */
 /* $B%l%C%I%Q%l%C%H%+%i!<%k%C%/%"%C%W%F!<%V%k%G!<%?(B */
 Red_Palette_Color_Lookup_Table_Data = 0x00281201,  
 
 /* Green Palette Color Lookup Table Data */
 /* $B%0%j!<%s%Q%l%C%H%+%i!<%k%C%/%"%C%W%F!<%V%k%G!<%?(B */
 Green_Palette_Color_Lookup_Table_Data = 0x00281202,
 
 /* Blue Palette Color Lookup Table Data */
 /* $B%V%k!<%Q%l%C%H%+%i!<%k%C%/%"%C%W%F!<%V%k%G!<%?(B */
 Blue_Palette_Color_Lookup_Table_Data = 0x00281203,
 
 /* Modality LUT Sequence $B%b%@%j%F%#(BLUT$B%7!<%1%s%9(B */
 Modality_LUT_Sequence = 0x00283000, 
 
 /* LUT Descriptor LUT$B5-=R;R(B */
 LUT_Descriptor = 0x00283002, 
 
 /* LUT Explanation LUT$B@bL@(B */
 LUT_Explanation = 0x00283003, 
 
 /* Modality LUT Type $B%b%@%j%F%#(BLUT$B%?%$%W(B */
 Madality_LUT_Type = 0x00283004,
 
 /* LUT Data LUT$B%G!<%?(B */
 LUT_Data = 0x00283006, 
 
 /* VOI LUT Sequence VOI LUT$B%7!<%1%s%9(B */
 VOI_LUT_Sequence = 0x00283010,
 
 /* Comments  */
 Group_0028_Comments = 0x00284000,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0032_Length = 0x00320000,
 
 /* Study Status ID $B8!::>uBV(BID */
 Study_Status_ID = 0x0032000A,
 
 /* Study Priority ID $B8!::M%@hEY(BID */
 Study_Priority_ID = 0x0032000C, 
 
 /* Study ID Issuer $B8!::(BID$BH/9T<T(B */
 Study_ID_Issuer = 0x00320012,
 
 /* Study Verified Date $B8!::3NG'F|IU(B */
 Study_Verified_Date = 0x00320032,  
 
 /* Study Verified Time $B8!::3NG';~9o(B */
 Study_Verified_Time = 0x00320033, 
 
 /* Study Read Date $B8!::FI1FF|IU(B */
 Study_Read_Date = 0x00320034,
 
 /* Study Read Time $B8!::FI1F;~9o(B */
 Study_Read_Time = 0x00320035, 
 
 /* Scheduled Study Start Date $BM=Ls8!::3+;OF|IU(B */
 Scheduled_Study_Start_Date = 0x00321000,
 
 /* Scheduled Study Start Time $BM=Ls8!::3+;O;~9o(B */
 Scheduled_Study_Start_Time = 0x00321001,
 
 /* Scheduled Study Stop Date $BM=Ls8!::=*N;F|IU(B */
 Scheduled_Study_Stop_Date = 0x00321010,
 
 /* Scheduled Study Stop Time $BM=Ls8!::=*N;;~9o(B */
 Scheduled_Study_Stop_Time = 0x00321011,
 
 /* Scheduled Study Location $BM=Ls8!::>l=j(B */
 Scheduled_Study_Location = 0x00321020,
 
 /* Scheduled Study Location AE Title(s) $BM=Ls8!::>l=j(BAE$BL>>N(B */
 Scheduled_Study_Location_AE_Title = 0x00321021,
 
 /* Reason for Study $B8!::$K$D$$$F$NM}M3(B */
 Reason__for_Study = 0x00321030,
 
 /* Requesting Physician $B0MMjB&0e;U(B */
 Requesting_Physician = 0x00321032, 
 
 /* Requesting Service $B0MMjB&%5!<%S%9(B */
 Requesting_Service = 0x00321033,
 
 /* Study Arrival Date $B8!::E~CeF|IU(B */
 Study_Arrival_Date = 0x00321040,
 
 /* Study Arrival Time $B8!::E~Ce;~9o(B */
 Study_Arrival_Time = 0x00321041,
 
 /* Study Completion Date $B8!::40N;F|IU(B */
 Study_Completion_Date = 0x00321050,
 
 /* Study Completion Time $B8!::40N;;~9o(B */
 Study_Completion_Time = 0x00321051,
 
 /* Study Component Status ID $B8!::9=@.MWAG>uBV(BID */
 Study_Component_Status_ID = 0x00321055,
 
 /* Requested Procedure Description $B0MMj=hCV5-=R(B */
 Requested_Procedure_Description = 0x00321060,
 
 /* Requested Procedure Code Sequence $B0MMj=hCV%3!<%I%7!<%1%s%9(B */
 Requested_Procedure_Code_Sequence = 0x00321064,
 
 /* Requested Contrast Agent $B0MMjB$1F:^(B */
 Requested_Contrast_Agent = 0x00321070,
 
 /* Study Comments $B8!::%3%a%s%H(B */
 Study_Comments = 0x00324000, 
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0038_Length = 0x00380000, 
 
 /* Referenced Patient Alias Sequence $B;2>H45<TJLL>%7!<%1%s%9(B */
 Referenced_Patient_Alias_Sequence = 0x00380004, 
 
 /* Visit Status ID $BMh1!>uBV(BID */
 Visit_Status_ID = 0x00380008,
 
 /* Admission ID $B<u?G(BID */
 Admissin_ID = 0x00380010, 
 
 /* Issuer of Admission ID $B<u?G(BID$B$NH/9T<T(B */
 Issuer_of_Admission_ID = 0x00380011,
 
 /* Route of Admissions $B<u?G$N7PO)(B */
 Route_of_Admissions = 0x00380016,
 
 /* Scheduled Admission Date $BM=Ls<u?GF|IU(B */
 Scheduled_Admissin_Date = 0x0038001A,  
 
 /* Scheduled Admission Time $BM=Ls<u?G;~9o(B */
 Scheduled_Adission_Time = 0x0038001B,
 
 /* Scheduled Discharge Date $BM=LsB`1!F|IU(B */
 Scheduled_Discharge_Date = 0x0038001C, 
 
 /* Scheduled Discharge Time $BM=LsB`1!;~4V(B */
 Scheduled_Discharge_Time = 0x0038001D,
 
 /* Scheduled Patient Institution Residence $BM=Ls45<T;\@_Fb5o=;>l=j(B */
 Scheduled_Patient_Institution_Residence = 0x0038001E,
 
 /* Admitting Date $B<u?GF|IU(B */
 Admitting_Date = 0x00380020, 
 
 /* Admitting Time $B<u?G;~9o(B */
 Admitting_Time = 0x00380021, 
 
 /* Discharge Date $BB`1!F|IU(B */
 Discharge_Date = 0x00380030,  
 
 /* Discharge Time $BB`1!;~9o(B */
 Discharge_Time = 0x00380032,  
 
 /* Discharge Diagnosis Description $BB`1!;~?GCG5-=R(B */
 Discharge_Diagnosis_Description = 0x00380040,  
 
 /* Discharge Diagnosis Code Sequence $BB`1!;~?GCG%3!<%I%7!<%1%s%9(B */
 Discharge_Diagnosis_Code_Sequence = 0x00380044,
 
 /* Special Needs $BFCJL$J2p=u(B */
 Special_Needs = 0x00380050,
 
 /* Current Patient Location $B8=:_$N45<T$N=;=j(B */
 Current_Patient_Location = 0x00380300, 
 
 /* Patient's Institution Residence $B45<T$N;\@_Fb5o=;(B */
 Patients_Institution_Residence = 0x00380400,  
 
 /* Patient State $B45<T$N>uBV(B */
 Patient_State = 0x00380500,
 
 /* Visit Comments $BMh1!%3%a%s%H(B */
 Visit_Comments = 0x00384000,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_0088_Length = 0x00880000,
 
 /* Storage Media File-set ID $BJ]B8G^BN%U%!%$%k=89g(BID */
 Storage_Media_File_set_ID = 0x00880130, 
 
 /* Storage Media File-set UID $BJ]B8G^BN%U%!%$%k=89g(BUID */
 Storage_Media_File_set_UID = 0x00880140, 
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_2000_Length = 0x20000000,
 
 /* Number Of Copies $B%3%T!<$N?t(B */
 Number_of_Copies = 0x20000010,
 
 /* Print Priority $B%W%j%s%HM%@hEY(B */
 Print_Priority = 0x20000020,
 
 /* Medium Type $BG^BN%?%$%W(B */
 Medium_Type = 0x20000030,
 
 /* Film Destination $B%U%#%k%`$"$F@h(B */
 Film_Destination = 0x20000040,
 
 /* Film Session Label $B%U%#%k%`%;%C%7%g%s%i%Y%k(B */
 Film_Session_Label = 0x20000050,
 
 /* Memory Allocation $B%a%b%j3d$jEv$F(B */
 Memory_Allocation = 0x20000060,
 
 /* Referenced Film Box Sequence $B;2>H%U%#%k%`%\%C%/%9%7!<%1%s%9(B */
 Referenced_Film_Box_Sequence = 0x20000500,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_2010_Length = 0x20100000,
 
 /* Image Display Format $B2hA|I=<(%U%)!<%^%C%H(B */
 Image_Display_Format = 0x20100010,
 
 /* Annotation Display Format ID $BCm<aI=<(%U%)!<%^%C%H(BID */
 Annotation_Display_Format_ID = 0x20100030,
 
 /* Film Orientation $B%U%#%k%`J}8~(B */
 Film_Orientation = 0x20100040,
 
 /* Film Size ID $B%U%#%k%`%5%$%:(BID */
 Film_Size_ID = 0x20100050,
 
 /* Magnification Type $B3HBg%?%$%W(B */
 Magnification_Type = 0x20100060,
 
 /* Smoothing Type $BJ?3j%?%$%W(B */
 Smoothing_Type = 0x20100080,
 
 /* Border Density $B1o<h$jG;EY(B */
 Border_Density = 0x20100100,
 
 /* Empty Image Density $B6u2hA|G;EY(B */
 Empty_Image_Density = 0x20100110, 
 
 /* Min Density $B:GDcG;EY(B */
 Min_Density = 0x20100120,
 
 /* Max Density $B:G9bG;EY(B */
 Max_Density = 0x20100130, 
 
 /* Trim $B$U$A>~$j(B */
 Trim = 0x20100140, 
 
 /* Configuration Information $B9=@.>pJs(B */
 Configuration_Information = 0x20100150,
 
 /* Referenced Film Session Sequence $B;2>H%U%#%k%`%;%C%7%g%s%7!<%1%s%9(B */
 Referenced_Film_Session_Sequence = 0x20100500,
 
 /* Referenced Image Box Sequence $B;2>H2hA|%\%C%/%9%7!<%1%s%9(B */
 Referenced_Basic_Image_Box_Sequence = 0x20100510, 
 
 /* Referenced Basic Annotation Box Sequence $B;2>H4pK\Cm<a%\%C%/%9%7!<%1%s%9(B */
 Referenced_Basic_Annotation_Box_Sequence = 0x20100520,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_2020_Length = 0x20200000, 
 
 /* Image Position $B2hA|0LCV(B */
 Image_Position = 0x20200010, 
 
 /* Polarity $B6K@-(B */
 Polarity = 0x20200020,  
 
 /* Requested Image Size $B0MMj2hA|@#K!(B */
 Requested_Image_Size = 0x20200030,
 
 /* Preformatted Grayscale Image Sequence */
 /* $B%U%)!<%^%C%H:Q%0%l%$%9%1!<%k2hA|%7!<%1%s%9(B */
 Preformatted_Greyscale_Image_Sequence = 0x20200110, 
 
 /* Preformatted Color Image Sequence*/
 /* $B%U%)!<%^%C%H:Q%+%i!<2hA|%7!<%1%s%9(B */
 Preformatted_Color_Image_Sequence = 0x20200111,
 
 /* Referenced Image Overlay Box Sequence */
 /* $B;2>H2hA|%*!<%P%l%$%\%C%/%9%7!<%1%s%9(B */
 Referenced_Image_Overlay_Box_Sequence = 0x20200130, 
 
 /* Referenced VOI LUT Box Sequence $B;2>H(BVOI LUT$B%\%C%/%9%7!<%1%s%9(B */
 Referenced_VOI_LUT_Sequence = 0x20200140,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_2030_Length = 0x20300000, 
 
 /* Annotation Position $BCm<a0LCV(B */
 Annotation_Position = 0x20300010,
 
 /* Text String $B%F%-%9%HNs(B */
 Text_String = 0x20300020,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_2040_Length = 0x20400000,
 
 /* Referenced Overlay Plane Sequence $B;2>H%*!<%P%l%$LL%7!<%1%s%9(B */
 Referenced_Overlay_Plane_Sequence = 0x20400010, 
 
 /* Referenced Overlay Plane Groups $B;2>H%*!<%P%l%$LL%0%k!<%W(B */
 Refenced_Overlay_Plane_Groups = 0x20400011,
 
 /* Overlay Magnification Type $B%*!<%P%l%$3HBg%?%$%W(B */
 Overlay_Magnification_Type = 0x20400060, 
 
 /* Overlay Smoothing Type $B%*!<%P%l%$J?3j%?%$%W(B */
 Overlay_Smoothing_Type = 0x20400070, 
 
 /* Overlay Foreground Density $B%*!<%P%l%$A0LLG;EY(B */
 Overlay_Foreground_Density = 0x20400080,
 
 /* Overlay Mode $B%*!<%P%l%$%b!<%I(B */
 overlay_Mode = 0x20400090,
 
 /* Threshold Density $BogCMG;EY(B */
 Threshold_Density = 0x20400100,
 
 /* Referenced Image Box Sequence $B;2>H2hA|%\%C%/%9%7!<%1%s%9(B */
 Referenced_Image_Box_Sequence = 0x20400500,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_2100_Length = 0x21000000,
 
 /* Execution Status $B<B9T>uBV(B */
 Execution_Status = 0x21000020, 
 
 /* Execution Status Info $B<B9T>uBV>pJs(B */
 Execution_Status_Info = 0x21000030,
 
 /* Creation Date $B:n@.F|IU(B */
 Creation_Date = 0x21000040,
 
 /* Creation Time $B:n@.;~9o(B */
 Creation_Time = 0x21000050,
 
 /* Originator $BH/9T85(B */
 Originator = 0x21000070,
 
 /* Referenced Print Job Sequence $B;2>H%W%j%s%H%8%g%V%7!<%1%s%9(B */
 Referenced_Print_Job_Sequence = 0x21000500,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_2110_Length = 0x21100000,
 
 /* Printer Status $B%W%j%s%?>uBV(B */
 Printer_Status = 0x21100010, 
 
 /* Printer Status Info $B%W%j%s%?>uBV>pJs(B */
 Printer_Status_Info = 0x21100020,
 
 /* Printer Name $B%W%j%s%?L>(B */
 Printer_Name = 0x21100030,
 
 /*************************************************************************/
 /* Group Length  */
 Group_4000_Length = 0x40000000,  
 
 /* Arbitrary  */
 Arbitray = 0x40000010,
 
 /* Comments  */
 Group_4000_Comments = 0x40004000, 
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_4008_Length = 0x40080000,
 
 /* Results ID $B7k2L(BID */
 Results_ID = 0x40080040,
 
 /* Results ID Issuer $B7k2L(BID$BH/9T<T(B */
 Results_ID_Issuer = 0x40080042,
 
 /* Referenced Interpretation Sequence $B;2>H2r<a%7!<%1%s%9(B */
 Referenced_Interpretation_Sequence = 0x40080050, 
 
 /* Interpretation Recorded Date $B2r<a5-O?F|IU(B */
 Interpretation_Recorded_Date = 0x40080100,
 
 /* Interpretation Recorded Time $B2r<a5-O?;~9o(B */
 Interpretation_Recorded_Time = 0x40080101, 
 
 /* Interpretation Recorder $B2r<a5-O?<T(B */
 Interpretation_Recorder = 0x40080102, 
 
 /* Reference to Recorded Sound $BO?2;$5$l$?2;@<$X$N;2>H(B */
 Reference_to_Recorded_Sound = 0x40080103,  
 
 /* Interpretation Transcription Date $B2r<aE><LF|IU(B */
 Interpretation_Transcription_Date = 0x40080108, 
 
 /* Interpretation Transcription Time $B2r<aE><L;~9o(B */
 Interpretation_Transcription_Time = 0x40080109,
 
 /* Interpretation Transcriber $B2r<aE><L<T(B */
 Interpretation_Transcriber = 0x4008010A,
 
 /* Interpretation Text $B2r<a%F%-%9%H(B */
 Interpretation_Text = 0x4008010B, 
 
 /* Interpretation Author $B2r<a:n@.<T(B */
 Interpretation_Author = 0x4008010C,
 
 /* Interpretation Approver Sequence $B2r<a>5G'<T%7!<%1%s%9(B */
 Interpretation_Approver_Sequence = 0x40080111,  
 
 /* Interpretation Approval Date $B2r<a>5G'F|IU(B */
 Interpretation_Approval_Date = 0x40080112,
 
 /* Interpretation Approval Time $B2r<a>5G';~9o(B */
 Interpretation_Approval_Time = 0x40080113,
 
 /* Physician Approving Interpretation $B2r<a>5G'0e;U(B */
 Physician_Approving_Interpretation = 0x40080114, 
 
 /* Interpretation Diagnosis Description $B2r<a?GCG5-=R(B */
 Interpretation_Diagnosis_Description = 0x40080115, 
 
 /* Diagnosis Code Sequence $B2r<a?GCG%3!<%I%7!<%1%s%9(B */
 Diagnosis_Code_Sequence = 0x40080117, 
 
 /* Results Distribution List Sequence $B7k2LG[I[%j%9%H%7!<%1%s%9(B */
 Results_Distribution_List_Sequence = 0x40080118,  
 
 /* Distribution Name $BG[I[L>A0(B */
 Distribution_Name = 0x40080119, 
 
 /* Distribution Address $BG[I[=;=j(B */
 Distribution_Address = 0x4008011A,  
 
 /* Interpretation ID $B2r<a(BID */
 Interpretation_ID = 0x40080200,
 
 /* Interpretation ID Issuer $B2r<a(BID$BH/9T<T(B */
 Interpretation_ID_Issuer = 0x40080202,
 
 /* Interpretation Type ID $B2r<a%?%$%W(BID */
 Interpretation_Type_ID = 0x40080210,
 
 /* Interpretation Status ID $B2r<a>uBV(BID */
 Interpretation_Status_ID = 0x40080212, 
 
 /* Impressions $B0u>](B */
 Impression = 0x40080300, 
 
 /* Results Comments $B7k2L%3%a%s%H(B */
 Group_4008_Comments = 0x40084000,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_5000_Length = 0x50000000,

 /* Curve Dimensions $B%+!<%VDj5A(B */
 Curve_Dimensions = 0x50000005, 
 
 /* Number of Points $BE@$N?t(B */
 Number_of_Points = 0x50000010,
 
 /* Type of Data $B%G!<%?$N%?%$%W(B */
 Type_of_Data = 0x50000020,
 
 /* Curve Description $B%+!<%V5-=R(B */
 Curve_Description = 0x50000022, 
 
 /* Axis Units $B<4C10L(B */
 Axis_Units = 0x50000030,
 
 /* Axis Labels $B<4%i%Y%k(B */
 Axis_Labels = 0x50000040,  
 
 /* Data Value Representation $B%G!<%?CMI=8=(B */
 Data_Value_Representation = 0x50000103, 
 
 /* Minimum Coordinate Value $B:G>.:BI8CM(B */
 Minimum_Coordinate_Value = 0x50000104,
 
 /* Maximum Coordinate Value $B:GBg:BI8CM(B */
 Maximum_Coordinate_Value = 0x50000105,
 
 /* Curve Range $B%+!<%VHO0O(B */
 Curve_Range = 0x50000106,
 
 /* Curve Data Descriptor $B%+!<%V%G!<%?5-=R(B */
 Curve_Data_Descriptor = 0x50000110,
 
 /* Coordinate Start Value $B:BI83+;OCM(B */
 Coordinate_Start_Value = 0x50000112,
 
 /* Coordinate Step Value $B:BI84V3VCM(B */
 Coordinate_Step_Value = 0x50000114,
 
 /* Audio Type $B%*!<%G%#%*%?%$%W(B */
 Audio_Type = 0x50002000, 
 
 /* Audio Sample Format $B%*!<%G%#%*%5%s%W%k%U%)!<%^%C%H(B */
 Audio_Sample_Format = 0x50002002, 
 
 /* Number of Channels $B%A%c%M%k$N?t(B */
 Number_of_Channels = 0x50002004,
 
 /* Number of Samples $B%5%s%W%k$N?t(B */
 Number_of_Samples = 0x50002006, 
 
 /* Sample Rate $B%5%s%W%kN((B */
 Sample_Rate = 0x50002008,
 
 /* Total Time $BAm;~4V(B */
 Total_Time = 0x5000200A,
 
 /* Audio Sample Data $B%*!<%G%#%*%5%s%W%k%G!<%?(B */
 Audio_Sample_Data = 0x5000200C, 
 
 /* Audio Comments $B%*!<%G%#%*%3%a%s%H(B */
 Audio_Comments = 0x5000200E,
 
 /* Curve Data $B%+!<%V%G!<%?(B */
 Curve_Data = 0x50003000,
 
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9(B */
 Group_6000_Length = 0x60000000,
 
 /* Overlay_Rows $B%*!<%P!<%l%$9T(B */
 Overlay_Rows = 0x60000010, 
 
 /* Overlay_Columns $B%*!<%P!<%l%$Ns(B */
 Overlay_Columns = 0x60000011,
 
 /* Number of Frames in Overlay $B%*!<%P!<%l%$$NCf$N%U%l!<%`$N?t(B */
 Number_of_Frames_in_Overlay = 0x60000015, 
 
 /* Overlay Type $B%*!<%P!<%l%$%?%$%W(B */
 Overlay_Type = 0x60000040, 
 
 /* Origin $B86E@(B */
 Origin = 0x60000050,
 
 /* Compression Code */
 Compression_Code = 0x60000060, 
 
 /* Bits Allocated $B3d$jEv$F%S%C%H(B */
 Bits_Allocated_6000 = 0x60000100,
 
 /* Bit Position $B%S%C%H0LCV(B */
 Bit_Position = 0x60000102,
 
 /* Overlay Format */
 Overlay_Format = 0x60000110,
 
 /* Overlay Location */
 Overlay_Location = 0x60000200, 
 
 /* Overlay Descriptor Gray $B%*!<%P!<%l%$5-=R;R(B-$B%0%l%$(B */
 Overlay_Descriptor_Gray = 0x60001100, 
 
 /* Overlay Descriptor Red $B%*!<%P!<%l%$5-=R;R(B-$B@V(B */
 Overlay_Descriptor_Red = 0x60001101,
 
 /* Overlay Descriptor Green $B%*!<%P!<%l%$5-=R;R(B-$BNP(B */
 Overlay_Descriptor_Green = 0x60001102, 
 
 /* Overlay Descriptor Blue $B%*!<%P!<%l%$5-=R;R(B-$B@D(B */
 Overlay_Descriptor_Blue = 0x60001103,
 
 /* Overlays Gray $B%*!<%P!<%l%$(B-$B%0%l%$(B */
 Overlays_Gray = 0x60001200, 
 
 /* Overlays Red $B%*!<%P!<%l%$(B-$B@V(B */
 Overlays_Red = 0x60001201,
 
 /* Overlays Green $B%*!<%P!<%l%$(B-$BNP(B */
 Overlays_Green = 0x60001202,
 
 /* Overlays Blue $B%*!<%P!<%l%$(B-$B@D(B */
 Overlays_Blue = 0x60001203,
 
 /* ROI Area ROI$BLL@Q(B */
 ROI_Area = 0x60001301,
 
 /* ROI Mean ROI$BJ?6Q(B */
 ROI_Mean = 0x60001302,
 
 /* ROI Standard Deviation ROI$BI8=`JP:9(B */
 ROI_Standard_Deviation = 0x60001303,
 
 /* Overlay Data $B%*!<%P!<%l%$%G!<%?(B */
 Overlay_Data = 0x60003000,
 
 /* Comments */
 Group_6000_Comments = 0x60004000,
  
 /*************************************************************************/
 /* Group Length $B%0%k!<%WD9$5(B */
 Group_7FE0_Length = 0x7FE00000,
 
 /* Pixel Data $B2hAG%G!<%?(B */
 Pixel_Data = 0x7FE00010, 
 
 /* Item $B9`L\(B */
 Item = 0xFFFEE000,
 
 /* Item Delimitation Item $B9`L\6h@Z$j9`L\(B */
 Item_Delimitation_Item = 0xFFFEE00D,
 
 /* Sequence Delimitation Item $B%7!<%1%s%96h@Z$j9`L\(B */
 Sequence_Delimitation_Item = 0xFFFEE0DD
} DICOM_TAG;

/* 
* Patient_Orientation		= 0x00200020 ,  $B45<TJ}8~(B
*  $B2hA|LL$K4XO"$7$?45<TJ}8~$O!"@5$N9T<4$*$h$S@5$NNs<4$N2rK63XE*$JJ}8~$r(B
*  $B<($9(B2 $B$D$NCM$K$h$C$F;XDj$5$l$k!#:G=i$NEPO?$O9T$NJ}8~$G$"$j!":G=i$N9T(B
*  $B$NCf$N:G=i$N2hAG$+$i$=$N9T$NCf$N:G8e$N2hAG$NJ}8~$K$h$C$FM?$($i$l$k!#(B
*  2 $BHVL\$NEPO?$ONs$NJ}8~$G$"$j!":G=i$NNs$N:G=i$N2hAG$+$i$=$NNs$NCf$N:G(B
*  $B8e$N2hAG$NJ}8~$K$h$C$FM?$($i$l$k!#(B
*  $B2rK63XE*J}8~$OF,J8;z$K$h$C$F;XDj$5$l$k!#(B
*   A(anterior:$BA0LL(B) P(posterior:$BGXLL(B)
*   R(right:$B1&(B) L(left:$B:8(B)
*	H(head:$BF,(B)	F(foot:$BB-(B)
*  $BJ}8~$NB0@-$N3F!9$NCM$O!">/$J$/$H$b$3$l$i$NJ8;z$N(B1$B$D$r4^$`!#(B
*/

/********************************************************************* */
 
RC read_dicom(FILE *fp, DICOM_HEADER *dicom, GRAPHIC *gr);
RC read_dicom_header(FILE *fp, DICOM_HEADER *dicom, int *pimage);
RC read_dicom_image(FILE *fp, DICOM_HEADER dicom, GRAPHIC *gr, int pimage);
void init_dicom(DICOM_HEADER *dicom);
void print_dicom_header(FILE *fp, DICOM_HEADER dicom);
RC print_dicom_tag(FILE *fp, DICOM_HEADER dicom);
RC set_dicom_header(DICOM_HEADER dicom, GRAPHIC *gr);
RC conv_2complement_to_mono16(GRAPHIC src, GRAPHIC *dst, DICOM_HEADER dicom);

#endif /*  DICOM_UTL_H  */
 
