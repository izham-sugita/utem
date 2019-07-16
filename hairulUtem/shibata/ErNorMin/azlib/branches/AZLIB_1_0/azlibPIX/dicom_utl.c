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
	{0x00000000, UL, "Group Length 0000"},/* グループ長さ 0000 */

	{0x00020000, UL, "Group_0002_Length"},            /* グループ長さ 0002 */
	{0x00020001, OB, "File_Meta_Information_Version"},/* ファイルメタ情報版 */
	{0x00020002, UI, "Media Stored SOP Class UID "},  /* 媒体保存SOPクラスUID */
	{0x00020003, UI, "Media Stored SOP Instance UID "},
	                                          /* 媒体保存SOPインスタンスUID */
	{0x00020010, UI, "Transfer_Syntax_UID "},        /* 転送構文UID */
	{0x00020012, UI, "Implementation_Class_UID"},    /* 実装クラスUID */
	{0x00020013, SH, " Implementation_Version_Name"},/* 実装版名 */
	
	{0x00080000, UL, "Group_0008_Length"},           /* グループ長さ 0008 */ 
	{0x00080001, NONE, "Length to End"},
	{0x00080005, CS, "Specific_Character_Set"},      /* 特定文字集合 */
	{0x00080008, CS, "Image_Type" },                 /* 画像タイプ */
	{0x00080010, NONE, "Recognition Code"},
	{0x00080012, DA, "Instance_Creation_Date"},      /* インスタンス作成日付 */
	{0x00080013, TM, "Instance_Creation_Time"},      /* インスタンス作成時刻 */
	{0x00080014, UI, "Instance_Creator_UID"},        /* インスタンス作成者UID */
	{0x00080016, UI, "SOP_Class_UID"},               /* SOPクラスUID */
	{0x00080018, UI, "SOP_Instance_UID"},            /* SOPインスタンスUID */
	{0x00080020, DA, "Study_Date"},                  /* 検査日付 */
	{0x00080021, DA, "Series_Date"},                 /* シリーズ日付 */
	{0x00080022, DA, "Acquisition_Date"},            /* 収集日付 */
	{0x00080023, DA, "Image_Date"},                  /* 画像日付 */
	{0x00080024, DA, "Overlay_Date"},                /* オーバーレイ日付 */
	{0x00080025, DA, "Curve_Date"},                  /* カーブ日付 */
	{0x00080030, TM, "Study_Time"},                  /* 検査時刻 */
	{0x00080031, TM, "Series_Time"},                 /* シリーズ時刻 */
	{0x00080032, TM, "Acquisition_Time"},            /* 収集時刻 */
	{0x00080033, TM, "Image_Time"},                  /* 画像時刻 */
	{0x00080034, TM, "Overlay_Time"},                /* オーバーレイ時刻 */
	{0x00080035, TM, "Curve_Time"},                  /* カーブ時刻 */
	{0x00080040, NONE, "Data Set Type"},
	{0x00080041, NONE, "Data Set Subtype"},
	{0x00080042, CS, "Nuclear_Medicine_Series_Type"},/* 核医学シリーズタイプ */
	{0x00080050, SH, "Accession_Number"},            /* 受付番号 */
	{0x00080052, CS, "Query_Retrieve_Level"},        /* 間合せ/取得レベル */
	{0x00080054, AE, "Retrieve_AE_Title"},           /* 取得AE名称 */
	{0x00080058, UI, "Failed_SOP_Instance_UID_List"},
	                                          /* 失敗SOPインスタンスUIDリスト */
	{0x00080060, CS, "Modality"},                     /* モダリティ */
	{0x00080061, CS, "Modality in Study"},            /* 検査の中のモダリティ */
	{0x00080064, CS, "Conversion_Type"},              /* 変換形式 */
	{0x00080068, CS, "Presentation Intent Type"},     /* 提示意図タイプ */
	{0x00080070, LO, "Manufacturer"},                 /* 製造者 */
	{0x00080080, LO, "Institution_Name"},             /*  施設名 */
	{0x00080081, ST, "Institution_Address"},          /* 施設住所 */
	{0x00080082, SQ, "Institution_Code_Sequence"},    /* 施設コードシーケンス */
	{0x00080090, PN, "Referring_Physicians_Name"},      /* 照会医師の名前 */
	{0x00080092, ST, "Referring_Physicians_Address"},   /* 紹介医師の住所 */
	{0x00081010, SH, "Station_Name"},                   /* ステーション名 */
	{0x00081030, LO, "Study_Description"},              /* 検査記述 */
	{0x00081040, LO, "Institutional_Department_Name"},  /* 施設部門名 */
	{0x00081050, PN, "Attending_Physicians_Name"},      /* 実施医師の名前 */
	{0x00081060, PN, "Name_of_Physician_Reading_Study"},/* 検査読影医師の名前 */
	{0x00081070, PN, "Operators_Name"},                 /* 操作者の名前 */
	{0x00081080, LO, "Admitting_Diagnoses_Description"},/* 受診時診断記述 */
	{0x00081090, LO, "Manufacturers_Model_Name"},       /* 製造者のモデル名 */
	{0x00082120, SH, "Stage_Name"},                     /* ステージ名 */
	{0x00082122, IS, "Stage_Number"},                   /* ステージ番号 */
	{0x00082124, IS, "Number_of_Stages"},               /* ステージの数 */
	{0x00082128, IS, "View_Number"},                    /* ビュー番号 */
	{0x0008212a, IS, "Number_of_Views_in_Stage"},/* ステージの中のビューの数 */

	{0x00100000, UL, "Group_0010_Length"},     /* グループ長さ 0010 */
	{0x00100010, PN, "Patients_Name"},         /* 患者の名前 */
	{0x00100020, LO, "Patient_ID"},            /* 患者ID */
	{0x00100030, DA, "Patients_Birth_Date"},   /* 患者の誕生日 */
	{0x00100032, DA, "Patients_Birth_Time"},   /* 患者の誕生時刻 */
	{0x00100040, CS, "Patients_Sex"},          /* 患者の性別 */
	{0x00101000, LO, "Other_Patient_IDs"},     /* 患者の他のID */
	{0x00101010, AS, "Patients_Age"},          /* 患者の年令 */
	{0x00101020, DS, "Patients_Size"},         /* 患者の身長 */
	{0x00101030, DS, "Patients_Weight"},       /* 患者の体重 */
	{0x00104000, LT, "Patient_Comments"},      /* 患者コメント */

	{0x00180000, UL, "Group_0018_Length"},       /* グループ長さ 0018 */
	{0x00180010, LO, "Contrast_Bolus_Agent"},    /* 造影剤/ボーラス剤 */
	{0x00180015, CS, "Body_Part_Examined"},      /* 検査される部位 */
	{0x00180020, CS, "Scanning_Sequence"},       /* スキャニングシーケンス */
	{0x00180021, CS, "Sequence_Variant"},        /* シーケンス変形 */
	{0x00180022, CS, "Scan_Options"},            /* スキャンオプション */
	{0x00180023, CS, "MR_Acquisition_Type"},     /* MR収集タイプ */
	{0x00180024, SH, "Sequence_Name"},           /* シーケンス名 */
	{0x00180050, DS, "Slice_Thickness"},         /* スライス厚さ */
	{0x00180060, DS, "KVP"},
	{0x00180081, DS, "Echo_Time"},               /* エコー時間 */
	{0x00180083, DS, "Number_of_Averages"},      /* 平均化の回数 */
	{0x00180086, IS, "Echo_Numbers"},            /* エコー番号 */
	{0x00180091, IS, "Echo_Train_Length"},       /* エコートレン長さ */
	{0x00181000, LO, "Device_Serial_Number"},    /* 装置シリアル番号 */
	{0x00181020, LO, "Software_Versions "},      /* ソフトウェア版 */
	{0x00181030, LO, "Protocol_Name "},          /* プロトコル名 */
	{0x00181100, DS, "Reconstruction_Diameter "},/* 再構成直径 */
	{0x00181120, DS, "Gantry_Detector_Tilt "},   /* 架台 検出器傾 */
	{0x00181130, DS, "Table_Height"},            /* テーブル高さ */
	{0x00181140, CS, "Rotation_Direction"},      /* 回転方向 */
	{0x00181150, DS, "Exposure_Time"},           /* 曝射時間 */
	{0x00181151, IS, "X_ray_Tube_Current "},     /* X線管電流 */
	{0x00181152, IS, "Exposure "},               /* 曝射量 */
	{0x00181210, DS, "Convolution_Kernel"},      /* コンボリーションカーネル */
	{0x00181250, SH, "Receiving_Coil"},          /* 受信コイル */
	{0x00181251, SH, "Transmitting_Coil"},       /* 送信コイル */
	{0x00181314, DS, "Flip_Angle"},              /* フリップ角 */
	{0x00185100, CS, "Patient_Position"},        /* 患者位置 */
	{0x00185101, CS, "View_Position"},           /* 視野位置 */

	{0x00200000, UL, "Group_0020_Length"},          /* グループ長さ 2000 */
	{0x0020000d, UI, "Study_Instance_UID"},         /* 検査インスタンスID */
	{0x0020000e, UI, "Series_Instance_UID"},        /* シリーズインスタンスID */
	{0x00200010, SH, "Study_ID"},                   /* 検査ID */
	{0x00200011, IS, "Series_Number"},              /* シリーズ番号 */
	{0x00200012, IS, "Acquisition_Number"},         /* 収集番号 */
	{0x00200013, IS, "Image_Number"},               /* インスタンス番号 */
	{0x00200020, CS, "Patient_Orientation"},        /* 患者方向 */
	{0x00200032, DS, "Image_Position_Patient"},     /* 画像位置(患者) */
	{0x00200037, DS, "Image_Orientation_Patient"},  /* 画像方向(患者) */
	{0x00200052, UI, "Frame_of_Reference_UID"},     /* 基準座標系UID */
	{0x00200060, CS, "Laterality"},                 /* 左右 */
	{0x00201002, IS, "Images_in_Acquisition"},      /* 収集の中の画像 */
	{0x00201040, LO, "Position_Reference_Indicator"},/* 位置基準標識 */
	{0x00201041, DS, "Slice_Location"},              /* スライス位置 */
	{0x00204000, LT, "Image_Comments"},              /* 画像コメント */

	{0x00280000, UL, "Group_0028_Length"},         /* グループ長さ 0028 */
	{0x00280002, US, "Samples_per_Pixel"},         /* 画素あたりサンプル */
	{0x00280004, CS, "Photometric_Interpretation"},/* 光度測定解釈 */
	{0x00280006, US, "Planar_Configuration"},      /* 面構成 */
	{0x00280010, US, "Rows"},                      /* 横行 */
	{0x00280011, US, "Columns"},                   /* 縦行 */
	{0x00280030, DS, "Pixel_Spacing"},             /* 画素間隔 */
	{0x00280034, IS, "Pixel_Aspect_Ratio"},        /* アスペクト比 */
	{0x00280100, US, "Bits_Allocated"},            /* 割り当てビット */
	{0x00280101, US, "Bits_Stored"},               /* 格納ビット */
	{0x00280102, US, "High_Bit"},                  /* 高位ビット */
	{0x00280103, US, "Pixel_Representation"},      /* 画素表現 */
	{0x00280106, US, "Smallest_Image_Pixel_Value"},/* 最小画像画素値 */
	{0x00280107, US, "Largest_Image_Pixel_Value"}, /* 最大画像画素値 */
	{0x00281050, DS, "Window_Center"},             /* ウィンドウ中心 */
	{0x00281051, DS, "Window_Width"},              /* ウィンドウ幅 */
	{0x00281052, DS, "Rescale_Intercept"},         /* リスケール切片 */
	{0x00281053, DS, "Rescale_Slope"},             /* リスケール傾斜 */
                                                
	{0x7FE00000, UL, "Group_7FE0_Length"},         /* グループ長さ 7FE0 */
	{0x7FE00010, OW, "Pixel_Data"}                 /* 画素データ */
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
				/* 奇数番号を持つ私的データ要素に対する暫定措置 */
				if(((ul>>16) % 2) != 0){
					RC_TRY( skip_tag_elemAE(fp, vr_flag) );
				/* 不明番号タグを表示して終了 */
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
				fprintf(stderr, "Unknown planar_configuration type(面構成)\n");
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


/* OB,OW,SQ,またはUNの明示的VRをもつデータ要素 */
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


/* OB,OW,SQ,またはUN以外の明示的VRをもつデータ要素 */
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


/* 暗黙的VRをもつデータ要素 */
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
	fprintf(fp, "検査日付 = %s\n", dicom.study_date );
	fprintf(fp, "検査時間 = %s\n", dicom.study_time );
	fprintf(fp, "モダリティ = %s\n", dicom.modality);
	fprintf(fp, "製造者名 = %s\n", dicom.manufacturer);
	fprintf(fp, "ステーション名 = %s\n", dicom.station_name);
	fprintf(fp, "患者性別 = %s\n", dicom.patients_sex);
	fprintf(fp, "患者年令 = %s\n", dicom.patients_age);
	fprintf(fp, "スライス厚さ = %8.3f [mm]\n", dicom.slice_thickness);
	fprintf(fp, "架台/検出器傾き = %8.3f [deg]\n", dicom.gantry_detector);
	fprintf(fp, "テーブル高さ = %8.3f [mm]\n", dicom.table_height);
	fprintf(fp, "回転方向 = %s: CW:時計回り  CC:反時計回り＼n",
	             dicom.rotation_direction);
	fprintf(fp, "患者位置 = %s: HFP:Head First-Prone(うつむき)\n"
	            "                 HFS:Head First-Supinei(仰向け)\n"
	            "                 HFDR:Head First-Decubitus Right(臥床姿勢)\n"
	            "                 HFDL:Head First-Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Left\n"
	            "                 FFDL:Feet First_Decubitus Right\n"
	            "                 FFP:Feet First_Prone\n"
	            "                 FFS:Feet First_Supine\n",
	             dicom.pation_position);
	fprintf(fp, "患者方向 = %s: A(anterior:前面), P(posterior:背面)\n"
	            "                 R(right:右), L(left:左)\n"
	            "                 H(head:頭), F(foot:足)\n",
	             dicom.patient_orientation);
	fprintf(fp, "スライス位置 = %8.3f\n", dicom.slice_location);
	fprintf(fp, "画素あたりのサンプル = %ld\n", dicom.samples_per_pixel);
	fprintf(fp, "光度測定解釈 = %s\n", dicom.photometric_interpretation);
	fprintf(fp, "面構成 = %d\n", dicom.planar_configuration);
	fprintf(fp, "横行X size = %ld\n", dicom.rows);
	fprintf(fp, "縦行Y size = %ld\n", dicom.columns);
	fprintf(fp, "画素間隔X = %6.4g [mm]\n", dicom.pixcel_spacing_x);
	fprintf(fp, "画素間隔Y = %6.4g [mm]\n", dicom.pixcel_spacing_y);
	fprintf(fp, "割り当てビット = %ld [bit]\n", dicom.bits_allocated);
	fprintf(fp, "格納ビット = %ld [bit]\n", dicom.bits_stored);
	fprintf(fp, "画素表現 = %ld\n", dicom.pixel_representation);
	fprintf(fp, "最小画素数 = %ld\n", dicom.smallest_image_pixel_value);
	fprintf(fp, "最大画素数 = %ld\n", dicom.largest_image_pixel_value);
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
	fprintf(fp, " Rotation Direction = %s: CW  CC＼n",
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
				fprintf(fp, "%0#10lx [私的要素]\n", dicom.tag[ii1]);
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
 * Pixel_Representation(画素表現)がセットされていると，そのモノクロ
 * 画像は 2の補数表現をしている．それをazlibPIX向けのmono16に変換する
 * Pixel_Representationがセットされていても，2の補数でないピクセルが
 * 含まれていることがあるので注意．各ピクセルごとに調べるしかない．
 */ 
RC conv_2complement_to_mono16(GRAPHIC src, GRAPHIC *dest, DICOM_HEADER dicom)
{
	RC_NULL_CHK(src.data.mono16)
	RC_TRY( copy_graphic_info(src, dest) );
	RC_TRY( allocate_graphic_data(dest) );

	GRAPHIC_SIZE_LOOP_3D(src.size.z, src.size.y, src.size.x)
	{
		/* src 2の補数表現かどうか2重のチェック */
		if( (dicom.pixel_representation != 0) &&
			((src.data.mono16[ii1][ii2][ii3] & (1 << 15)) != 0) ){
			dest->data.mono16[ii1][ii2][ii3] =
			  -( (~(src.data.mono16[ii1][ii2][ii3]-1))&(0xFFFF) );
		}else
			dest->data.mono16[ii1][ii2][ii3] = src.data.mono16[ii1][ii2][ii3];
	}GRAPHIC_SIZE_LOOP_3D_END;
	return(NORMAL_RC);
}

