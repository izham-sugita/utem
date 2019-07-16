/*********************************************************************
 * dicom_element.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA> <Yoshihito KURIZUKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: dicom_element.h 864 2007-01-17 10:02:04Z sasaoka $ */

#ifndef DICOM_ELEMENT_H
#define DICOM_ELEMENT_H

typedef enum{
/* OB,OW,SQまたはUN */
	OB, /* その他のバイト列 */
	OW, /* その他のワード列 */
	SQ, /* 項目のシーケンス */
	UN, /* 未知 */
/* OB,OW,SQまたはUN以外の明示的VR */
/* (VR名), (VR名の意味)  (値の長さ) */
	AE, /* 応用エンディティ  (max 16 bytes) */
	AS, /* 年令列            (4bytes char) */
	AT, /* 属性タグ          (4bytes) */
	CS, /* コード列          (max 16 bytes) */
	DA, /* 日付              (8bytes) */
	DS, /* 10進数列          (max 16 bytes) */
	DT, /* 日時              (max 16 bytes) */
	FL, /* 単精度浮動小数点  (4bytes) */
	IS, /* 倍精度浮動小数点  (8bytes) */
	LO, /* 長列              (max 16 bytes) */
	LT, /* 長テキスト        (max 10240 char) */
	PN, /* 人名              (max 64 char) */
	SH, /* 短列              (max 16 char) */
	SL, /* 符号付長整数      (4bytes) */
	SS, /* 符号付短整数      (2bytes) */
	ST, /* 短テキスト        (max 1024 bytes) */
	TM, /* 時間              (max 16 bytes) */
	UI, /* 個有識別子        (max 64 bytes) */
	UL, /* 符号無し長整数    (4bytes) */
	US, /* 符号無し短整数    (2bytes) */
	UT, /* 無制限テキスト    (2^32-2 ?) */
/* NONE は VR無し(暗示的) */
	NONE /* VR無し */
} VR;

typedef struct {
	unsigned long tag;
	VR vr;
	char name[DICOM_STR_SIZE];
} DICOM_ELEM;

DICOM_ELEM Dicom_elem[] = {
{0x00000000, UL, "Group Length 0000"}, /* グループ長さ 0000 */

{0x00020000, UL, "Group_0002_Length"},             /* グループ長さ 0002 */
{0x00020001, OB, "File_Meta_Information_Version"}, /* ファイルメタ情報版 */
{0x00020002, UI, "Media Stored SOP Class UID "},   /* 媒体保存SOPクラスUID */
{0x00020003, UI, "Media Stored SOP Instance UID "},
                                          /* 媒体保存SOPインスタンスUID */
{0x00020010, UI, "Transfer_Syntax_UID "},        /* 転送構文UID */
{0x00020012, UI, "Implementation_Class_UID"},    /* 実装クラスUID */
{0x00020013, SH, " Implementation_Version_Name"},/* 実装版名 */

{0x00080000, UL, "Group_0008_Length"},           /* グループ長さ 0008 */ 
{0x00080001, NONE, "Length to End"},             /* */
{0x00080005, CS, "Specific_Character_Set"},      /* 特定文字集合 */
{0x00080008, CS, "Image_Type" },                 /* 画像タイプ */
{0x00080010, NONE, "Recognition Code"},          /* */
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
{0x00080040, NONE, "Data Set Type"},             /* */
{0x00080041, NONE, "Data Set Subtype"},          /* */
{0x00080042, CS, "Nuclear_Medicine_Series_Type"},/* 核医学シリーズタイプ */
{0x00080050, SH, "Accession_Number"},            /* 受付番号 */
{0x00080052, CS, "Query_Retrieve_Level"},        /* 間合せ/取得レベル */
{0x00080054, AE, "Retrieve_AE_Title"},           /* 取得AE名称 */
{0x00080058, UI, "Failed_SOP_Instance_UID_List"},
                                          /* 失敗SOPインスタンスUIDリスト */
{0x00080060, CS, "Modality"},                        /* モダリティ */
{0x00080061, CS, "Modality in Study"},               /* 検査の中のモダリティ */
{0x00080064, CS, "Conversion_Type"},                 /* 変換形式 */
{0x00080068, CS, "Presentation Intent Type"},        /* 提示意図タイプ */
{0x00080070, LO, "Manufacturer"},                    /* 製造者 */
{0x00080080, LO, "Institution_Name"},                /*  施設名 */
{0x00080081, ST, "Institution_Address"},             /* 施設住所 */
{0x00080082, SQ, "Institution_Code_Sequence"},       /* 施設コードシーケンス */
{0x00080090, PN, "Referring_Physicians_Name"},       /* 照会医師の名前 */
{0x00080092, ST, "Referring_Physicians_Address"},    /* 紹介医師の住所 */
{0x00081010, SH, "Station_Name"},                    /* ステーション名 */
{0x00081030, LO, "Study_Description"},               /* 検査記述 */
{0x00081040, LO, "Institutional_Department_Name"},   /* 施設部門名 */
{0x00081050, PN, "Attending_Physicians_Name"},       /* 実施医師の名前 */
{0x00081060, PN, "Name_of_Physician_Reading_Study"}, /* 検査読影医師の名前 */
{0x00081070, PN, "Operators_Name"},                  /* 操作者の名前 */
{0x00081080, LO, "Admitting_Diagnoses_Description"}, /* 受診時診断記述 */
{0x00081090, LO, "Manufacturers_Model_Name"},        /* 製造者のモデル名 */
{0x00082120, SH, "Stage_Name"},                      /* ステージ名 */
{0x00082122, IS, "Stage_Number"},                    /* ステージ番号 */
{0x00082124, IS, "Number_of_Stages"},                /* ステージの数 */
{0x00082128, IS, "View_Number"},                     /* ビュー番号 */
{0x0008212a, IS, "Number_of_Views_in_Stage"},    /* ステージの中のビューの数 */

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

{0x00180000, UL, "Group_0018_Length"},        /* グループ長さ 0018 */
{0x00180010, LO, "Contrast_Bolus_Agent"},     /* 造影剤/ボーラス剤 */
{0x00180015, CS, "Body_Part_Examined"},       /* 検査される部位 */
{0x00180020, CS, "Scanning_Sequence"},        /* スキャニングシーケンス */
{0x00180021, CS, "Sequence_Variant"},         /* シーケンス変形 */
{0x00180022, CS, "Scan_Options"},             /* スキャンオプション */
{0x00180023, CS, "MR_Acquisition_Type"},      /* MR収集タイプ */
{0x00180024, SH, "Sequence_Name"},            /* シーケンス名 */
{0x00180050, DS, "Slice_Thickness"},          /* スライス厚さ */
{0x00180060, DS, "KVP"},                      /* */
{0x00180081, DS, "Echo_Time"},                /* エコー時間 */
{0x00180083, DS, "Number_of_Averages"},       /* 平均化の回数 */
{0x00180086, IS, "Echo_Numbers"},             /* エコー番号 */
{0x00180091, IS, "Echo_Train_Length"},        /* エコートレン長さ */
{0x00181000, LO, "Device_Serial_Number"},     /* 装置シリアル番号 */
{0x00181020, LO, "Software_Versions "},       /* ソフトウェア版 */
{0x00181030, LO, "Protocol_Name "},           /* プロトコル名 */
{0x00181100, DS, "Reconstruction_Diameter "}, /* 再構成直径 */
{0x00181120, DS, "Gantry_Detector_Tilt "},    /* 架台 検出器傾 */
{0x00181130, DS, "Table_Height"},             /* テーブル高さ */
{0x00181140, CS, "Rotation_Direction"},       /* 回転方向 */
{0x00181150, DS, "Exposure_Time"},            /* 曝射時間 */
{0x00181151, IS, "X_ray_Tube_Current "},      /* X線管電流 */
{0x00181152, IS, "Exposure "},                /* 曝射量 */
{0x00181210, DS, "Convolution_Kernel"},       /* コンボリーションカーネル */
{0x00181250, SH, "Receiving_Coil"},           /* 受信コイル */
{0x00181251, SH, "Transmitting_Coil"},        /* 送信コイル */
{0x00181314, DS, "Flip_Angle"},               /* フリップ角 */
{0x00185100, CS, "Patient_Position"},         /* 患者位置 */
{0x00185101, CS, "View_Position"},            /* 視野位置 */

{0x00200000, UL, "Group_0020_Length"},            /* グループ長さ 2000 */
{0x0020000d, UI, "Study_Instance_UID"},           /* 検査インスタンスID */
{0x0020000e, UI, "Series_Instance_UID"},          /* シリーズインスタンスID */
{0x00200010, SH, "Study_ID"},                     /* 検査ID */
{0x00200011, IS, "Series_Number"},                /* シリーズ番号 */
{0x00200012, IS, "Acquisition_Number"},           /* 収集番号 */
{0x00200013, IS, "Image_Number"},                 /* インスタンス番号 */
{0x00200020, CS, "Patient_Orientation"},          /* 患者方向 */
{0x00200032, DS, "Image_Position_Patient"},       /* 画像位置(患者) */
{0x00200037, DS, "Image_Orientation_Patient"},    /* 画像方向(患者) */
{0x00200052, UI, "Frame_of_Reference_UID"},       /* 基準座標系UID */
{0x00200060, CS, "Laterality"},                   /* 左右 */
{0x00201002, IS, "Images_in_Acquisition"},        /* 収集の中の画像 */
{0x00201040, LO, "Position_Reference_Indicator"}, /* 位置基準標識 */
{0x00201041, DS, "Slice_Location"},               /* スライス位置 */
{0x00204000, LT, "Image_Comments"},               /* 画像コメント */

{0x00280000, UL, "Group_0028_Length"},          /* グループ長さ 0028 */
{0x00280002, US, "Samples_per_Pixel"},          /* 画素あたりサンプル */
{0x00280004, CS, "Photometric_Interpretation"}, /* 光度測定解釈 */
{0x00280006, US, "Planar_Configuration"},       /* 面構成 */
{0x00280010, US, "Rows"},                       /* 横行 */
{0x00280011, US, "Columns"},                    /* 縦行 */
{0x00280030, DS, "Pixel_Spacing"},              /* 画素間隔 */
{0x00280034, IS, "Pixel_Aspect_Ratio"},         /* アスペクト比 */
{0x00280100, US, "Bits_Allocated"},             /* 割り当てビット */
{0x00280101, US, "Bits_Stored"},                /* 格納ビット */
{0x00280102, US, "High_Bit"},                   /* 高位ビット */
{0x00280103, US, "Pixel_Representation"},       /* 画素表現 */
{0x00280106, US, "Smallest_Image_Pixel_Value"}, /* 最小画像画素値 */
{0x00280107, US, "Largest_Image_Pixel_Value"},  /* 最大画像画素値 */
{0x00281050, DS, "Window_Center"},              /* ウィンドウ中心 */
{0x00281051, DS, "Window_Width"},               /* ウィンドウ幅 */
{0x00281052, DS, "Rescale_Intercept"},          /* リスケール切片 */
{0x00281053, DS, "Rescale_Slope"},              /* リスケール傾斜 */

{0x7FE00000, UL, "Group_7FE0_Length"},         /* グループ長さ 7FE0 */
{0x7FE00010, OW, "Pixel_Data"}                 /* 画素データ */
};

#endif /* DICOM_ELEMENT_H */

