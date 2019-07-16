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
/* OB,OW,SQまたはUN */
	OB, /* その他のバイト列 */
	OW, /* その他のワード列 */
	SQ, /* 項目のシーケンス */
	UN, /* 未知 */
/* OB,OW,SQまたはUN以外の明示的VR */
/* (VR名), (VR名の意味)  (値の長さ)			*/
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
	UT, /* 無制限テキスト    (2^32-2 ?)*/
/* NONE は VR無し(暗示的) */
	NONE /* VR無し */
} VR;

typedef struct {
	unsigned long tag;
	VR vr;
	char name[DICOM_STR_SIZE];
} DICOM_ELEM;

/* 
*  参考文献: *  DICOM(医療におけるディジタル画像と通信)
*    巻5: データ構造と符号化  表6.2-1 DICOM値表現(p.p.13-17)
*    巻6: データ辞書    6. DICOMデータ要素の登録(p.p.5-21)
*/
/*****************************************************************************
 * DICOM データ要素                                                          *
 *****************************************************************************/
typedef enum {
 /* Group Length グループ長さ */
 Group_0000_Length = 0x00000000,
 
 /* Group Length to End 終わりまでの長さ */
 Group_0000_Length_to_End = 0x00000001,
 
 /* Affected SOP Class UID 影響されるSOPクラスUID */
 Affected_SOP_Class_UID = 0x00000002,
 
 /* Requested SOP Class UID 要求されたSOPクラスUID */
 Requested_SOP_Class_UID = 0x00000003,
 
 /* Recognition Code 認識コード */
 Recognition_Code = 0x00000010,
 
 /* Command Field コマンド領域 */
 Command_Field = 0x00000100,
 
 /* Message ID メッセージID */
 Message_ID = 0x00000110,
 
 /* Message Id being Responded to 応答されているメッセージID */
 Message_Id_being_Responded_to = 0x00000120,
 
 /* Initiator 起動側 */
 Initiator = 0x00000200,
 
 /* Receiver 受信側 */
 Receiver = 0x00000300,
 
 /* Find Location FIND位置 */
 Find_Location = 0x00000400,
 
 /* Move Destination MOVE宛先 */
 Move_Destination = 0x00000600,
 
 /* Priority 優先度 */
 Priority = 0x00000700,
 
 /* Data Set Type データ集合タイプ */
 Data_Set_Type = 0x00000800,
 
 /* Number of Matches 比較の数 */
 Number_of_Matches = 0x00000850, 
 
 /* Response Sequence Number 応答シーケンス番号 */
 Response_Sequence_Number = 0x00000860,
 
 /* Status 状態 */
 Status = 0x00000900,
 
 /* Offending Element 違反要素 */
 Offending_Element = 0x00000901,
 
 /* Error Comment エラーコメント */
 Error_Comment = 0x00000902, 
 
 /* Error ID エラーID */
 Error_ID = 0x00000903,
 
 /* Affected SOP Instance UID 影響されるSOPインスタンスUID */
 Affected_SOP_Instance_UID = 0x00001000,
 
 /* Requested SOP Instance UID 要求されたsopインスタンスUID */
 Requested_SOP_Instance_UID = 0x00001001,
 
 /* Event Type ID イベントタイプID */
 Event_Type_ID = 0x00001002,
 
 /* Attribute Identifier List 属性識別子リスト */
 Attribute_Identifier_List = 0x00001005,
 
 /* Action Type ID 動作タイプID */
 Action_Type_ID = 0x00001008,
 
 /* Requested SOP Instance UID List  */
 Requested_SOP_Instance_UID_List = 0x00001012,
 
 /* Number of Remaining Sub operations 残存副操作の数 */
 Number_of_Remaining_Sub_operations = 0x00001020,
 
 /* Number of Completed Sub operations 完了副操作の数 */
 Number_of_Completed_Sub_operations = 0x00001021,
 
 /* Number of Failed Sub operations 失敗副操作の数 */
 Number_of_Failed_Sub_operations = 0x00001022,
 
 /* Number of Warning Sub operations 警告副操作の数 */
 Number_of_Warning_Sub_operations = 0x00001023,
 
 /* Move Originator Application Entity Title MOVE発行元応用エンティティ名称 */
 Move_Originator_Application_Entity_Title = 0x00001030,
 
 /* Move Originator Message ID MOVE発行元メッセージID */
 Move_Originator_Message_ID = 0x00001031,
 
 /* Message Set ID メッセージ集合ID */
 Message_Set_ID = 0x00005010,
 
 /* End Message Set ID 最後のメッセージ集合ID */
 End_Message_Set_ID = 0x00005020,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0002_Length = 0x00020000,
 
 /* File Meta Information Version ファイルメタ情報版 */
 File_Meta_Information_Version = 0x00020001,
 
 /* Media Stored SOP Class UID 媒体保存sopクラスUID */
 Media_Stored_SOP_Class_UID = 0x00020002,
 
 /* Media Stored SOP Instance UID 媒体保存sopインスタンスUID */
 Media_Stored_SOP_Instance_UID = 0x00020003,
 
 /* Transfer Syntax UID 転送構文UID */
 Transfer_Syntax_UID = 0x00020010,
 
 /* Implementation Class UID 実装クラスUID */
 Implementation_Class_UID = 0x00020012,
 
 /* Implementation Version Name 実装版名 */
 Implementation_Version_Name = 0x00020013,
 
 /* Source Application Entity Title 発生元応用エンティティ名称 */
 Source_Application_Entity_Title = 0x00020016,
 
 /* Private Information Creator UID 私的情報生成者UID */
 Private_Information_Creator_UID = 0x00020100,
 
 /* Private Information 私的情報 */
 Private_Information = 0x00020102,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0004_Length = 0x00040000,
 
 /* File set ID ファイル集合ID */
 File_set_ID = 0x00041130, 
 
 /* File set Desctiptor File File ID ファイル集合記述ファイルID */
 File_set_Descriptor_File_File_ID = 0x00041141,
 
 /* File set Descriptor File Format ファイル集合記述ファイルの特定文字集合 */
 File_set_Descriptor_File_Format = 0x00041142,
 
 /* Root Directory Entitys First Directory Record Offset */
 /* ルートディレクトリエンティティの最初のディレクトリレコードのオフセット */
 Root_Directory_Entitys_First_Directory_Record_Offset = 0x00041200,
 
 /* Root Directory Entitys Last Directory Record Offset */
 /* ルートディレクトリエンティティの最後のディレクトリレコードのオフセット */
 Root_Directory_Entitys_Last_Directory_Record_Offset = 0x00041202,
 
 /* File set Consistence Flag ファイル集合一慣性フラグ */
 File_set_Consistence_Flag = 0x00041212,
 
 /* Directory Record Sequence ディレクトリレコードシーケンス */
 Directory_Record_Sequence = 0x00041220,
 
 /* Next Directory Record Offset 次のディレクトリレコードのオフセット */
 Next_Directory_Record_Offset = 0x00041400,
 
 /* Record In use Flag レコード使用中フラグ */
 Record_In_use_Flag = 0x00041410,
 
 /* Referenced Lower level Directory Entity Offset */
 /* 参照下位ディレクトリエンティティのオフセット */
 Referenced_Lower_level_Directory_Entity_Offset = 0x00041420,
 
 /* Directory Record Type ディレクトリレコードタイプ */
 Directory_Record_Type = 0x00041430,
 
 /* Private Record UID 私的レコードUID */
 Private_Record_UID = 0x00041432,
 
 /* Referenced File ID 参照ファイルID */
 Referenced_File_ID = 0x00041500,
 
 /* Referenced SOP Class UID in File ファイルの中の参照sopクラスUID */
 Referenced_SOP_Class_UID_in_File = 0x00041510,
 
 /* Referenced SOP Instance UID in File */
 /* ファイルの中の参照sopインスタンスUID */
 Referenced_SOP_Instance_UID_in_File = 0x00041511,
 
 /* Number of References 参照の数 */
 Number_of_References = 0x00041600,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0008_Length = 0x00080000,
 
 /* Length to End  */
 Group_0008_Length_to_End = 0x00080001,
 
 /* Specific Character Set 特定文字集合 */
 Specific_Character_Set = 0x00080005,
 
 /* Image Type 画像タイプ */
 Image_Type = 0x00080008,
 
 /* Recognition Code  */
 Recognition_Code_0008 = 0x00080010,
 
 /* Instance Creation Date インスタンス生成日付 */
 Instance_Creation_Date = 0x00080012,
 
 /* Instance Creation Time インスタンス生成時刻 */
 Instance_Creation_Time = 0x00080013,
 
 /* Instance Creator UID インスタンス生成者UID */
 Instance_Creator_UID = 0x00080014,
 
 /* SOP Class UID sopクラスUID */
 SOP_Class_UID = 0x00080016,
 
 /* SOP Instance UID sopインスタンスUID */
 SOP_Instance_UID = 0x00080018,
 
 /* Study Date 検査日付 */
 Study_Date = 0x00080020,
 
 /* Series Date シリーズ日付 */
 Series_Date = 0x00080021,
 
 /* Acquisition Date 収集日付 */
 Acquisition_Date = 0x00080022,
 
 /* Image Date 画像日付 */
 Image_Date = 0x00080023,
 
 /* Overlay Date オーバレイ日付 */
 Overlay_Date = 0x00080024,
 
 /* Curve Date カーブ日付 */
 Curve_Date = 0x00080025,
 
 /* Study Time 検査時刻 */
 Study_Time = 0x00080030,
 
 /* Series Time シリーズ時刻 */
 Series_Time = 0x00080031,
 
 /* Acquisition Time 収集時刻 */
 Acquisition_Time = 0x00080032,
 
 /* Image Time 画像時刻 */
 Image_Time = 0x00080033,
 
 /* Overlay Time オーバレイ時刻 */
 Overlay_Time = 0x00080034,
 
 /* Curve Time カーブ時刻 */
 Curve_Time = 0x00080035,
 
 /* Data Set Type */
 Data_Set_Type_0008 = 0x00080040,
 
 /* Data Set Subtype */
 Data_Set_Subtype = 0x00080041,
 
 /* Nuclear Medicine Series Type 核医学シリーズタイプ */
 Nuclear_Medicine_Series_Type = 0x00080042,
 
 /* Accession Number 受付番号 */
 Accession_Number = 0x00080050,
 
 /* Query/Retrieve Level 問合せ／取得レベル */
 Query_Retrieve_Level = 0x00080052,
 
 /* Retrieve AE Title 取得AE名称 */
 Retrieve_AE_Title = 0x00080054,
 
 /* Failed SOP Instance UID List 失敗sopインスタンスUIDリスト */
 Failed_SOP_Instance_UID_List = 0x00080058,
 
 /* Modality モダリティ */
 Modality = 0x00080060,
 
 /* Conversion Type 変換形式 */
 Conversion_Type = 0x00080064,
 
 /* Manufacturer 製造者 */
 Manufacturer = 0x00080070,
 
 /* Institution Name 施設名 */
 Institution_Name = 0x00080080,
 
 /* Institution Address 施設の住所 */
 Institution_Address = 0x00080081,
 
 /* Institution Code Sequence 施設コードシーケンス */
 Institution_Code_Sequence = 0x00080082,
 
 /* Referring Physician's Name 照会医師名 */
 Referring_Physicians_Name = 0x00080090,
 
 /* Referring Physician's Address 照会医師住所 */
 Referring_Physicians_Address = 0x00080092,
 
 /* Referring Physician's Telephone Numbers 照会医師電話番号 */
 Referring_Physicians_Telephone_Numbers = 0x00080094,
 
 /* Code Value コード値 */
 Code_Value = 0x00080100,
 
 /* Coding Scheme Designator 符号化体系指定子 */
 Coding_Scheme_Designator = 0x00080102,
 
 /* Code Meaning コード意味 */
 Code_Meaning = 0x00080104,
 
 /* Network ID  */
 Network_ID = 0x00081000,
 
 /* Station Name ステーション名 */
 Station_Name = 0x00081010,
 
 /* Study Description 検査記述 */
 Study_Description = 0x00081030,
 
 /* Procedure Code Sequence 処置コードシーケンス */
 Procedure_Code_Sequence = 0x00081032,
 
 /* Series Description シリーズ記述 */
 Series_Description = 0x0008103E,
 
 /* Institutional Department Name 施設部門名 */
 Institutional_Department_Name = 0x00081040,
 
 /* Performing Physician's Name 実施医師の名前 */
 Attending_Physicians_Name = 0x00081050,
 
 /* Name of  Physician(s) Reading Study 検査読影医師名 */
 Name_of_Physician_Reading_Study = 0x00081060,
 
 /* Operators' Name 操作者の名前 */
 Operators_Name = 0x00081070,
 
 /* Admitting Diagnoses Description 受診時診断記述 */
 Admitting_Diagnoses_Description = 0x00081080,
 
 /* Admitting Diagnosis Code Sequence 受診時診断コードシーケンス */
 Admitting_Diagnosis_Code_Sequence = 0x00081084,
 
 /* Manufacturer's Model Name 製造者モデル名 */
 Manufacturers_Model_Name = 0x00081090,
 
 /* Referenced Results Sequence 参照結果シーケンス */
 Referenced_Results_Sequence = 0x00081100,
 
 /* Referenced Study Sequence 参照検査シーケンス */
 Referenced_Study_Sequence = 0x00081110,
 
 /* Referenced Study Component Sequence 参照検査構成要素シーケンス */
 Referenced_Study_Component_Sequence = 0x00081111,
 
 /* Referenced Series Sequence 参照シリーズシーケンス */
 Referenced_Series_Sequence = 0x00081115,
 
 /* Referenced Patient Sequence 参照患者シーケンス */
 Referenced_Patient_Sequence = 0x00081120,
 
 /* Referenced Visit Sequence 参照来院シーケンス */
 Referenced_Visit_Sequence = 0x00081125,
 
 /* Referenced Overlay Sequence 参照オーバレイシーケンス */
 Referenced_Overlay_Sequence = 0x00081130,
 
 /* Referenced Image Sequence 参照画像シーケンス */
 Referenced_Image_Sequence = 0x00081140,
 
 /* Referenced Curve Sequence 参照カーブシーケンス */
 Referenced_Curve_Sequence = 0x00081145,
 
 /* Referenced SOP Class UID 参照sopクラスUID */
 Referenced_SOP_Class_UID = 0x00081150,
 
 /* Referenced SOP Instance UID 参照sopインスタンスUID */
 Referenced_SOP_Instance_UID = 0x00081155,
 
 /* Derivation Description 導出記述 */
 Derivation_Description = 0x00082111,
 
 /* Source Image Sequence 発生元画像シーケンス */
 Source_Image_Sequence = 0x00082112,
 
 /* Stage Name ステージ名 */
 Stage_Name = 0x00082120,
 
 /* Stage Number ステージ番号 */
 Stage_Number = 0x00082122,
 
 /* Number of Stages ステージ数 */
 Number_of_Stages = 0x00082124,
 
 /* Number of Event Timers イベントタイマの数 */
 Number_of_Event_Timers = 0x00082129,
 
 /* View Number ビュー番号 */
 View_Number = 0x00082128,
 
 /* Number of Views in Stage ステージの中のビュー数 */
 Number_of_Views_in_Stage = 0x0008212A,
 
 /* Event Elapsed Time(s) イベント経過時間 */
 Event_Elapsed_Time = 0x00082130,
 
 /* Event Timer Name(s) イベントタイマ名 */
 Event_Timer_Name = 0x00082132,
 
 /* Start Trim 開始トリム */
 Start_Trim = 0x00082142,
 
 /* Stop Trim 停止トリム */
 Stop_Trim = 0x00082143,
 
 /* Recommended Display Frame Rate 推奨表示フレーム率 */
 Recommended_Display_Frame_Rate = 0x00082144,
 
 /* Transducer Position 探触子の位置 */
 Transducer_Position = 0x00082200,
 
 /* Transducer Orientation 探触子の方向 */
 Transducer_Orientation = 0x00082204,
 
 /* Anatomic Structure 解剖学的構造 */
 Anatomic_Structure = 0x00082208, 
 
 /* Comments  */
 Group_0008_Comments = 0x00084000,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0010_Length = 0x00100000,
 
 /* Patient's Name 患者の名前 */
 Patients_Name = 0x00100010,
 
 /* Patient ID 患者ID */
 Patient_ID = 0x00100020,
 
 /* Issuer of Patient ID 患者IDの発行者 */
 Issuer_of_Patient_ID = 0x00100021,
 
 /* Patient's Birth Date 患者の誕生日 */
 Patients_Birth_Date = 0x00100030,
 
 /* Patient's Birth Time 患者の誕生時刻 */
 Patients_Birth_Time = 0x00100032,
 
 /* Patient's Sex 患者の性別 */
 Patients_Sex = 0x00100040,
 
 /*  Patient's Social Security Number 患者の社会保障番号 */
 Patients_Social_Security_Number = 0x00100042,
 
 /* Patient's Insurance Plan Code Sequence */
 /* 患者の保険計画コードシーケンス */
 Patients_Insurance_Plan_Code_Sequence = 0x00100050,
 
 /* Other Patient IDs 患者の他のID */
 Other_Patient_IDs = 0x00101000,
 
 /* Other Patient Names 患者の他の名前 */
 Other_Patient_Names = 0x00101001,
 
 /* Patient's Birth Name 患者の誕生名 */
 Patients_Maiden_Name = 0x00101005,
 
 /* Patient's Age 患者の年齢 */
 Patients_Age = 0x00101010,
 
 /* Patient's Size 患者の身長 */
 Patients_Size = 0x00101020,
 
 /* Patient's Weight 患者の体重 */
 Patients_Weight = 0x00101030,
 
 /* Patient's Address 患者の住所 */
 Patients_Address = 0x00101040,
 
 /* Insurance Plan Identification  */
 Insurance_Plan_Identification = 0x00101050,
 
 /* Patient's Mother's Birth Name 患者の母の誕生名 */
 Patients_Mothers_Maiden_Name = 0x00101060,
 
 /* Military Rank 軍の階級 */
 Military_Rank = 0x00101080,
 
 /* Branch of Service サービス部門 */
 Branch_of_Service = 0x00101081,
 
 /* Medical Record Locator 医療記録所在識別子 */
 Medical_Record_Locator = 0x00101090,
 
 /* Medical Alerts 医療警告 */
 Medical_Alerts = 0x00102000,
 
 /* Contrast Allergies 造影剤アレルギー */
 Contrast_Allergies = 0x00102110,
 
 /* Country of Residence 居住の国 */
 Country_of_Residence = 0x00102150,
 
 /* Region of Residence 居住の地域 */
 Region_of_Residence = 0x00102152, 
 
 /* Patient's Telephone Numbers 患者の電話番号 */
 Patients_Telephone_Numbers = 0x00102154,
 
 /* Ethnic Group 民族グループ */
 Ethnic_Group = 0x00102160,
 
 /* Occupation 職業 */
 Occupation = 0x00102180,
 
 /* Smoking Status 喫煙の状態 */
 Smoking_Status = 0x001021A0,
 
 /* Additional Patient History 患者の追加病歴 */
 Additional_Patient_History = 0x001021B0, 
 
 /* Pregnancy Status 妊娠の状態 */
 Pregnancy_Status = 0x001021C0, 
 
 /* Last Menstrual Date 最終月経日 */
 Last_Menstrual_Date = 0x001021D0,
 
 /* Patient's Religious Preference 患者の信じる宗教 */
 Patients_Religious_Preference = 0x001021F0,
 
 /* Patient Comments 患者コメント */
 Patient_Comments = 0x00104000,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0018_Length = 0x00180000,
 
 /* Contrast/Bolus Agent 造影剤／ボーラス薬剤 */
 Contrast_Bolus_Agent = 0x00180010,
 
 /* Body Part Examined 検査される部位 */
 Body_Part_Examined = 0x00180015,
 
 /* Scanning Sequence スキャニングシーケンス */
 Scanning_Sequence = 0x00180020,
 
 /* Sequence Variant シーケンス変形 */
 Sequence_Variant = 0x00180021,
 
 /* Scan Options スキャンオプション */
 Scan_Options = 0x00180022,
 
 /* MR Acquisition Type ＭＲ収集タイプ */
 MR_Acquisition_Type = 0x00180023, 
 
 /* Sequence Name シーケンス名 */
 Sequence_Name = 0x00180024,
 
 /* Angio Flag アンジオフラグ */
 Angio_Flag = 0x00180025,
 
 /* Radionuclide 放射線核種 */
 Radionuclide = 0x00180030,
 
 /* Radio-pharmaceutical 放射性医薬品 */
 Radiopharmaceutical = 0x00180031,
 
 /* Energy Window Centerline エネルギウインドウ中心線 */
 Energy_Window_Centerline = 0x00180032,
 
 /* Energy Window Total Width エネルギウインドウ全幅 */
 Energy_Window_Total_Width = 0x00180033,
 
 /* Intervention Drug Name インタベンション薬品名 */
 Intervention_Drug_Name = 0x00180034,
 
 /* Intervention Drug Start Time インタベンション薬品開始時刻 */
 Intervention_Drug_Start_Time = 0x00180035,
 
 /* Cine Rate シネのフレーム率 */
 Cine_Rate = 0x00180040, 
 
 /* Slice Thickness スライス厚さ */
 Slice_Thickness = 0x00180050,
 
 /* KVP ＫＶＰ */
 KVP = 0x00180060,
 
 /* Counts Accumulated 積算カウント */
 Counts_Accumulated = 0x00180070, 
 
 /* Acquisition Termination Condition 収集終了条件 */
 Acquisition_Termination_Condition = 0x00180071,
 
 /* Effective Series Duration 実効シリーズ持続時間 */
 Effective_Series_Duration = 0x00180072,
 
 /* Echo Time エコー時間 */
 Echo_Time = 0x00180081, 
 
 /* Inversion Time 反転時間 */
 Inversion_Time = 0x00180082,
 
 /* Number of Averages 平均化回数 */
 Number_of_Averages = 0x00180083,
 
 /* Imaging Frequency 画像周波数 */
 Imaging_Frequency = 0x00180084,
 
 /* Imaged Nucleus 画像原子核 */
 Imaged_Nucleus = 0x00180085, 
 
 /* Echo Number(s) エコー番号 */
 Echo_Numbers = 0x00180086, 
 
 /* Magnetic Field Strength 磁場強度 */
 Magnetic_Field_Strength = 0x00180087,
 
 /* Spacing Between Slices スライス間間隔 */
 Spacing_Between_Slices = 0x00180088,
 
 /* Number of Phase Encoding Steps 位相符号化ステップの数 */
 Number_of_Phase_Encoding_Steps = 0x00180089,
 
 /* Data Collection Diameter データ採集直径 */
 Data_Collection_Diameter = 0x00180090,
 
 /* Echo Train Length エコートレン長さ */
 Echo_Train_Length = 0x00180091,
 
 /* Percent Sampling サンプリング百分率 */
 Percent_Sampling = 0x00180093, 
 
 /* Percent Phase Field of View 位相視野百分率 */
 Percent_Phase_Field_of_View = 0x00180094,
 
 /* Pixel Bandwidth 画素バンド幅 */
 Pixel_Bandwidth = 0x00180095,
 
 /* Device Serial Number 装置シリアル番号 */
 Device_Serial_Number = 0x00181000,
 
 /* Plate ID プレートID */
 Plate_ID = 0x00181004, 
 
 /* Secondary Capture Device ID 二次取得装置ID */
 Secondary_Capture_Device_ID = 0x00181010,
 
 /* Date of Secondary Capture 二次取得の日付 */
 Date_of_Secondary_Capture = 0x00181012,
 
 /* Time of Secondary Capture 二次取得の時刻 */
 Time_of_Secondary_Capture = 0x00181014,
 
 /* Secondary Capture Device Manufacturer 二次取得装置製造者 */
 Secondary_Capture_Device_Manufacturer = 0x00181016,
 
 /* Secondary Capture Device Manufacturer's Model Name */
 /* 二次取得装置製造者型式名 */
 Secondary_Capture_Device_Manufacturers_Model_Name = 0x00181018,
 
 /* Secondary Capture Device Software Version(s)*/
 /* 二次取得装置ソフトウェア版 */
 Secondary_Capture_Device_Software_Version = 0x00181019,
 
 /* Software Version(s) ソフトウェア版 */
 Software_Versions = 0x00181020,
 
 /* Video Image Format Acquired 取得したビデオ画像の形式 */
 Video_Image_Format_Acquired = 0x00181022,
 
 /* Digital Image Format Acquired 取得したデジタル画像の形式 */
 Digital_Image_Format_Acquired = 0x00181023,
 
 /* Protocol Name プロトコル名 */
 Protocol_Name = 0x00181030,
 
 /* Contrast/Bolus Route 造影剤／ボーラス経路 */
 Contrast_Bolus_Route = 0x00181040,
 
 /* Contrast/Bolus Volume 造影剤／ボーラス容量 */
 Contrast_Bolus_Volume = 0x00181041, 
 
 /* Contrast/Bolus Start Time 造影剤／ボーラス開始時刻 */
 Contrast_Bolus_Start_Time = 0x00181042,
 
 /* Contrast/Bolus Stop Time 造影剤／ボーラス停止時刻 */
 Contrast_Bolus_Stop_Time = 0x00181043,
 
 /* Contrast/Bolus Total Dose 造影剤／ボーラス全投与量 */
 Contrast_Bolus_Total_Dose = 0x00181044,
 
 /* Syringe counts 注入カウント */
 Syringe_Counts = 0x00181045,
 
 /* Spatial Resolution 空間分解能 */
 Spatial_Resolution = 0x00181050, 
 
 /* Trigger Time トリガ時間 */
 Trigger_Time = 0x00181060,
 
 /* Trigger Source or Type トリガ源またはタイプ */
 Trigger_Source_or_Type = 0x00181061,
 
 /* Nominal Interval 公称Ｒ−Ｒ間隔 */
 Nominal_Interval = 0x00181062,
 
 /* Frame Time フレーム時間 */
 Frame_Time = 0x00181063,
 
 /* Framing Type フレーミングタイプ */
 Framing_Type = 0x00181064,
 
 /* Frame Time Vector フレーム時間ベクトル */
 Frame_Time_Vector = 0x00181065,
 
 /* Frame Delay フレーム遅れ */
 Frame_Delay = 0x00181066, 
 
 /* Radionuclide Route 放射性核種経路 */
 Radionuclide_Route = 0x00181070, 
 
 /* Radionuclide Volume 放射性核種容量 */
 Radionuclide_Volume = 0x00181071,
 
 /* Radionuclide Start Time 放射性核種開始時間 */
 Radionuclide_Start_Time = 0x00181072,
 
 /* Radionuclide Stop Time 放射性核種停止時間 */
 Radionuclide_Stop_Time = 0x00181073,
 
 /* Radionuclide Total Dose 放射性核種総投与量 */
 Radionuclide_Total_Dose = 0x00181074, 
 
 /* Beat Rejection Flag 拍動除去フラグ */
 Beat_Rejection_Flag = 0x00181080, 
 
 /* Low R-R Value 下限Ｒ−Ｒ値 */
 Low_R_R_Value = 0x00181081,
 
 /* High R-R Value 上限Ｒ−Ｒ値 */
 High_R_R_Value = 0x00181082, 
 
 /* Intervals Acquired 取得間隔 */
 Intervals_Acquired = 0x00181083,
 
 /* Intervals Rejected 除去間隔 */
 Intervals_Rejected = 0x00181084, 
 
 /* PVC Rejection ＰＶＣの除去 */
 PVC_Rejection = 0x00181085, 
 
 /* Skip Beats スキップ拍動 */
 Skip_Beats = 0x00181086,
 
 /* Heart Rate 心拍数 */
 Heart_Rate = 0x00181088,
 
 /* Cardiac Number of Images 心画像の数 */
 Cardiac_Number_of_Images = 0x00181090,
 
 /* Trigger Window トリガウインドウ */
 Trigger_Window = 0x00181094,
 
 /* Reconstruction Diameter 再構成直径 */
 Reconstruction_Diameter = 0x00181100,
 
 /* Distance Source to Detector 線源検出器間距離 */
 Distance_Source_to_Detector = 0x00181110, 
 
 /* Distance Source to Patient 線源患者間距離 */
 Distance_Source_to_Patient = 0x00181111,
 
 /* Gantry/Detector Tilt 架台／検出器傾き */
 Gantry_Detector_Tilt = 0x00181120,
 
 /* Table Height テーブル高さ */
 Table_Height = 0x00181130,
 
 /* Table Traverse テーブルトラバース */
 Table_Traverse = 0x00181131,
 
 /* Rotation Direction 回転方向 */
 Rotation_Direction = 0x00181140,
 
 /* Angular Position 角度位置 */
 Angular_Position = 0x00181141,
 
 /* Radial Position 径方向位置 */
 Radial_Position = 0x00181142,
 
 /* Scan Arc スキャンアーク */
 Scan_Arc = 0x00181143,
 
 /* Angular Step 角度ステップ */
 Angular_Step = 0x00181144,
 
 /* Center of Rotation Offset 回転オフセットの中心 */
 Center_of_Rotation_Offset = 0x00181145, 
 
 /* Rotation Offset 回転オフセット */
 Rotation_Offset = 0x00181146,
 
 /* Field of View Shape 視野の形状 */
 Field_of_View_Shape = 0x00181147,
 
 /* Field of View Dimension(s) 視野の寸法 */
 Field_of_View_Dimensions = 0x00181149, 
 
 /* Exposure Time 曝射時間 */
 Exposure_Time = 0x00181150,
 
 /* X-ray Tube Current Ｘ線管電流 */
 X_ray_Tube_Current = 0x00181151,
 
 /* Exposure 曝射量 */
 Exposure = 0x00181152, 
 
 /* Filter Type フィルタタイプ */
 Filter_Type = 0x00181160,
 
 /* Generator Power 発生装置出力 */
 Generator_Power = 0x00181170,
 
 /* Collimator/grid Name コリメータ／グリッド名 */
 Collimator_grid_Name = 0x00181180, 
 
 /* Collimator Type コリメータタイプ */
 Collimator_Type = 0x00181181,
 
 /* Focal Distance 焦点距離 */
 Focal_Distance = 0x00181182,
 
 /* X Focus Center X 焦点中心 */
 X_Focus_Center = 0x00181183,
 
 /* Y Focus Center Y 焦点中心 */
 Y_Focus_Center = 0x00181184,
 
 /* Focal Spot(s) 焦点 */
 Focal_Spot = 0x00181190, 
 
 /* Date of Last Calibration 最終較正日付 */
 Date_of_Last_Calibration = 0x00181200,
 
 /* Time of Last Calibration 最終較正時刻 */
 Time_of_Last_Calibration = 0x00181201, 
 
 /* Convolution Kernel コンボリューションカーネル */
 Convolution_Kernel = 0x00181210, 
 
 /* Upper/Lower Pixel Values  */
 Upper_Lower_Pixel_Values = 0x00181240,
 
 /* Actual Frame Duration 実フレーム持続時間 */
 Actual_Frame_Duration = 0x00181242, 
 
 /* Count Rate カウント率 */
 Count_Rate = 0x00181243,
 
 /* Receiving Coil 受信コイル */
 Receiving_Coil = 0x00181250,
 
 /* X-ray Tube Current Ｘ線管電流 */
 Transmitting_Coil = 0x00181151,
 
 /* Filter Type フィルタタイプ */
 Screen_Type = 0x00181160, 
 
 /* Phosphor Type 蛍光体タイプ */
 Phosphor_Type = 0x00181261,
 
 /* Scan Velocity スキャン速度 */
 Scan_Velocity = 0x00181300, 
 
 /* Whole Body Technique 全身技術 */
 Whole_Body_Technique = 0x00181301,
 
 /* Scan Length スキャン長さ */
 Scan_Length = 0x00181302,
 
 /* Acquisition Matrix 収集マトリックス */
 Acquisition_Matrix = 0x00181310,
 
 /* Phase Encoding Direction 位相符号化方向 */
 Phase_Encoding_Direction = 0x00181312,
 
 /* Flip Angle フリップ角 */
 Flip_Angle = 0x00181314,
 
 /* Variable Flip Angle Flag 可変フリップ角フラグ */
 Variable_Flip_Angle_Flag = 0x00181315,
 
 /* SAR */
 SAR = 0x00181316, 
 
 /* dB/dt dB/dt */
 dB_dt = 0x00181318,
 
 /* Acquisition Device Processing Description 収集装置処理記述 */
 Acquisition_Device_Processing_Description = 0x00181400, 
 
 /* Acquisition Device Processing Code 収集装置処理コード */
 Acquisition_Device_Processing_Code = 0x00181401, 
 
 /* Cassette Orientation カセッテ方位 */
 Cassette_Orientation = 0x00181402, 
 
 /* Cassette Size カセッテサイズ */
 Cassette_Size = 0x00181403,
 
 /* Exposures on Plate プレート上の曝射量 */
 Exposures_on_Plate = 0x00181404, 
 
 /* Relative X-ray Exposure 相対Ｘ線曝射量 */
 Relative_X_ray_Exposure = 0x00181405,
 
 /* Comments  */
 Group_0018_Comments = 0x00184000,
 
 /* Output Power 出力 */
 Output_Power = 0x00185000,
 
 /* Transducer Data 探触子データ */
 Transducer_Data = 0x00185010,
 
 /* Focus Depth 焦点深さ */
 Focus_Depth = 0x00185012,
 
 /* Preprocessing Function 前処理関数 */
 Preprocessing_Function = 0x00185020,
 
 /* Postprocessing Function 後処理関数 */
 Postprocessing_Function = 0x00185021,
 
 /* Mechanical Index MI */
 Mechanical_Index = 0x00185022,
 
 /* Thermal Index 骨TI */
 Thermal_Index = 0x00185024, 
 
 /* Cranial Thermal Index 頭蓋TI */
 Cranial_Thermal_Index = 0x00185026, 
 
 /* Soft Tissue Thermal Index 軟組織TI */
 Soft_Tissue_Thermal_Index = 0x00185027, 
 
 /* Soft Tissue-focus Thermal Index 軟組織焦点TI */
 Soft_Tissue_focus_Thermal_Index = 0x00185028,
 
 /* Soft Tissue-surface Thermal Index 軟組織表面TI */
 Soft_Tissue_surface_Thermal_Index = 0x00185029,
 
 /* Dynamic Range  */
 Dynamic_Range = 0x00185030, 
 
 /* Total Gain  */
 Total_Gain = 0x00185040, 
 
 /* Depth of Scan Field 表示の深さ */
 Depth_of_Scan_Field = 0x00185050,
 
 /* Patient Position 患者位置 */
 Patient_Position = 0x00185100,
 
 /* View Position 視野位置 */
 View_Position = 0x00185101,
 
 /* Image Transformation Matrix 画像変換マトリックス */
 Image_Transformation_Matrix = 0x00185210,
 
 /* Image Translation Vector 画像変換ベクトル */
 Image_Translation_Vector = 0x00185212, 
 
 /* Sensitivity 感度 */
 Sensitivity = 0x00186000, 
 
 /* Sequence of Ultrasound Regions 超音波領域シーケンス */
 Sequence_of_Ultrasound_Regions = 0x00186011, 
 
 /* Region Spatial Format 領域空間フォーマット */
 Region_Spatial_Format = 0x00186012, 
 
 /* Region Data Type 領域のデータタイプ */
 Region_Data_Type = 0x00186014,
 
 /* Region Flags 領域フラグ */
 Region_Flags = 0x00186016, 
 
 /* Region Location Min X0 領域位置 Min x0 */
 Region_Location_Min_X0 = 0x00186018,
 
 /* Region Location Min Y0 領域位置 Min y0 */
 Region_Location_Min_Y0 = 0x0018601A, 
 
 /* Region Location Max X1 領域位置 Max x1 */
 Region_Location_Max_X1 = 0x0018601C,
 
 /* Region Location Max Y1 領域位置 Max y1 */
 Region_Location_Max_Y1 = 0x0018601E,
 
 /* Reference Pixel X0 参照画素 x0 */
 Reference_Pixel_X0 = 0x00186020, 
 
 /* Reference Pixel Y0 参照画素 y0 */
 Reference_Pixel_Y0 = 0x00186022, 
 
 /* Physical Units X Direction 物理単位 X 方向 */
 Physical_Units_X_Direction = 0x00186024,
 
 /* Physical Units Y Direction 物理単位 Y 方向 */
 Physical_Units_Y_Direction = 0x00186026,
 
 /* Reference Pixel Physical Value X 参照画素物理値 X */
 Reference_Pixel_Physical_Value_X = 0x00186028, 
 
 /* Reference Pixel Physical Value Y 参照画素物理値 Y */
 Reference_Pixel_Physical_Value_Y = 0x0018602A,
 
 /* Physical Delta X 物理変化量 X */
 Physical_Delta_X = 0x0018602C,
 
 /* Physical Delta Y 物理変化量 Y */
 Physical_Delta_Y = 0x0018602E,
 
 /* Transducer Frequency 探触子周波数 */
 Transducer_Frequency = 0x00186030,
 
 /* Transducer Type 探触子タイプ */
 Transducer_Type = 0x00186031,
 
 /* Pulse Repetition Frequency パルス繰返し周波数 */
 Pulse_Repetition_Frequency = 0x00186032, 
 
 /* Doppler Correction Angle ドプラ角度補正 */
 Doppler_Correction_Angle = 0x00186034, 
 
 /* Steering Angle ステアリング角度 */
 Sterring_Angle = 0x00186036,
 
 /* Doppler Sample Volume X Position ドプラサンプル容積 X 位置 */
 Doppler_Sample_Volume_X_Position = 0x00186038, 
 
 /* Doppler Sample Volume Y Position ドプラサンプル容積 Y 位置 */
 Doppler_Sample_Volume_Y_Position = 0x0018603A,
 
 /* TM-Line Position X0 TM線の位置 x0 */
 TM_Line_Position_X0 = 0x0018603C,
 
 /* TM-Line Position Y0 TM線の位置 y0 */
 TM_Line_Position_Y0 = 0x0018603E,
 
 /* TM-Line Position X1 TM線の位置 x1 */
 TM_Line_Position_X1 = 0x00186040, 
 
 /* TM-Line Position Y1 TM線の位置 y1 */
 TM_Line_Position_Y1 = 0x00186042,
 
 /* Pixel Component Organization 画素構成要素の方式 */
 Pixel_Component_Organization = 0x00186044,
 
 /* Pixel Component Mask 画素構成要素マスク */
 Pixel_Component_Mask = 0x00186046,
 
 /* Pixel Component Range Start 画素構成要素範囲の始点 */
 Pixel_Component_Range_Start = 0x00186048,
 
 /* Pixel Component Range Stop 画素構成要素範囲の終点 */
 Pixel_Component_Range_Stop = 0x0018604A, 
 
 /* Pixel Component Physical Units 画素構成要素物理単位 */
 Pixel_Component_Physical_Units = 0x0018604C,
 
 /* Pixel Component Data Type 画素構成要素データタイプ */
 Pixel_Component_Data_Type = 0x0018604E, 
 
 /* Number of Table Break Points 表ブレイクポイントの数 */
 Number_of_Table_Break_Points = 0x00186050,
 
 /* Table of X Break Points X ブレイクポイントの表 */
 Table_of_X_Break_Points = 0x00186052,
 
 /* Table of Y Break Points Y ブレイクポイントの表 */
 Table_of_Y_Break_Points = 0x00186054,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0020_Length = 0x00200000, 
 
 /* Study Instance UID 検査インスタンスUID */
 Study_Instance_UID = 0x0020000D, 
 
 /* Series Instance UID シリーズインスタンスUID */
 Series_Instance_UID = 0x0020000E,
 
 /* Study ID 検査ID */
 Study_ID = 0x00200010, 
 
 /* Series Number シリーズ番号 */
 Series_Number = 0x00200011, 
 
 /* Acquisition Number 収集番号 */
 Acquisition_Number = 0x00200012,
 
 /* Image Number 画像番号 */
 Image_Number = 0x00200013,
 
 /* Isotope Number 同位元素番号 */
 Isotope_Number = 0x00200014, 
 
 /* Phase Number 位相番号 */
 Phase_Number = 0x00200015,  
 
 /* Interval Number 間隔番号 */
 Interval_Number = 0x00200016,
 
 /* Time Slot Number 時間スロット番号 */
 Time_Slot_Number = 0x00200017,  
 
 /* Angle Number 角度番号 */
 Angle_Number = 0x00200018,  
 
 /* Patient Orientation 患者方向 */
 Patient_Orientation = 0x00200020,  
 
 /* Overlay Number オーバレイ番号 */
 Overlay_Number = 0x00200022,
 
 /* Curve Number カーブ番号 */
 Curve_Number = 0x00200024,
 
 /* Image Position  */
 Image_Position_0020 = 0x00200030,
 
 /* Image Position (Patient) 画像位置（患者） */
 Image_Position_Patient = 0x00200032,
 
 /* Image Orientation  */
 Image_Orientation = 0x00200035,
 
 /* Image Orientation (Patient) 画像方向（患者） */
 Image_Orientation_Patient = 0x00200037,
 
 /* Location  */
 Location = 0x00200050,
 
 /* Frame of Reference UID 基準座標系UID */
 Frame_of_Reference_UID = 0x00200052,
 
 /* Laterality 左右 */
 Laterality = 0x00200060,
 
 /* Image Geometry Type  */
 Image_Geometry_Type = 0x00200070,
 
 /* Masking Image  */
 Masking_Image_UID = 0x00200080,
 
 /* Temporal Position Identifier 時間位置識別子 */
 Temporal_Position_Identifier = 0x00200100, 
 
 /* Number of Temporal Positions 時間的位置の数 */
 Number_of_Temporal_Positions = 0x00200105,
 
 /* Temporal Resolution 時間分解能 */
 Temporal_Resolution = 0x00200110,
 
 /* Series in Study 検査の中のシリーズ */
 Series_in_Study = 0x00201000,
 
 /* Acquisitions in Series  */
 Acquisitions_in_Series = 0x00201001,
 
 /* Images in Acquisition 収集の中の画像 */
 Images_in_Acquisition = 0x00201002,
 
 /* Acquisition in Study 検査の中の収集 */
 Acquisition_in_Study = 0x00201004,
 
 /* Reference  */
 Reference = 0x00201020, 
 
 /* Position Reference Indicator 位置参照インジケータ */
 Position_Reference_Indicator = 0x00201040, 
 
 /* Slice Location スライス位置 */
 Slice_Location = 0x00201041, 
 
 /* Other Study Numbers 他の検査番号 */
 Other_Study_Numbers = 0x00201070,
 
 /* Number of Patient Related Studies 患者に関係した検査の数 */
 Number_of_Patient_Related_Studies = 0x00201200,
 
 /* Number of Patient Related Series 患者に関係したシリーズの数 */
 Number_of_Patient_Related_Series = 0x00201202,
 
 /* Number of Patient Related Images 患者に関係した画像の数 */
 Number_of_Patient_Related_Images = 0x00201204,
 
 /* Number of Study Related Series 検査に関係したシリーズの数 */
 Number_of_Study_Related_Series = 0x00201206, 
 
 /* Number of Study Related Images 検査に関係した画像の数 */
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
 
 /* Image Comments 画像コメント */
 Image_Comments = 0x00204000,
 
 /* Original Image Identification  */
 Original_Image_Identification = 0x00205000,
 
 /* Original Image Identification Nomenclature  */
 Original_Image_Identification_Nomenclature = 0x00205002, 
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0028_Length = 0x00280000,
 
 /* Samples per Pixel 画素あたりサンプル */
 Samples_per_Pixel = 0x00280002, 
 
 /* Photometric Interpretation 光度測定解釈 */
 Photometric_Interpretation = 0x00280004, 
 
 /* Image Dimensions  */
 Image_Dimensions = 0x00280005,
 
 /* Planar Configuration 面構成 */
 Planar_Configuration = 0x00280006, 
 
 /* Number of Frames フレームの数 */
 Number_of_Frames = 0x00280008,
 
 /* Frame Increment Pointer フレーム増分ポインタ */
 Frame_Increment_Pointer = 0x00280009,
 
 /* Rows 行 */
 Rows = 0x00280010,
 
 /* Columns 列 */
 Columns = 0x00280011,
 
 /* Pixel Spacing 画素間隔 */
 Pixel_Spacing = 0x00280030,
 
 /* Zoom Factor ズーム率 */
 Zoom_Factor = 0x00280031,
 
 /* Zoom Center ズーム中心 */
 Zoom_Center = 0x00280032,
 
 /* Pixel Aspect Ratio 画素アスペクト比 */
 Pixel_Aspect_Ratio = 0x00280034, 
 
 /* Image Format  */
 Image_Format = 0x00280040,
 
 /* Manipulated Image  */
 Manipulated_Image = 0x00280050,
 
 /* Corrected Image 補正画像 */
 Corrected_Image = 0x00280051,
 
 /* Compression Code  */
 Compression_Code_0028 = 0x00280060,
 
 /* Bits Allocated 割り当てビット */
 Bits_Allocated = 0x00280100,
 
 /* Bits Stored 格納ビット */
 Bits_Stored = 0x00280101,
 
 /* High Bit 高位ビット */
 High_Bit = 0x00280102,
 
 /* Pixel Representation 画素表現 */
 Pixel_Representation = 0x00280103,
 
 /* Smallest Valid Pixel Value  */
 Smallest_Valid_Pixel_Value = 0x00280104, 
 
 /* Largest Valid Pixel Value  */
 Largest_Valid_Pixel_Value = 0x00280105,
 
 /* Smallest Image Pixel Value 最小画像画素値 */
 Smallest_Image_Pixel_Value = 0x00280106,
 
 /* Largest Image Pixel Value 最大画像画素値 */
 Largest_Image_Pixel_Value = 0x00280107, 
 
 /* Smallest Pixel Value in Series シリーズの中の最小画素値 */
 Smallest_Pixel_Value_in_Series = 0x00280108, 
 
 /* Largest Pixel Value in Series シリーズの中の最大画素値 */
 Largest_Pixel_Value_in_Series = 0x00280109, 
 
 /* Pixel Padding Value 画素パディング値 */
 Pixel_Padding_Value = 0x00280120,
 
 /* Image Location  */
 Image_Location = 0x00280200,
 
 /* Window Center ウィンドウ中心 */
 Window_Center = 0x00281050, 
 
 /* Window Width ウィンドウ幅 */
 Window_Width = 0x00281051,
 
 /* Rescale Intercept リスケール切片 */
 Rescale_Intercept = 0x00281052, 
 
 /* Rescale Slope リスケール傾斜 */
 Rescale_Slope = 0x00281053,
 
 /* Rescale type リスケールタイプ */
 Rescale_Type = 0x00281054,
 
 /*  */
 Window_Center_Width_Explanation = 0x00281055, 
 
 /* Gray Scale  */
 Gray_Scale = 0x00281080,
 
 /* Gray Lookup Table Descriptor  */
 Gray_Lookup_Table_Descriptor = 0x00281100,
 
 /* Red Palette Color Lookup Table Descriptor */
 /* レッドパレットカラールックアップテーブル記述子 */
 Red_Palette_Color_Lookup_Table_Descriptor = 0x00281101, 
 
 /* Green Palette Color Lookup Table Descriptor */
 /* グリーンパレットカラールックアップテーブル記述子 */
 Green_Palette_Color_Lookup_Table_Descriptor = 0x00281102,
 
 /* Blue Palette Color Lookup Table Descriptor */
 /*ブルーパレットカラールックアップテーブル記述子 */
 Blue_Palette_Color_Lookup_Table_Descriptor = 0x00281103,
 
 /* Gray Lookup Table Data */
 Gray_Lookup_Table_Data = 0x00281200,
 
 /* Red Palette Color Lookup Table Data */
 /* レッドパレットカラールックアップテーブルデータ */
 Red_Palette_Color_Lookup_Table_Data = 0x00281201,  
 
 /* Green Palette Color Lookup Table Data */
 /* グリーンパレットカラールックアップテーブルデータ */
 Green_Palette_Color_Lookup_Table_Data = 0x00281202,
 
 /* Blue Palette Color Lookup Table Data */
 /* ブルーパレットカラールックアップテーブルデータ */
 Blue_Palette_Color_Lookup_Table_Data = 0x00281203,
 
 /* Modality LUT Sequence モダリティLUTシーケンス */
 Modality_LUT_Sequence = 0x00283000, 
 
 /* LUT Descriptor LUT記述子 */
 LUT_Descriptor = 0x00283002, 
 
 /* LUT Explanation LUT説明 */
 LUT_Explanation = 0x00283003, 
 
 /* Modality LUT Type モダリティLUTタイプ */
 Madality_LUT_Type = 0x00283004,
 
 /* LUT Data LUTデータ */
 LUT_Data = 0x00283006, 
 
 /* VOI LUT Sequence VOI LUTシーケンス */
 VOI_LUT_Sequence = 0x00283010,
 
 /* Comments  */
 Group_0028_Comments = 0x00284000,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0032_Length = 0x00320000,
 
 /* Study Status ID 検査状態ID */
 Study_Status_ID = 0x0032000A,
 
 /* Study Priority ID 検査優先度ID */
 Study_Priority_ID = 0x0032000C, 
 
 /* Study ID Issuer 検査ID発行者 */
 Study_ID_Issuer = 0x00320012,
 
 /* Study Verified Date 検査確認日付 */
 Study_Verified_Date = 0x00320032,  
 
 /* Study Verified Time 検査確認時刻 */
 Study_Verified_Time = 0x00320033, 
 
 /* Study Read Date 検査読影日付 */
 Study_Read_Date = 0x00320034,
 
 /* Study Read Time 検査読影時刻 */
 Study_Read_Time = 0x00320035, 
 
 /* Scheduled Study Start Date 予約検査開始日付 */
 Scheduled_Study_Start_Date = 0x00321000,
 
 /* Scheduled Study Start Time 予約検査開始時刻 */
 Scheduled_Study_Start_Time = 0x00321001,
 
 /* Scheduled Study Stop Date 予約検査終了日付 */
 Scheduled_Study_Stop_Date = 0x00321010,
 
 /* Scheduled Study Stop Time 予約検査終了時刻 */
 Scheduled_Study_Stop_Time = 0x00321011,
 
 /* Scheduled Study Location 予約検査場所 */
 Scheduled_Study_Location = 0x00321020,
 
 /* Scheduled Study Location AE Title(s) 予約検査場所AE名称 */
 Scheduled_Study_Location_AE_Title = 0x00321021,
 
 /* Reason for Study 検査についての理由 */
 Reason__for_Study = 0x00321030,
 
 /* Requesting Physician 依頼側医師 */
 Requesting_Physician = 0x00321032, 
 
 /* Requesting Service 依頼側サービス */
 Requesting_Service = 0x00321033,
 
 /* Study Arrival Date 検査到着日付 */
 Study_Arrival_Date = 0x00321040,
 
 /* Study Arrival Time 検査到着時刻 */
 Study_Arrival_Time = 0x00321041,
 
 /* Study Completion Date 検査完了日付 */
 Study_Completion_Date = 0x00321050,
 
 /* Study Completion Time 検査完了時刻 */
 Study_Completion_Time = 0x00321051,
 
 /* Study Component Status ID 検査構成要素状態ID */
 Study_Component_Status_ID = 0x00321055,
 
 /* Requested Procedure Description 依頼処置記述 */
 Requested_Procedure_Description = 0x00321060,
 
 /* Requested Procedure Code Sequence 依頼処置コードシーケンス */
 Requested_Procedure_Code_Sequence = 0x00321064,
 
 /* Requested Contrast Agent 依頼造影剤 */
 Requested_Contrast_Agent = 0x00321070,
 
 /* Study Comments 検査コメント */
 Study_Comments = 0x00324000, 
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0038_Length = 0x00380000, 
 
 /* Referenced Patient Alias Sequence 参照患者別名シーケンス */
 Referenced_Patient_Alias_Sequence = 0x00380004, 
 
 /* Visit Status ID 来院状態ID */
 Visit_Status_ID = 0x00380008,
 
 /* Admission ID 受診ID */
 Admissin_ID = 0x00380010, 
 
 /* Issuer of Admission ID 受診IDの発行者 */
 Issuer_of_Admission_ID = 0x00380011,
 
 /* Route of Admissions 受診の経路 */
 Route_of_Admissions = 0x00380016,
 
 /* Scheduled Admission Date 予約受診日付 */
 Scheduled_Admissin_Date = 0x0038001A,  
 
 /* Scheduled Admission Time 予約受診時刻 */
 Scheduled_Adission_Time = 0x0038001B,
 
 /* Scheduled Discharge Date 予約退院日付 */
 Scheduled_Discharge_Date = 0x0038001C, 
 
 /* Scheduled Discharge Time 予約退院時間 */
 Scheduled_Discharge_Time = 0x0038001D,
 
 /* Scheduled Patient Institution Residence 予約患者施設内居住場所 */
 Scheduled_Patient_Institution_Residence = 0x0038001E,
 
 /* Admitting Date 受診日付 */
 Admitting_Date = 0x00380020, 
 
 /* Admitting Time 受診時刻 */
 Admitting_Time = 0x00380021, 
 
 /* Discharge Date 退院日付 */
 Discharge_Date = 0x00380030,  
 
 /* Discharge Time 退院時刻 */
 Discharge_Time = 0x00380032,  
 
 /* Discharge Diagnosis Description 退院時診断記述 */
 Discharge_Diagnosis_Description = 0x00380040,  
 
 /* Discharge Diagnosis Code Sequence 退院時診断コードシーケンス */
 Discharge_Diagnosis_Code_Sequence = 0x00380044,
 
 /* Special Needs 特別な介助 */
 Special_Needs = 0x00380050,
 
 /* Current Patient Location 現在の患者の住所 */
 Current_Patient_Location = 0x00380300, 
 
 /* Patient's Institution Residence 患者の施設内居住 */
 Patients_Institution_Residence = 0x00380400,  
 
 /* Patient State 患者の状態 */
 Patient_State = 0x00380500,
 
 /* Visit Comments 来院コメント */
 Visit_Comments = 0x00384000,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_0088_Length = 0x00880000,
 
 /* Storage Media File-set ID 保存媒体ファイル集合ID */
 Storage_Media_File_set_ID = 0x00880130, 
 
 /* Storage Media File-set UID 保存媒体ファイル集合UID */
 Storage_Media_File_set_UID = 0x00880140, 
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_2000_Length = 0x20000000,
 
 /* Number Of Copies コピーの数 */
 Number_of_Copies = 0x20000010,
 
 /* Print Priority プリント優先度 */
 Print_Priority = 0x20000020,
 
 /* Medium Type 媒体タイプ */
 Medium_Type = 0x20000030,
 
 /* Film Destination フィルムあて先 */
 Film_Destination = 0x20000040,
 
 /* Film Session Label フィルムセッションラベル */
 Film_Session_Label = 0x20000050,
 
 /* Memory Allocation メモリ割り当て */
 Memory_Allocation = 0x20000060,
 
 /* Referenced Film Box Sequence 参照フィルムボックスシーケンス */
 Referenced_Film_Box_Sequence = 0x20000500,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_2010_Length = 0x20100000,
 
 /* Image Display Format 画像表示フォーマット */
 Image_Display_Format = 0x20100010,
 
 /* Annotation Display Format ID 注釈表示フォーマットID */
 Annotation_Display_Format_ID = 0x20100030,
 
 /* Film Orientation フィルム方向 */
 Film_Orientation = 0x20100040,
 
 /* Film Size ID フィルムサイズID */
 Film_Size_ID = 0x20100050,
 
 /* Magnification Type 拡大タイプ */
 Magnification_Type = 0x20100060,
 
 /* Smoothing Type 平滑タイプ */
 Smoothing_Type = 0x20100080,
 
 /* Border Density 縁取り濃度 */
 Border_Density = 0x20100100,
 
 /* Empty Image Density 空画像濃度 */
 Empty_Image_Density = 0x20100110, 
 
 /* Min Density 最低濃度 */
 Min_Density = 0x20100120,
 
 /* Max Density 最高濃度 */
 Max_Density = 0x20100130, 
 
 /* Trim ふち飾り */
 Trim = 0x20100140, 
 
 /* Configuration Information 構成情報 */
 Configuration_Information = 0x20100150,
 
 /* Referenced Film Session Sequence 参照フィルムセッションシーケンス */
 Referenced_Film_Session_Sequence = 0x20100500,
 
 /* Referenced Image Box Sequence 参照画像ボックスシーケンス */
 Referenced_Basic_Image_Box_Sequence = 0x20100510, 
 
 /* Referenced Basic Annotation Box Sequence 参照基本注釈ボックスシーケンス */
 Referenced_Basic_Annotation_Box_Sequence = 0x20100520,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_2020_Length = 0x20200000, 
 
 /* Image Position 画像位置 */
 Image_Position = 0x20200010, 
 
 /* Polarity 極性 */
 Polarity = 0x20200020,  
 
 /* Requested Image Size 依頼画像寸法 */
 Requested_Image_Size = 0x20200030,
 
 /* Preformatted Grayscale Image Sequence */
 /* フォーマット済グレイスケール画像シーケンス */
 Preformatted_Greyscale_Image_Sequence = 0x20200110, 
 
 /* Preformatted Color Image Sequence*/
 /* フォーマット済カラー画像シーケンス */
 Preformatted_Color_Image_Sequence = 0x20200111,
 
 /* Referenced Image Overlay Box Sequence */
 /* 参照画像オーバレイボックスシーケンス */
 Referenced_Image_Overlay_Box_Sequence = 0x20200130, 
 
 /* Referenced VOI LUT Box Sequence 参照VOI LUTボックスシーケンス */
 Referenced_VOI_LUT_Sequence = 0x20200140,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_2030_Length = 0x20300000, 
 
 /* Annotation Position 注釈位置 */
 Annotation_Position = 0x20300010,
 
 /* Text String テキスト列 */
 Text_String = 0x20300020,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_2040_Length = 0x20400000,
 
 /* Referenced Overlay Plane Sequence 参照オーバレイ面シーケンス */
 Referenced_Overlay_Plane_Sequence = 0x20400010, 
 
 /* Referenced Overlay Plane Groups 参照オーバレイ面グループ */
 Refenced_Overlay_Plane_Groups = 0x20400011,
 
 /* Overlay Magnification Type オーバレイ拡大タイプ */
 Overlay_Magnification_Type = 0x20400060, 
 
 /* Overlay Smoothing Type オーバレイ平滑タイプ */
 Overlay_Smoothing_Type = 0x20400070, 
 
 /* Overlay Foreground Density オーバレイ前面濃度 */
 Overlay_Foreground_Density = 0x20400080,
 
 /* Overlay Mode オーバレイモード */
 overlay_Mode = 0x20400090,
 
 /* Threshold Density 閾値濃度 */
 Threshold_Density = 0x20400100,
 
 /* Referenced Image Box Sequence 参照画像ボックスシーケンス */
 Referenced_Image_Box_Sequence = 0x20400500,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_2100_Length = 0x21000000,
 
 /* Execution Status 実行状態 */
 Execution_Status = 0x21000020, 
 
 /* Execution Status Info 実行状態情報 */
 Execution_Status_Info = 0x21000030,
 
 /* Creation Date 作成日付 */
 Creation_Date = 0x21000040,
 
 /* Creation Time 作成時刻 */
 Creation_Time = 0x21000050,
 
 /* Originator 発行元 */
 Originator = 0x21000070,
 
 /* Referenced Print Job Sequence 参照プリントジョブシーケンス */
 Referenced_Print_Job_Sequence = 0x21000500,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_2110_Length = 0x21100000,
 
 /* Printer Status プリンタ状態 */
 Printer_Status = 0x21100010, 
 
 /* Printer Status Info プリンタ状態情報 */
 Printer_Status_Info = 0x21100020,
 
 /* Printer Name プリンタ名 */
 Printer_Name = 0x21100030,
 
 /*************************************************************************/
 /* Group Length  */
 Group_4000_Length = 0x40000000,  
 
 /* Arbitrary  */
 Arbitray = 0x40000010,
 
 /* Comments  */
 Group_4000_Comments = 0x40004000, 
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_4008_Length = 0x40080000,
 
 /* Results ID 結果ID */
 Results_ID = 0x40080040,
 
 /* Results ID Issuer 結果ID発行者 */
 Results_ID_Issuer = 0x40080042,
 
 /* Referenced Interpretation Sequence 参照解釈シーケンス */
 Referenced_Interpretation_Sequence = 0x40080050, 
 
 /* Interpretation Recorded Date 解釈記録日付 */
 Interpretation_Recorded_Date = 0x40080100,
 
 /* Interpretation Recorded Time 解釈記録時刻 */
 Interpretation_Recorded_Time = 0x40080101, 
 
 /* Interpretation Recorder 解釈記録者 */
 Interpretation_Recorder = 0x40080102, 
 
 /* Reference to Recorded Sound 録音された音声への参照 */
 Reference_to_Recorded_Sound = 0x40080103,  
 
 /* Interpretation Transcription Date 解釈転写日付 */
 Interpretation_Transcription_Date = 0x40080108, 
 
 /* Interpretation Transcription Time 解釈転写時刻 */
 Interpretation_Transcription_Time = 0x40080109,
 
 /* Interpretation Transcriber 解釈転写者 */
 Interpretation_Transcriber = 0x4008010A,
 
 /* Interpretation Text 解釈テキスト */
 Interpretation_Text = 0x4008010B, 
 
 /* Interpretation Author 解釈作成者 */
 Interpretation_Author = 0x4008010C,
 
 /* Interpretation Approver Sequence 解釈承認者シーケンス */
 Interpretation_Approver_Sequence = 0x40080111,  
 
 /* Interpretation Approval Date 解釈承認日付 */
 Interpretation_Approval_Date = 0x40080112,
 
 /* Interpretation Approval Time 解釈承認時刻 */
 Interpretation_Approval_Time = 0x40080113,
 
 /* Physician Approving Interpretation 解釈承認医師 */
 Physician_Approving_Interpretation = 0x40080114, 
 
 /* Interpretation Diagnosis Description 解釈診断記述 */
 Interpretation_Diagnosis_Description = 0x40080115, 
 
 /* Diagnosis Code Sequence 解釈診断コードシーケンス */
 Diagnosis_Code_Sequence = 0x40080117, 
 
 /* Results Distribution List Sequence 結果配布リストシーケンス */
 Results_Distribution_List_Sequence = 0x40080118,  
 
 /* Distribution Name 配布名前 */
 Distribution_Name = 0x40080119, 
 
 /* Distribution Address 配布住所 */
 Distribution_Address = 0x4008011A,  
 
 /* Interpretation ID 解釈ID */
 Interpretation_ID = 0x40080200,
 
 /* Interpretation ID Issuer 解釈ID発行者 */
 Interpretation_ID_Issuer = 0x40080202,
 
 /* Interpretation Type ID 解釈タイプID */
 Interpretation_Type_ID = 0x40080210,
 
 /* Interpretation Status ID 解釈状態ID */
 Interpretation_Status_ID = 0x40080212, 
 
 /* Impressions 印象 */
 Impression = 0x40080300, 
 
 /* Results Comments 結果コメント */
 Group_4008_Comments = 0x40084000,
 
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_5000_Length = 0x50000000,

 /* Curve Dimensions カーブ定義 */
 Curve_Dimensions = 0x50000005, 
 
 /* Number of Points 点の数 */
 Number_of_Points = 0x50000010,
 
 /* Type of Data データのタイプ */
 Type_of_Data = 0x50000020,
 
 /* Curve Description カーブ記述 */
 Curve_Description = 0x50000022, 
 
 /* Axis Units 軸単位 */
 Axis_Units = 0x50000030,
 
 /* Axis Labels 軸ラベル */
 Axis_Labels = 0x50000040,  
 
 /* Data Value Representation データ値表現 */
 Data_Value_Representation = 0x50000103, 
 
 /* Minimum Coordinate Value 最小座標値 */
 Minimum_Coordinate_Value = 0x50000104,
 
 /* Maximum Coordinate Value 最大座標値 */
 Maximum_Coordinate_Value = 0x50000105,
 
 /* Curve Range カーブ範囲 */
 Curve_Range = 0x50000106,
 
 /* Curve Data Descriptor カーブデータ記述 */
 Curve_Data_Descriptor = 0x50000110,
 
 /* Coordinate Start Value 座標開始値 */
 Coordinate_Start_Value = 0x50000112,
 
 /* Coordinate Step Value 座標間隔値 */
 Coordinate_Step_Value = 0x50000114,
 
 /* Audio Type オーディオタイプ */
 Audio_Type = 0x50002000, 
 
 /* Audio Sample Format オーディオサンプルフォーマット */
 Audio_Sample_Format = 0x50002002, 
 
 /* Number of Channels チャネルの数 */
 Number_of_Channels = 0x50002004,
 
 /* Number of Samples サンプルの数 */
 Number_of_Samples = 0x50002006, 
 
 /* Sample Rate サンプル率 */
 Sample_Rate = 0x50002008,
 
 /* Total Time 総時間 */
 Total_Time = 0x5000200A,
 
 /* Audio Sample Data オーディオサンプルデータ */
 Audio_Sample_Data = 0x5000200C, 
 
 /* Audio Comments オーディオコメント */
 Audio_Comments = 0x5000200E,
 
 /* Curve Data カーブデータ */
 Curve_Data = 0x50003000,
 
 /*************************************************************************/
 /* Group Length グループ長 */
 Group_6000_Length = 0x60000000,
 
 /* Overlay_Rows オーバーレイ行 */
 Overlay_Rows = 0x60000010, 
 
 /* Overlay_Columns オーバーレイ列 */
 Overlay_Columns = 0x60000011,
 
 /* Number of Frames in Overlay オーバーレイの中のフレームの数 */
 Number_of_Frames_in_Overlay = 0x60000015, 
 
 /* Overlay Type オーバーレイタイプ */
 Overlay_Type = 0x60000040, 
 
 /* Origin 原点 */
 Origin = 0x60000050,
 
 /* Compression Code */
 Compression_Code = 0x60000060, 
 
 /* Bits Allocated 割り当てビット */
 Bits_Allocated_6000 = 0x60000100,
 
 /* Bit Position ビット位置 */
 Bit_Position = 0x60000102,
 
 /* Overlay Format */
 Overlay_Format = 0x60000110,
 
 /* Overlay Location */
 Overlay_Location = 0x60000200, 
 
 /* Overlay Descriptor Gray オーバーレイ記述子-グレイ */
 Overlay_Descriptor_Gray = 0x60001100, 
 
 /* Overlay Descriptor Red オーバーレイ記述子-赤 */
 Overlay_Descriptor_Red = 0x60001101,
 
 /* Overlay Descriptor Green オーバーレイ記述子-緑 */
 Overlay_Descriptor_Green = 0x60001102, 
 
 /* Overlay Descriptor Blue オーバーレイ記述子-青 */
 Overlay_Descriptor_Blue = 0x60001103,
 
 /* Overlays Gray オーバーレイ-グレイ */
 Overlays_Gray = 0x60001200, 
 
 /* Overlays Red オーバーレイ-赤 */
 Overlays_Red = 0x60001201,
 
 /* Overlays Green オーバーレイ-緑 */
 Overlays_Green = 0x60001202,
 
 /* Overlays Blue オーバーレイ-青 */
 Overlays_Blue = 0x60001203,
 
 /* ROI Area ROI面積 */
 ROI_Area = 0x60001301,
 
 /* ROI Mean ROI平均 */
 ROI_Mean = 0x60001302,
 
 /* ROI Standard Deviation ROI標準偏差 */
 ROI_Standard_Deviation = 0x60001303,
 
 /* Overlay Data オーバーレイデータ */
 Overlay_Data = 0x60003000,
 
 /* Comments */
 Group_6000_Comments = 0x60004000,
  
 /*************************************************************************/
 /* Group Length グループ長さ */
 Group_7FE0_Length = 0x7FE00000,
 
 /* Pixel Data 画素データ */
 Pixel_Data = 0x7FE00010, 
 
 /* Item 項目 */
 Item = 0xFFFEE000,
 
 /* Item Delimitation Item 項目区切り項目 */
 Item_Delimitation_Item = 0xFFFEE00D,
 
 /* Sequence Delimitation Item シーケンス区切り項目 */
 Sequence_Delimitation_Item = 0xFFFEE0DD
} DICOM_TAG;

/* 
* Patient_Orientation		= 0x00200020 ,  患者方向
*  画像面に関連した患者方向は、正の行軸および正の列軸の解剖学的な方向を
*  示す2 つの値によって指定される。最初の登録は行の方向であり、最初の行
*  の中の最初の画素からその行の中の最後の画素の方向によって与えられる。
*  2 番目の登録は列の方向であり、最初の列の最初の画素からその列の中の最
*  後の画素の方向によって与えられる。
*  解剖学的方向は頭文字によって指定される。
*   A(anterior:前面) P(posterior:背面)
*   R(right:右) L(left:左)
*	H(head:頭)	F(foot:足)
*  方向の属性の各々の値は、少なくともこれらの文字の1つを含む。
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
 
