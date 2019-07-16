/*********************************************************************
 * messages.h
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: messages.h 842 2006-12-04 06:21:20Z sasaoka $ */

#ifndef MESSAGES_H
#define MESSAGES_H

static STRING_TABLE ErrorStr[] = {
/*** AzLib Specific (00000-09999) ***/
/* File IO */
{ 1001, "This file is not %s file.\n"},
{ 1002, "%s file is not supported.\n"},
{ 1003, "This file format is not supported.\n"},
{ 1004, "Large file support error. (%s)\n"},
/* Solver */
{ 3001, "Pivoting Error\n"},
/* Graphic */
{ 7001, "This graphic type is not supported.\n"},
{ 7002, "%s size of this image is incorrect.\n"},
{ 7003, "Color Bit [%d] is not supported.\n"},
{ 7004, "Non-uniform resolution image is not supported.\n"},
{ 7301, "DICOM tag [%0#10lx] is not implremented.\n"},
{ 7401, "TIFF Edition Number [%d] is not supported.\n"},
{ 7402, "TIFF Compression Method [%ld] is not supported.\n"},
{ 7403, "TIFF Photometric interpretation [%ld] is not supported.\n"},
{ 7404, "Some IFD found. (Next IFD = %ld)\n"
        "Tiff image that includes some images is not supported\n"},
/* Memory Allocation */
{ 9001, "Allocation of %d [MBytes] memory was failed.\n"},

/*** Shared (10000-99999) ***/
/* File IO */
{10001, "Can not open file [%s].\n"},
{10002, "Failed processing while reading file [%s] line %d.\n"},
{10003, "Failed reading file [%s] line %d.\n"},
{10004, "Opening scratch file was failed. "
        "Please confirm SCRATCH_DIR parameter.\n"},
{10005, "Scratch file I/O failure is detected. "
        "Disk space for scratch files may be insufficient.\n"},
{10006, "File output failed. Disk space may be insufficient.\n"},
{10007, "File %s output failed. Disk space may be insufficient.\n"},
{10008, "Reading a file failed.\n"},
{10009, "Reading file %s failed.\n"},
/* Parameter File */
{11001, "Description of array ({...}) on Parameter "
        "File line %d is incorrect.\n"},
{11002, "Description of character string (\"...\") on Parameter "
        "File line %d is incorrect.\n"},
{11003, "On Parameter File line %d, %s is not defined.\n"},
{11004, "Parameter %s is not integer.\n"},
{11005, "Parameter %s is not real number.\n"},
{11006, "Parameter OBJECT_FUNC does not exist on Parameter File. "
        "It is required to define objective function or analysis type.\n"},
{11007, "Type of function of parameter %s is incorrect.\n"},
{11008, "Subcase ID refered by Parameter %s does not exist.\n"},
{11009, "Parameter %s refers Local coordinate ID %d. But it does not exist.\n"},
{11010, "Parameter %s is not character string(\"...\").\n"},
{11011, "Description of parenthesis or brace on Parameter File line %d "
        "is incorrect.\n"},
{11012, "Description of arithmetical progression (*:*:*) in Parameter File "
        "line %d is incorrect.\n"},
{11013, "In parameter %s, size of array exceeded the maximum value %d.\n"},
{11014, "Parameter %s is incorrect.\n"},
{11015, "Parameter %s refers Node ID %d. But the Node does not exist.\n"},
/* Nastran Bulk Data File */
{12001, "Maximum number of nesting levels of INCLUDE statments on "
        "Nastran Bulk Data is exceeded limit of %d.\n"},
{12002, "Nastran Bulk Data[%s] has too many fields. "
        "Maximum number of fields is %d.\n"},
{12003, "Element %d refers node ID %d. But it does not exist.\n"},
{12004, "Element %d refers property ID %d. But it does not exist.\n"},
{12005, "Property %d refers material ID %d. But it does not exist.\n"},
{12006, "Node %d refers Local coordinate ID %d. But it does not exist.\n"},
{12007, "%s refers node ID %d. But it does not exist.\n"},
{12008, "Description of Nastran Bulk Data [%s] is incorrect.\n"},
{12009, "There is no BEGIN BULK line of Nastran Bulk Data file.\n"},
{12010, "Local coordinate %d is defined by coordinate %d. "
        "But the coordinate does not exist. Or it is not supported.\n"},
{12011, "Some elements have same ID %d.\n"},
{12012, "In Subcase %d, Direct Matrix Input DMIG:%s refers node ID %d. "
        "But the node does not exist.\n"},
{12013, "Direct Matrix Input DMIG:%s on scalar point is not supported.\n"},
{12014, "Some nodes have same ID %d.\n"},
{12015, "Shell Element %d refers property ID %d. "
        "But it is not PSHELL format.\n"},
{12016, "Solid Element %d refers property ID %d. "
        "But it is not PSOLID format.\n"},
{13001, "Writing scratch file failed. Please confirm free area of the disk.\n"},
{13002, "Reading scratch file failed.\n"},
/* Solver */
{30001, "%d [MBytes] memory is required to execute of sparse/eigen solver. "
        "Please increase the value of MEM_SIZE.\n"},
{30002, "Failed in preprocessing of CG Method. Problem may be in condition "
        "of constraint.\n"},
{30003, "Failed in processing of boundary surface of the model. "
        "Element connectivity may be incorrect.\n"},
{30004, "Specification of the independent DOF of RBAR Element %d "
        "is incorrect.\n"},
{30005, "Fatal error of sparse/eigen solver occurred. Error message : %s\n"},
{30006, "Some master/slave DOFs of MPC or rigid element are incorrect. "
        "Please check Node ID %d DOF %d.\n"},
{30007, "Eigen mode tracking failed.\n"},
{30008, "Any independent DOFs are not exist in this analysis. "
        "Please confirm the condition of this analysis.\n"},
{40001, "Starting a child module for parallel processing failed. "
        "Please confirm the executable file %s and "
        "PATH environment variable.\n"},
/* Memory Allocation */
{90001, "Allocation of %d [MBytes] memory was failed. "
        "Please check the value of MEM_SIZE in Parameter File.\n"},
{90002, "Available memory is insufficient. "
       	"Please check the value of MEM_SIZE in Parameter File.\n"},
{90003, "Available memory is insufficient under the 32bit environment.\n"},
{99001, "License error is detected. (Error ID = %d)\n"},
{99002, "A license server does not respond. "
        "Please check a setup of your license server or condition of "
        "your network, etc.\n"},
{99003, "There is no license of this program. "
        "It may be that the license has expired or the program " 
        "under execution has reached the number of licenses.\n"},
{99004, "There is no license for 64bit version.(Error ID = %d)\n"},
{99005, "There is no license for SE version.(Error ID = %d)\n"},
/* Unknown Error */
{99999, "Internal error is detectd.\n"}
};

static STRING_TABLE WarningStr[] = {
/*** AzLib Specific (00000-09999) ***/
/* Solver */
{ 3001, "This element has Negative volume."
        "(element ID [%d], Gauss Point [%d])\n"},
{ 3002, "Unused node [%d] is detected.\n"},
{ 3003, "CG iteration was NOT converged!!\n"},
{ 3004, "MATRIX[%d-%d]: %15.7e -> %15.7e\n"}, /* nonzero_cg.c : fix please!!*/

/*** Shared (10000-99999) ***/
/* Parameter File */
{11001, "While processing the Parameter File line %d, overflow is occurred. "
        "The number is processed with the real number.\n"},
{11002, "In Parameter File of the %d line, Subsutitution for %s is invalid.\n"},
{11003, "Invalid sentence on Parameter File line %d id detected. [%s] "
        "This sentence is egnored.\n"},
{11004, "There is a parameter which is not substituted in Parameter File of "
        "the %d line.\n"},
{11005, "Parameter[%s] is over setting range(%d - %d). "
        "It changes into %d.\n"},
{11006, "Parameter[%s]is over setting range(%11.3e - %11.3e). "
        "It changes into %11.3e.\n"},
{11007, "The Parameter CONST_FUNC(Constraint Function) is not defined.\n"},
{11008, "The element ID %d was  specified on parameter DISABLE_ELEMENT(SUBCASE "
        "ID %d). But it does not exist.\n"},
{11009, "In DISABLE_ELEMENT(SUBCASE ID %d), the element ID %s was specified. "
	    "But it cannot be disabled because it is in design domain.\n"},
{11010, "Parameter OPTIMIZE_TYPE = ANALYSIS is specified. "
        "But maximum iteration of optimization (IT_MAX) is not zero. "
        "Optimization process is not performed.\n"},
{11011, "Parameter [%s = %d] is incorrect. It changes into %d.\n"},
{11012, "Parameter %s refers element property ID %d. "
        "But it does not exist.\n"},
{11013, "Parameter %s refers element property ID %d. This property ID must be "
        "included into NONDESIGN_PROP parameter. This property ID will be "
        "added to NONDESIGN_PROP parameter.\n"},
{11014, "Unsupported Objective/Constraint function for shell elements was "
        "specified. The shell elements will be specified as "
        "non-design domain.\n"},
{11015, "Parameter %s will be ignored.\n"},
{11016, "Parameter NORMAL_FACT_MAX and NORMAL_FACT_MIN will be corrected.\n"},
{11017, "Parameter NORMAL_FACT_INIT and NORMAL_FACT_MIN will be corrected.\n"},
{11018, "Parameter NORMAL_FACT_INIT and NORMAL_FACT_MAX will be corrected.\n"},
{11019, "Parameter MEM_SIZE is limited to 1800 under the 32bit environment.\n"},
/* Nastran Bulk Data File */
{12001, "Bulk Data FIle [%s] line %d is regarded as free field format.\n"},
{12002, "%d inputs of Bulk Data are ignored. "
        "Because they are not supported by this program. "
        "Please refer to file [%s] to check the ignored Bulk Data input.\n"},
{12003, "Element ID %d is refered by PLOAD4 Bulk Data. "
        "But the element ID is not exist.\n"},
{12004, "Multiple Bulk Data Input [%s] in Bulk Data File is detected. "
        "Only the last input is effective.\n"},
{12005, "CTETRA ID %d was modified the connectivity.\n"},
{12006, "Some mid nodes of a quadratic element ID %d are not described. "
        "This element is treated as a linear element.\n"},
{12007, "In the Bulk Data Input %s, field %s(%d) is incorrect.\n"},
{12008, "Bulk Data FIle [%s] line %d is ignored.\n"},
/* Pre-Post IO */
{13001, "The data type is ignored. [internal number is %d]\n"},
{13002, "The thickness of design plate refered by property %d is not "
        "defined.\n"},
{13003, "Although angle acceleration is specified by RFORCE, it is not "
        "supported\n"},
{13004, "Maximum six number of master nodes supported "
        "by MPC\n"},
{13005, "The thickness of base plate is over the thickness of "
        "the plate itself. The thickness of base plate is set as 0.\n"},
{13006, "The data type which is not supported in the output of "
        "Neutral File of FEMAP is ignored. [internal number %d]\n"},
{13007, "The data type which is not supported in the output of "
        "Universal File of IDEAS is ignored. [internal number %d]\n"},
{13009, "Constraint value of node %d set ID %d will be ignored for writing "
        "FEMAP Neutral File\n"},
{13010, "On writing FEMAP Neutral File, six master D.O.Fs of RBAR element "
        "must be established on same node. On RBAR element %d, node %d will be "
        "a master node.\n"},
/* Solver */
{30001, "Jacobian of the element ID %d is 0 or negative. "
        "Please check distortion of element."},
{30002, "Subspace iteration of eigenvalue analysis was not converged. "
        "There is the possibility with the inadequate accuracy of "
        "eigenvalue and eigenvector. "
        "It may be that a problem is in the stability of a model.\n"},
{30003, "Iteration of CG mthoed was not counverged. There is the possibility "
        "with the inadequate accuracy of displacement. "
        "It may be that a problem is in the constraint conditions "
        "of the model.\n"},
{30004, "The similar mode was not detected by mode tracking. "
        "Please increase the value of parameter EIGEN_MODE_INC.\n"},
{30005, "Since available memory is insufficient, frequent disk I/O may couse "
        "slow performance. If possible, please increase the value of MEM_SIZE "
        "parameter more %d[MB] to improve the performance. "
        "If you perform static analysis, you can use Iterative Solver "
        "(SOLVER_TYPE = ICCG) which requires small memory.\n"},
{30006, "Since available memory is insufficient, slow performance Iterative "
        "Solver was selected. If possible, please increase the value of "
        "MEM_SIZE parameter more %d[MB] to improve the performance.\n"},
{30007, "Slave D.O.F (Node ID: %d, D.O.F ID: %d) of MPC of Rigid Element "
        "is constrainted. The condition will be ignored.\n"},
{30008, "Node ID: %d D.O.F ID: %d is configured as slave D.O.F of multiple "
        "MPCs or Rigid elements.\n"},
{30009, "Eigen mode %d and %d are regarded as duplicate root. The difference "
        "of the eigen value is %15.7e.\n"},
{30010, "Illegal property ID or material ID is specified for the element %d. "
        "This element is ignored for the analysis.\n"},
{30011, "Since available memory is insufficient, frequent disk I/O may couse "
        "slow performance. If possible, please increase the value of MEM_SIZE "
        "parameter more %d[MB] to improve the performance.\n"},
{30012, "Because of EXTRA_EIGEN_MODE_INC parameter, the eingen value analysis "
        "will be restarted.\n"},
{30013, "The shell elements are not supported in the thermal "
        "conduction analysis.\n"},
{30014, "The beam elements are not supported in the thermal "
        "conduction analysis.\n"},
/* Error and Warning */
{90001, "Since there is too much ERROR %5.5d, The error message output will "
        "be omitted.\n"},
{90002, "Since there is too much WARNING %5.5d, The warning message output "
        "will be omitted.\n"},
{99001, "The License will expire within 15 days.\n"},
{99002, "The License will expire within 2 days.\n"},
{99003, "The license server connection failed. \n"},
{99004, "Many reconnection with the license server are occured "
        "in a relatively short period. \n"},
{99005, "The location of the license file or server is not configured. "
        "The default location (%s) will be used.\n"}
};

#endif /* MESSAGES_H */

