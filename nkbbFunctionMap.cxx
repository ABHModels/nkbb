//
// Code auto-generated by initpackage (XSPEC12 local 
// model package code generator).  Do not edit.
// Package: nkbb
// Function body: nkbbFunctionMap.cxx

#include    "nkbbFunctionMap.h"

#include    <XSFunctions/Utilities/XSModelFunction.h>

void 
createnkbbFunctionMap()
{


	XSFunctionMap["nkbb"]     = new XSCall<xsccCall>(lmodnkbb);

}
