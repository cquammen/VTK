/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPolyDataToPolyDataFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPolyDataToPolyDataFilter.h"

#include "vtkPolyData.h"

vtkCxxRevisionMacro(vtkPolyDataToPolyDataFilter, "1.17");

//----------------------------------------------------------------------------
vtkPolyDataToPolyDataFilter::vtkPolyDataToPolyDataFilter() 
{
  this->NumberOfRequiredInputs = 1;
}
//----------------------------------------------------------------------------
// Specify the input data or filter.
void vtkPolyDataToPolyDataFilter::SetInput(vtkPolyData *input)
{
  this->vtkProcessObject::SetNthInput(0, input);
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
vtkPolyData *vtkPolyDataToPolyDataFilter::GetInput()
{
  if (this->NumberOfInputs < 1)
    {
    return NULL;
    }
  
  return (vtkPolyData *)(this->Inputs[0]);
}

//----------------------------------------------------------------------------
void vtkPolyDataToPolyDataFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
