/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkUnstructuredGridToPolyDataFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkUnstructuredGridToPolyDataFilter.h"

#include "vtkUnstructuredGrid.h"

vtkCxxRevisionMacro(vtkUnstructuredGridToPolyDataFilter, "1.7");

//----------------------------------------------------------------------------
// Specify the input data or filter.
void vtkUnstructuredGridToPolyDataFilter::SetInput(vtkUnstructuredGrid *input)
{
  this->vtkProcessObject::SetNthInput(0, input);
}

//----------------------------------------------------------------------------
// Specify the input data or filter.
vtkUnstructuredGrid *vtkUnstructuredGridToPolyDataFilter::GetInput()
{
  if (this->NumberOfInputs < 1)
    {
    return NULL;
    }
  
  return (vtkUnstructuredGrid *)(this->Inputs[0]);
}


//----------------------------------------------------------------------------
// Copy the update information across
void vtkUnstructuredGridToPolyDataFilter::ComputeInputUpdateExtents(vtkDataObject *output)
{
  vtkDataObject *input = this->GetInput();

  if (input)
    {
    this->vtkPolyDataSource::ComputeInputUpdateExtents(output);
    input->RequestExactExtentOn();
    }
}

//----------------------------------------------------------------------------
void vtkUnstructuredGridToPolyDataFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
