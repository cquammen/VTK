/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageCursor3D.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageCursor3D.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageCursor3D, "1.17");
vtkStandardNewMacro(vtkImageCursor3D);

//----------------------------------------------------------------------------
void vtkImageCursor3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  
  int idx;
  
  os << indent << "Cursor Radius: " << this->CursorRadius << "\n";
  os << indent << "Cursor Value: " << this->CursorValue << "\n";
  os << indent << "Cursor Position: (" << this->CursorPosition[0];
  for (idx = 1; idx < 3; ++idx)
    {
    os << ", " << this->CursorPosition[idx];
    }
  os << ")\n";  
}

//----------------------------------------------------------------------------
vtkImageCursor3D::vtkImageCursor3D()
{
  this->CursorPosition[0] = 0;
  this->CursorPosition[1] = 0;
  this->CursorPosition[2] = 0;
  
  this->CursorRadius = 5;
  this->CursorValue = 255;
}



template <class T>
void vtkImageCursor3DExecute(vtkImageCursor3D *self,
                             vtkImageData *outData, T *ptr)
{
  int min0, max0, min1, max1, min2, max2;
  int c0, c1, c2;
  int idx;
  float value = 0.0;
  int rad = self->GetCursorRadius();
  
  c0 = (int)(self->GetCursorPosition()[0]);
  c1 = (int)(self->GetCursorPosition()[1]);
  c2 = (int)(self->GetCursorPosition()[2]);
  value = self->GetCursorValue();
  
  outData->GetExtent(min0, max0, min1, max1, min2, max2);
  
  if (c1 >= min1 && c1 <= max1 && c2 >= min2 && c2 <= max2)
    {
    for (idx = c0 - rad; idx <= c0 + rad; ++idx)
      {
      if (idx >= min0 && idx <= max0)
        {
        ptr = (T *)(outData->GetScalarPointer(idx, c1, c2));
        *ptr = (T)(value);
        }
      }
    }
  
  
  if (c0 >= min0 && c0 <= max0 && c2 >= min2 && c2 <= max2)
    {
    for (idx = c1 - rad; idx <= c1 + rad; ++idx)
      {
      if (idx >= min1 && idx <= max1)
        {
        ptr = (T *)(outData->GetScalarPointer(c0, idx, c2));
        *ptr = (T)(value);
        }
      }
    }
  
  
  if (c0 >= min0 && c0 <= max0 && c1 >= min1 && c1 <= max1)
    {
    for (idx = c2 - rad; idx <= c2 + rad; ++idx)
      {
      if (idx >= min2 && idx <= max2)
        {
        ptr = (T *)(outData->GetScalarPointer(c0, c1, idx));
        *ptr = (T)(value);
        }
      }
    }
}

//----------------------------------------------------------------------------
// Split up into finished and border datas.  Fill the border datas.
void vtkImageCursor3D::ExecuteData(vtkDataObject *out)
{
  void *ptr = NULL;
  
  // let superclass allocate data
  this->vtkImageInPlaceFilter::ExecuteData(out);

  vtkImageData *outData = this->GetOutput();
  
  switch (outData->GetScalarType())
    {
    vtkTemplateMacro3(vtkImageCursor3DExecute, this, 
                      outData, (VTK_TT *)(ptr));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }
  
}




