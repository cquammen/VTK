/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRectilinearWipe.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageRectilinearWipe - make a rectilinear combination of two images.
// .SECTION Description
// vtkImageRectilinearWipe makes a rectilinear combination of two
// images. The two input images must correspond in size, scalar type and
// number of components.
// The resulting image has four possible configurations
// called:
//   Quad - alternate input 0 and input 1 horizontally and
//     vertically. Select this with SetWipeModeToQuad. The Position
//     specifies the location of the quad intersection.
//   Corner - 3 of one input and 1 of the other. Select the location of
//     input 0 with with SetWipeModeToLowerLeft, SetWipeModeToLowerRight,
//     SetWipeModeToUpperLeft and SetWipeModeToUpperRight. The Position
//     selects the location of the corner.
//   Horizontal - alternate input 0 and input 1 with a vertical
//     split. Select this with SetWipeModeToHorizontal. Position[0]
//     specifies the location of the vertical transition between input 0
//     and input 1.
//   Vertical - alternate input 0 and input 1 with a horizontal
//     split. Only the y The intersection point of the rectilinear points
//     is controlled with the Point ivar.

// .SECTION Thanks
// This work was supported by PHS Research Grant No. 1 P41 RR13218-01
// from the National Center for Research Resources.

// .SECTION See Also
// vtkImageCheckerboard

#ifndef __vtkImageRectilinearWipe_h
#define __vtkImageRectilinearWipe_h

#include "vtkImageTwoInputFilter.h"

#define VTK_WIPE_QUAD 0
#define VTK_WIPE_HORIZONTAL 1
#define VTK_WIPE_VERTICAL 2
#define VTK_WIPE_LOWER_LEFT 3
#define VTK_WIPE_LOWER_RIGHT 4
#define VTK_WIPE_UPPER_LEFT 5
#define VTK_WIPE_UPPER_RIGHT 6

class VTK_IMAGING_EXPORT vtkImageRectilinearWipe : public vtkImageTwoInputFilter
{
public:
  static vtkImageRectilinearWipe *New();
  vtkTypeRevisionMacro(vtkImageRectilinearWipe,vtkImageTwoInputFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get the location of the image transition.
  vtkSetVector2Macro(Position,int);
  vtkGetVectorMacro(Position,int,2);

  // Description:
  // Specify the wipe mode. This mode determnis how input 0 and input
  // 1 are combined to produce the output. Each mode uses one or both
  // of the values stored in Position.
  //   SetWipeToQuad - alternate input 0 and input 1 horizontally and
  //     vertically. The Position specifies the location of the quad
  //     intersection.
  //   SetWipeToLowerLeft{LowerRight,UpperLeft.UpperRight} - 3 of one
  //     input and 1 of the other. Select the location of input 0 to the
  //     LowerLeft{LowerRight,UpperLeft,UpperRight}. Position
  //     selects the location of the corner.
  //   SetWipeToHorizontal - alternate input 0 and input 1 with a vertical
  //     split. Position[0] specifies the location of the vertical
  //     transition between input 0 and input 1.
  //   SetWipeToVertical - alternate input 0 and input 1 with a
  //     horizontal split. Position[1] specifies the location of the
  //     horizonal transition between input 0 and input 1.
  vtkSetClampMacro(Wipe,int,
                   VTK_WIPE_QUAD,VTK_WIPE_UPPER_RIGHT);
  vtkGetMacro(Wipe,int);
  void SetWipeToQuad()
    {this->SetWipe(VTK_WIPE_QUAD);}
  void SetWipeToHorizontal()
    {this->SetWipe(VTK_WIPE_HORIZONTAL);}
  void SetWipeToVertical()
    {this->SetWipe(VTK_WIPE_VERTICAL);}
  void SetWipeToLowerLeft()
    {this->SetWipe(VTK_WIPE_LOWER_LEFT);}
  void SetWipeToLowerRight()
    {this->SetWipe(VTK_WIPE_LOWER_RIGHT);}
  void SetWipeToUpperLeft()
    {this->SetWipe(VTK_WIPE_UPPER_LEFT);}
  void SetWipeToUpperRight()
    {this->SetWipe(VTK_WIPE_UPPER_RIGHT);}

protected:
  vtkImageRectilinearWipe();
  ~vtkImageRectilinearWipe() {};

  void ThreadedExecute(vtkImageData **inDatas, vtkImageData *outData,
                       int extent[6], int id);
  int Position[2];
  int Wipe;
private:
  vtkImageRectilinearWipe(const vtkImageRectilinearWipe&);  // Not implemented.
  void operator=(const vtkImageRectilinearWipe&);  // Not implemented.
};

#endif
