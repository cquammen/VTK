/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkHull.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkHull - produce an n-sided convex hull
// .SECTION Description
// vtkHull is a filter which will produce an n-sided convex hull given a
// set of n planes. (The convex hull bounds the input polygonal data.)
// The hull is generated by squeezing the planes towards the input
// vtkPolyData, until the planes just touch the vtkPolyData. Then,
// the resulting planes are used to generate a polyhedron (i.e., hull)
// that is represented by triangles.
//
// The n planes can be defined in a number of ways including 1) manually 
// specifying each plane; 2) choosing the six face planes of the input's
// bounding box; 3) choosing the eight vertex planes of the input's
// bounding box; 4) choosing the twelve edge planes of the input's
// bounding box; and/or 5) using a recursively subdivided octahedron.
// Note that when specifying planes, the plane normals should point
// outside of the convex region.
//
// The output of this filter can be used in combination with vtkLODActor 
// to represent a levels-of-detail in the LOD hierarchy. Another use of
// this class is to manually specify the planes, and then generate the
// polyhedron from the planes (without squeezing the planes towards the
// input). The method GenerateHull() is used to do this.

#ifndef __vtkHull_h
#define __vtkHull_h

#include "vtkPolyDataToPolyDataFilter.h"

class vtkCellArray;
class vtkPlanes;
class vtkPoints;
class vtkPolyData;

class VTK_GRAPHICS_EXPORT vtkHull : public vtkPolyDataToPolyDataFilter
{
public:
  static vtkHull *New();
  vtkTypeRevisionMacro(vtkHull,vtkPolyDataToPolyDataFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Remove all planes from the current set of planes.  
  void RemoveAllPlanes( void );

  // Description:
  // Add a plane to the current set of planes. It will be added at the
  // end of the list, and an index that can later be used to set this
  // plane's normal will be returned. The values A, B, C are from the
  // plane equation Ax + By + Cz + D = 0. This vector does not have to
  // have unit length (but it must have a non-zero length!). If a value
  // 0 > i >= -NumberOfPlanes is returned, then the plane is parallel
  // with a previously inserted plane, and |-i-1| is the index of the
  // plane that was previously inserted. If a value i < -NumberOfPlanes
  // is returned, then the plane normal is zero length.
  int  AddPlane( float A, float B, float C );
  int  AddPlane( float plane[3] );

  // Description:
  // Set the normal values for plane i. This is a plane that was already
  // added to the current set of planes with AddPlane(), and is now being
  // modified. The values A, B, C are from the plane equation 
  // Ax + By + Cz + D = 0. This vector does not have to have unit length.
  // Note that D is set to zero, except in the case of the method taking
  // a vtkPlanes* argument, where it is set to the D value defined there.
  void SetPlane( int i, float A, float B, float C );
  void SetPlane( int i, float plane[3] );

  // Description:
  // Variations of AddPlane()/SetPlane() that allow D to be set. These 
  // methods are used when GenerateHull() is used.
  int AddPlane( float A, float B, float C, float D );
  int AddPlane( float plane[3], float D );
  void SetPlane( int i, float A, float B, float C, float D );
  void SetPlane( int i, float plane[3], float D );

  // Description:
  // Set all the planes at once using a vtkPlanes implicit function.
  // This also sets the D value, so it can be used with GenerateHull().
  void SetPlanes( vtkPlanes *planes );

  // Description:
  // Get the number of planes in the current set of planes.
  vtkGetMacro( NumberOfPlanes, int );
  
  // Description:
  // Add the 8 planes that represent the vertices of a cube - the combination
  // of the three face planes connecting to a vertex - (1,1,1), (1,1,-1),
  // (1,-1,1), (1,-1,1), (-1,1,1), (-1,1,-1), (-1,-1,1), (-1,-1-1).
  void AddCubeVertexPlanes();

  // Description:
  // Add the 12 planes that represent the edges of a cube - halfway between
  // the two connecting face planes - (1,1,0), (-1,-1,0), (-1,1,0), (1,-1,0),
  // (0,1,1), (0,-1,-1), (0,1,-1), (0,-1,1), (1,0,1), (-1,0,-1),
  // (1,0,-1), (-1,0,1)
  void AddCubeEdgePlanes();

  // Description:
  // Add the six planes that make up the faces of a cube - (1,0,0),
  // (-1, 0, 0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)
  void AddCubeFacePlanes();

  // Description:
  // Add the planes that represent the normals of the vertices of a
  // polygonal sphere formed by recursively subdividing the triangles
  // in an octahedron.  Each triangle is subdivided by connecting the
  // midpoints of the edges thus forming 4 smaller triangles. The
  // level indicates how many subdivisions to do with a level of 0
  // used to add the 6 planes from the original octahedron, level 1
  // will add 18 planes, and so on.
  void AddRecursiveSpherePlanes( int level );

  // Description:
  // A special method that is used to generate a polyhedron directly
  // from a set of n planes. The planes that are supplied by the user
  // are not squeezed towards the input data (in fact the user need
  // not specify an input). To use this method, you must provide an
  // instance of vtkPolyData into which the points and cells defining
  // the polyhedron are placed. You must also provide a bounding box
  // where you expect the resulting polyhedron to lie. This can be
  // a very generous fit, it's only used to create the initial polygons
  // that are eventually clipped.
  void GenerateHull(vtkPolyData *pd, float *bounds);
  void GenerateHull(vtkPolyData *pd, float xmin, float xmax,
                    float ymin, float ymax, float zmin, float zmax);

protected:
  vtkHull();
  ~vtkHull();

  // The planes - 4 doubles per plane for A, B, C, D
  double     *Planes;

  // This indicates the current size (in planes - 4*sizeof(float)) of 
  // the this->Planes array. Planes are allocated in chunks so that the
  // array does not need to be reallocated every time a new plane is added
  int       PlanesStorageSize;

  // The number of planes that have been added
  int       NumberOfPlanes;

  // Internal method used to find the position of each plane
  void      ComputePlaneDistances();

  // Internal method used to create the actual polygons from the set 
  // of planes
  void      ClipPolygonsFromPlanes( vtkPoints *points, vtkCellArray *polys,
                                    float *bounds );

  // Internal method used to create the initial "big" polygon from the
  // plane equation. This polygon is clipped by all other planes to form
  // the final polygon (or it may be clipped entirely)
  void      CreateInitialPolygon( double *, int, float * );

  // The method that does it all...
  void      Execute();
private:
  vtkHull(const vtkHull&);  // Not implemented.
  void operator=(const vtkHull&);  // Not implemented.
};

#endif
