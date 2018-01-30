/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float *pt, const int *dims, 
                        const float *X, const float *Y, const float *F)
{
    // Some error checking!
    //if the x point is below the first (lowest) x dimension or past the last (greatest) x dimension
    if ((pt[0] < X[0]) || (pt[0] > X[dims[0]-1])){
        return 0;
    }
    //if the given y point is below the lowest y dimension or beyond the greatest y dimension
    if ((pt[1] < Y[0]) || (pt[1] > Y[dims[1]-1])){
        return 0;
    }

    int i, j;
    int coords[2] = {0, 0}; //the logical point index representing the cell
    for (i = 0; i < dims[0]-1; i++){ //iterate through x dimensions, if the point is bracketed between two x's
        if ((X[i] <= pt[0]) && (pt[0] < X[i+1])){
            coords[0] = i; //store the x index away as the logical x
            break;
        }
    }
    for (j = 0; j < dims[1]-1; j++){ //iterate through y dimensions
        if ((Y[j] <= pt[1]) && (pt[1] < Y[j+1])){ //if the pt is between the actual Y dimensions
            coords[1] = j; //store away the logical y index
            break;
        }
    }
    // placeholders for logical point indices for lowerleft, lowerright, upperleft, upperright
    int ll[2];
    int lr[2];
    int ul[2];
    int ur[2];
    // placeholders for point indices for ll, lr, ul, ur
    int llInd, lrInd, ulInd, urInd;
    /* ll = i,j
    // lr = i+1, j
    // ul = i, j+1
    // ur = i+1, j+1 */
    ll[0] = coords[0];
    ll[1] = coords[1];
    lr[0] = coords[0]+1;
    lr[1] = coords[1];
    ul[0] = coords[0];
    ul[1] = coords[1]+1;
    ur[0] = coords[0]+1;
    ur[1] = coords[1]+1;
    llInd = GetPointIndex(ll, dims);
    lrInd = GetPointIndex(lr, dims);
    ulInd = GetPointIndex(ul, dims);
    urInd = GetPointIndex(ur, dims);
    // store away the Scalar values at each point
    float Fll = F[llInd];
    float Flr = F[lrInd];
    float Ful = F[ulInd];
    float Fur = F[urInd];
    
    // t value - make sure to calculate actual X coordinates, NOT logical point index
    float xbottomT = (pt[0]-X[ll[0]]) / (X[lr[0]]-X[ll[0]]);
    float xbottomF = Fll + xbottomT * (Flr - Fll);
    
    // t value - this is unnecessary because it's the same as xbottomT
    float xtopT = (pt[0]-X[ul[0]]) / (X[ur[0]]-X[ul[0]]);
    float xtopF = Ful + xtopT * (Fur - Ful);

    float ymidT = (pt[1]-Y[ll[1]]) / (Y[ul[1]]-Y[ll[1]]);
    float ymidF = xbottomF + ymidT * (xtopF-xbottomF);

    return ymidF;
}

// ****************************************************************************
//  Function: BoundingBoxForCell
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     cellId: a cellIndex (I.e., between 0 and GetNumberOfCells(dims))
//     bbox (output): the bounding box of cellId.  Format should be
//                     bbox[0]: the minimum X value in cellId.
//                     bbox[1]: the maximum X value in cellId.
//                     bbox[2]: the minimum Y value in cellId.
//                     bbox[3]: the maximum Y value in cellId.
//
//  Returns:  None (argument bbox is output)
//
// ****************************************************************************

void
BoundingBoxForCell(const float *X, const float *Y, const int *dims,
                   int cellId, float *bbox)
{
    int* idx = new int[2];
    GetLogicalCellIndex(idx, cellId, dims);
    int x = idx[0];
    int y = idx[1];
    // numcells is |X|-1 * |Y|-1, but they are zero indexed, so it's numCells-1
    // Error check: if the given cellId is beyond the number of cells we have
    if (cellId > (dims[0]-1)*(dims[1]-1)-1){
        bbox[0] = -100;
        bbox[1] = 100; 
        bbox[2] = -100;
        bbox[3] = 100;
    } else {
        //make sure to get actual X and Y min/max, not logical
        bbox[0] = X[x]; 
        bbox[1] = X[x + 1]; 
        bbox[2] = Y[y];
        bbox[3] = Y[y + 1];
    }
    delete [] idx;
}

// ****************************************************************************
//  Function: CountNumberOfStraddingCells
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//  Returns:  the number of cells that straddle 0, i.e., the number of cells
//            that contains points who have F>0 and also have points with F<0.
//
// ****************************************************************************

int
CountNumberOfStraddlingCells(const float *X, const float *Y, const int *dims,
                             const float *F)
{
    int numCells = GetNumberOfCells(dims);
    int totalCount = 0; // total count of cells that are straddling
    int* coords = new int[2]; //logical cell index
    int ll[2]; // the four corners of the cell
    int lr[2];
    int ul[2];
    int ur[2];
    int llInd, lrInd, ulInd, urInd; //point indices for the four corners
    for (int i = 0; i < numCells; i++){
        GetLogicalCellIndex(coords, i, dims);
        ll[0] = coords[0];
        ll[1] = coords[1];
        lr[0] = coords[0]+1;
        lr[1] = coords[1];
        ul[0] = coords[0];
        ul[1] = coords[1]+1;
        ur[0] = coords[0]+1;
        ur[1] = coords[1]+1;
        llInd = GetPointIndex(ll, dims);
        lrInd = GetPointIndex(lr, dims);
        ulInd = GetPointIndex(ul, dims);
        urInd = GetPointIndex(ur, dims);
        int posCount = 0; // how many of the 4 corners have scalar value > 0
        int negCount = 0; // how many of the 4 corners have scalar value < 0
        // if any of the corners are above zero
        if (F[llInd] > 0 || F[lrInd] > 0 || F[ulInd] > 0 || F[urInd] > 0){
            // increments to 1, which is also True
            posCount++;
        }
        // if any of the corners are below zero
        if (F[llInd] < 0 || F[lrInd] < 0 || F[ulInd] < 0 || F[urInd] < 0){
            // increments to 1, which is also True
            negCount++;
        }
        // if 1 & 1 -> True
        if (posCount && negCount){
            totalCount++;
        }
    }
    delete [] coords;
    return totalCount;
}

int main()
{
    int  i;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj2_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int numCells = CountNumberOfStraddlingCells(X, Y, dims, F);
    cerr << "The number of cells straddling zero is " << numCells << endl;

    float bbox[4];
    const int ncells = 5;
    int cellIds[ncells] = { 0, 50, 678, 1000, 1200 };
    for (i = 0 ; i < ncells ; i++)
    {
        BoundingBoxForCell(X, Y, dims, cellIds[i], bbox);
        cerr << "The bounding box for cell " << cellIds[i] << " is " 
             << bbox[0] << "->" << bbox[1] << ", " << bbox[2] << "->" << bbox[3]
             << endl;
    }

    const int npts = 10;
    float pt[npts][3] = 
         {
            {1.01119, 0.122062, 0},
            {0.862376, 1.33839, 0},
            {0.155026, 0.126123, 0},
            {0.69736, 0.0653565, 0},
            {0.2, 0.274117, 0},
            {0.893699, 1.04111, 0},
            {0.608791, -0.0533753, 0},
            {1.00543, 0.138024, 0},
            {0.384128, -0.0768977, 0},
            {0.666757, 0.60259, 0},
         };

    

    for (i = 0 ; i < npts ; i++)
    {
        float f = EvaluateFieldAtLocation(pt[i], dims, X, Y, F);
        cerr << "Evaluated field at (" << pt[i][0] <<"," << pt[i][1] << ") as "
             << f << endl;
    }
    
   
    cerr << "Infinite loop here, else Windows people may have the terminal "
         << "disappear before they see the output."
         << " Remove these lines if they annoy you." << endl;
    cerr << "(press Ctrl-C to exit program)" << endl;
    while (1) ; 
}




