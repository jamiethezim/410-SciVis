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
#include <vtkPNGWriter.h>

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

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
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


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}



// *************************
// helper function
// args:
// cellId - the cell you want to classify
// int isovalue - the deciding value that determines presence of point values
// dims - int* the dimensions of the grid
//
// returns the point indices of the four corners of the cell

void GetFourCorners(int* four_corners, int cellId, const int* dims){
    int coords[2];
    GetLogicalCellIndex(coords, cellId, dims);

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

    four_corners[0] = llInd;
    four_corners[1] = lrInd;
    four_corners[2] = ulInd;
    four_corners[3] = urInd;

}

// ***********************************
// helper function
// identifies which case the cell isosurface is
// from case 0 to case 15
// args:
//   cellId - the cell we want to classify
//   isovalue - the breaking point for determining if the corner is > or <
//   dims - int array - dims[0] is number of x points, dims[1] number of y points. 
//         needed for getting four corners of cell
//   F - the float values at every point
// ***********************************
int IdentifyCase(int cellId, float isovalue, const int* dims, const float* F){
    int corners[4];
    GetFourCorners(corners, cellId, dims);
    // now corners[0] is pointidx of lower left
    // corners[1] is pointidx of lr
    // corners[2] is pointidx of ul
    // corners[3] is pointidx of ur
    int classification = 0x0;
    // ur is V3
    if (F[corners[3]] > isovalue){
        classification |= 0x8;
    }
    // ul is V2
    if (F[corners[2]] > isovalue){
        classification |= 0x4;
    }
    // lr is V1
    if (F[corners[1]] > isovalue){
        classification |= 0x2;
    }
    // ll is V0
    if (F[corners[0]] > isovalue){
        classification |= 0x1;
    }
    return classification;
}

// ******************************************
// CalculatePosition
// given FX, FA, FB, B, & A, we can find the position of the point! point X
float CalculatePosition(float FX, float FA, float FB, float A, float B){
    float res = 0.0;
    float FS = (FX - FA)/(FB - FA);
    float t = B - A;
    res = FS * t + A;
    return res;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    

    fprintf(stdout, "min_x is %f, max_x is %f, min_y is %f, max_y is %f\n", X[0], X[dims[0]-1], Y[0], Y[dims[1]-1]);
    int numSegments[16] = {0, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 0};
    int lup[16][4];
    lup[0][0] = lup[0][1] = lup[0][2] = lup[0][3] = -1;
    lup[1][0] = 0; lup[1][1] = 3; lup[1][2] = lup[1][3] = -1;
    lup[2][0] = 0; lup[2][1] = 1; lup[2][2] = lup[2][3] = -1;
    lup[3][0] = 1; lup[3][1] = 3; lup[3][2] = lup[3][3] = -1;
    lup[4][0] = 2; lup[4][1] = 3; lup[4][2] = lup[4][3] = -1;
    lup[5][0] = 0; lup[5][1] = 2; lup[5][2] = lup[5][3] = -1;
    lup[6][0] = 0; lup[6][1] = 1; lup[6][2] = 2; lup[6][3] = 3;
    lup[7][0] = 1; lup[7][1] = 2; lup[7][2] = lup[7][3] = -1;
    lup[8][0] = 1; lup[8][1] = 2; lup[8][2] = lup[8][3] = -1;
    lup[9][0] = 0; lup[9][1] = 3; lup[9][2] = 1; lup[9][3] = 2;
    lup[10][0] = 0; lup[10][1] = 2; lup[10][2] = lup[10][3] = -1;
    lup[11][0] = 2; lup[11][1] = 3; lup[11][2] = lup[11][3] = -1;
    lup[12][0] = 1; lup[12][1] = 3; lup[12][2] = lup[12][3] = -1;
    lup[13][0] = 0; lup[13][1] = 1; lup[13][2] = lup[13][3] = -1;
    lup[14][0] = 0; lup[14][1] = 3; lup[14][2] = lup[14][3] = -1;
    lup[15][0] = lup[15][1] = lup[15][2] = lup[15][3] = -1;

    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

    float IsoVal = 3.2;
    int icase;
    int nsegments;
    int fourcorns[4] = {0, 0, 0, 0};
    int llInd;
    int lrInd;
    int urInd;
    int ulInd;
    int ll[2] = {0, 0};
    int lr[2] = {0, 0};
    int ul[2] = {0, 0};
    int ur[2] = {0, 0};
    float Fll = 0.0;
    float Flr = 0.0;
    float Ful = 0.0;
    float Fur = 0.0;
    int edge1;
    int edge2;
    float pt1[2] = {0.0, 0.0};
    float pt2[2] = {0.0, 0.0};
    //iterate through cells
    fprintf(stdout, "num cells is %d\n", GetNumberOfCells(dims));
    for (i = 0; i < GetNumberOfCells(dims); i++){
        icase = IdentifyCase(i, IsoVal, dims, F);
        nsegments = numSegments[icase];
        fprintf(stdout, "cell %d is case %d and has %d segments and \n", i, icase, nsegments);
        for (j = 0; j < nsegments; j++){
            edge1 = lup[icase][2*j];
            edge2 = lup[icase][2*j+1];
            fprintf(stdout, "\thas edge1 at %d and edge2 at %d\n", edge1, edge2);
            GetFourCorners(fourcorns, i, dims);
            llInd = fourcorns[0];
            lrInd = fourcorns[1];
            ulInd = fourcorns[2];
            urInd = fourcorns[3];
            GetLogicalPointIndex(ll, llInd, dims);
            GetLogicalPointIndex(lr, lrInd, dims);
            GetLogicalPointIndex(ul, ulInd, dims);
            GetLogicalPointIndex(ur, urInd, dims);
            Fll = F[llInd];
            Flr = F[lrInd];
            Ful = F[ulInd];
            Fur = F[urInd];
            fprintf(stdout, "ll is %f,%f has val %f\n", X[ll[0]], Y[ll[1]], Fll);
            fprintf(stdout, "lr is %f,%f has val %f\n", X[lr[0]], Y[lr[1]], Flr);
            fprintf(stdout, "ul is %f,%f has val %f\n", X[ul[0]], Y[ul[1]], Ful);
            fprintf(stdout, "ur is %f,%f has val %f\n", X[ur[0]], Y[ur[1]], Fur);
            if (edge1 == 0){
                //interpolate ll -> lr
                pt1[0] = CalculatePosition(IsoVal, Fll, Flr, X[ll[0]], X[lr[0]]);
                pt1[1] = Y[ll[1]];
            }
            else if (edge1 == 1){
                //interpolate lr -> ur
                pt1[0] = X[lr[0]];
                pt1[1] = CalculatePosition(IsoVal, Fur, Flr, Y[lr[1]], Y[ur[1]]); // ****** comment this!
            }
            else if (edge1 == 2){
                //interpolate ul -> ur
                pt1[0] = Y[ul[1]];
                pt1[1] = CalculatePosition(IsoVal, Ful, Fur, X[ul[0]], X[ur[0]]);
            }
            else if (edge1 == 3){
                //interpolate ll -> ul
                pt1[0] = X[ll[0]];
                pt1[1] = CalculatePosition(IsoVal, F[ulInd], F[llInd], Y[ll[1]], Y[ul[1]]); //*** comment this!
            }
            if (edge2 == 0){
                //interpolate ll -> lr
                pt2[0] = CalculatePosition(IsoVal, F[llInd], F[lrInd], X[ll[0]], X[lr[0]]);
                pt2[1] = Y[ll[1]];
            }
            else if (edge2 == 1){
                //interpolate lr -> ur
                pt2[0] = X[lr[0]];
                pt2[1] = CalculatePosition(IsoVal, F[urInd], F[lrInd], Y[lr[1]], Y[ur[1]]); //*comment this!!!
            }
            else if (edge2 == 2){
                //interpolate ul -> ur
                pt2[0] = Y[ul[1]];
                pt2[1] = CalculatePosition(IsoVal, F[ulInd], F[urInd], X[ul[0]], X[ur[0]]);
            }
            else if (edge2 == 3){
                //interpolate ll -> ul
                pt2[0] = X[ll[0]];
                pt2[1] = CalculatePosition(IsoVal, F[ulInd], F[llInd], Y[ll[1]], Y[ul[1]]); ///**** comment this!
            }
            fprintf(stdout, "pt1x %f,%f, pt2x %f,%f\n", pt1[0], pt1[1], pt2[0], pt2[1]);
            sl.AddSegment(pt1[0], pt1[1], pt2[0], pt2[1]);
        }
    }
    fprintf(stdout, "i'm done!\n");
// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!

    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
