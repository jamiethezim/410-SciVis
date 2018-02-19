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

#include "TriangleList.h"

static int Edges[12][2] = 
    {{0, 1}, {1, 3}, {2, 3}, {0, 2},
     {4, 5}, {5, 7}, {6, 7}, {4, 6},
     {0, 4}, {1, 5}, {2, 6}, {3, 7} };

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
    return dims[0]*dims[1]*dims[2];
    // 2D
   // return dims[0]*dims[1];
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
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
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
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
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
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
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
    idx[0] = pointId%dim[0];
    idx[1] = (pointId/dims[0])%dims[1];
    idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
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
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}




// *************************
// helper function
// args:
// cellId - the cell you want to classify
// int isovalue - the deciding value that determines presence of point values
// dims - int* the dimensions of the grid
//
// returns the point indices of the eight vertices of the cell

void GetEightVertices(int* eight_vertices, int cellId, const int* dims){
    int coords[3]; //x, y, z of vertex 0
    GetLogicalCellIndex(coords, cellId, dims);

    //logical point index
    //coords is v0
    int v0[3];
    int v1[3];
    int v2[3];
    int v3[3];
    int v4[3];
    int v5[3];
    int v6[3];
    int v7[3];

    //coords [0] is x, coords[1] is y, coords[2] is z
    /* v0 is coords[0], coords[1], coords[2]
    v1 is coords[0]+1, coords[1], coords[2]
    v2 is coords[0], coords[1]+1, coords[2]
    v3 is coords[0]+1, coords[1]+1, coords[2]
    v4 is coords[0], coords[1], coords[2]+1
    v5 is coords[0]+1, coords[1], coords[2]+1
    v6 is coords[0], coords[1]+1, coords[2]+1
    v7 is coords[0]+1, coords[1]+1, coords[2]+1
    */
    v0 = {coords[0], coords[1], coords[2]};
    v1 = {coords[0]+1, coords[1], coords[2]};
    v2 = {coords[0], coords[1]+1, coords[2]};
    v3 = {coords[0]+1, coords[1]+1, coords[2]};
    v4 = {coords[0], coords[1], coords[2]+1};
    v5 = {coords[0]+1, coords[1], coords[2]+1};
    v6 = {coords[0], coords[1]+1, coords[2]+1};
    v7 = {coords[0]+1, coords[1]+1, coords[2]+1};

    eight_vertices[0] = GetPointIndex(v0, dims);
    eight_vertices[1] = GetPointIndex(v1, dims);
    eight_vertices[2] = GetPointIndex(v2, dims);
    eight_vertices[3] = GetPointIndex(v3, dims);
    eight_vertices[4] = GetPointIndex(v4, dims);
    eight_vertices[5] = GetPointIndex(v5, dims);
    eight_vertices[6] = GetPointIndex(v6, dims);
    eight_vertices[7] = GetPointIndex(v7, dims);

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
    int vertices[8];
    GetEightVertices(vertices, cellId, dims);
    int classification = 0x0;
    if (F[vertices[7]] > isovalue){
        classification |= 0x80; //128 in binary is 1000,0000 is 0x80 in hex
    }
    if (F[vertices[6]] > isovalue){
        classification |= 0x40; //64 in binary is 0100,0000 is 0x40 in hex
    }
    if (F[vertices[5]] > isovalue){
        classification |= 0x20; //32 in binary is 0010,0000 is 0x20 in hex
    }
    if (F[vertices[4]] > isovalue){
        classification |= 0x10; //16
    }
    if (F[vertices[3]] > isovalue){
        classification |= 0x8;
    }
    if (F[vertices[2]] > isovalue){
        classification |= 0x4;
    }
    if (F[vertices[1]] > isovalue){
        classification |= 0x2;
    }
    if (F[vertices[0]] > isovalue){
        classification |= 0x1;
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


// helper function
// tri is an input to the function that will contain the return value
// returns 3 vertices composed of 3 endpoints
// tri[0] is the point of the triang on edge1
    // tri[0][0] is the actual x dimension, tri[0][1] is the actual y-dimension, tri[0][2] is the actual z dimension
// tri[1] is the point of the triang on edge2
// tri[2] is the point of the triang on edge3
// int** vertices is the three edges described by vertex endpoints

float** InterpolateTriangle(float** tri, int** tri_by_edges, int cellId, float isoval, const int* dims, const float* F, const float* X, const float* Y, const float* Z){
                InterpolateTriangle(actual_triangle, new_triangle_by_edges, i, IsoVal, dims, F, X, Y, Z);


    int pt0[2] = Edges[tri_by_edges[0]]; // edge 0 between vertex 0 and vertex 1
    int pt1[2] = Edges[tri_by_edges[1]] // edge 3 between vertex 0 and vertex 2
    int pt2[2] = Edges[tri_by_edges[2]] // edge 8 between vertex 0 and vertex 4
    int all_vertices[8];
    GetEightVertices(all_vertices, cellId, dims); //now vertices is the point indices
    // contains ptind of first vertex, second vertext

    // logical pt indices of six points making up three point on edges
    int edge1pt1_ind = all_vertices[pt0[0]];
    int edge1pt2_ind = all_vertices[pt0[1]];
    int edge2pt1_ind = all_vertices[pt1[0]];
    int edge2pt2_ind = all_vertices[pt1[1]];
    int edge3pt1_ind = all_vertices[pt2[0]];
    int edge3pt2_ind = all_vertices[pt2[1]];

    int edge1pt1_logical[3];
    int edge1pt2_logical[3];
    int edge2pt1_logical[3];
    int edge2pt2_logical[3];
    int edge3pt1_logical[3];
    int edge3pt2_logical[3];
    GetLogicalPointIndex(edge1pt1_logical, edge1pt1_ind, dims);
    GetLogicalPointIndex(edge1pt2_logical, edge1pt2_ind, dims);
    GetLogicalPointIndex(edge2pt1_logical, edge2pt1_ind, dims);
    GetLogicalPointIndex(edge2pt2_logical, edge2pt2_ind, dims);
    GetLogicalPointIndex(edge3pt1_logical, edge3pt1_ind, dims);
    GetLogicalPointIndex(edge3pt2_logical, edge3pt2_ind, dims);

    // now we have point inds and logical point inds for all six vertices making up the three edges
    float CalculatePosition(float FX, float FA, float FB, float A, float B){

    float first = CalculatePosition(isoval, F[edge1pt1_ind], F[edge1pt2], X);
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("test_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    

    TriangleList tl;
    // do I need to add borders around the frame like in proj 5?

    //static int Edges[12][2] = 
    //{{0, 1}, {1, 3}, {2, 3}, {0, 2},
    // {4, 5}, {5, 7}, {6, 7}, {4, 6},
    // {0, 4}, {1, 5}, {2, 6}, {3, 7} };

    float IsoVal = 3.2;
    int icase;
    int edge0 = 0;
    int edge1 = 0;
    int edge2 = 0;
    int new_triangle_by_edges[3] = {0, 0, 0};
    int triangle_by_vertex[3][2] = {{0,0},{0,0},{0,0}};
    float actual_triangle[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    float pt1[2] = {0.0, 0.0};
    float pt2[2] = {0.0, 0.0};
    //iterate through cells
    fprintf(stdout, "num cells is %d\n", GetNumberOfCells(dims));
    for (i = 0; i < GetNumberOfCells(dims); i++){
        icase = IdentifyCase(i, IsoVal, dims, F);
        int* edges = TriCases[icase];
        while (*edges != -1){
            edge0 = *edges++;
            edge1 = *edges++;
            edge2 = *edges++;
            new_triangle_by_edges = {edge0, edge1, edge2}; //0, 3, 8


            InterpolateTriangle(actual_triangle, new_triangle_by_edges, i, IsoVal, dims, F, X, Y, Z);
            //get real points now
            // how do we want it?
            // actual xyz coords of each of three points
            AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)

        }






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
            fprintf(stdout, "\tll is %f,%f has val %f\n", X[ll[0]], Y[ll[1]], Fll);
            fprintf(stdout, "\tlr is %f,%f has val %f\n", X[lr[0]], Y[lr[1]], Flr);
            fprintf(stdout, "\tul is %f,%f has val %f\n", X[ul[0]], Y[ul[1]], Ful);
            fprintf(stdout, "\tur is %f,%f has val %f\n", X[ur[0]], Y[ur[1]], Fur);
            if (edge1 == 0){
                //interpolate ll -> lr
                pt1[0] = CalculatePosition(IsoVal, Fll, Flr, X[ll[0]], X[lr[0]]);
                pt1[1] = Y[ll[1]];
            }
            else if (edge1 == 1){
                //interpolate lr -> ur
                pt1[0] = X[lr[0]];
                pt1[1] = CalculatePosition(IsoVal, Flr, Fur, Y[lr[1]], Y[ur[1]]); // ****** comment this!
            }
            else if (edge1 == 2){
                //interpolate ul -> ur
                pt1[0] = CalculatePosition(IsoVal, Ful, Fur, X[ul[0]], X[ur[0]]);
                pt1[1] = Y[ul[1]];
            }
            else if (edge1 == 3){
                //interpolate ll -> ul
                pt1[0] = X[ll[0]];
                pt1[1] = CalculatePosition(IsoVal, F[llInd], F[ulInd], Y[ll[1]], Y[ul[1]]); //*** comment this!
            }
            if (edge2 == 0){
                //interpolate ll -> lr
                pt2[0] = CalculatePosition(IsoVal, F[llInd], F[lrInd], X[ll[0]], X[lr[0]]);
                pt2[1] = Y[ll[1]];
            }
            else if (edge2 == 1){
                //interpolate lr -> ur
                pt2[0] = X[lr[0]];
                pt2[1] = CalculatePosition(IsoVal, F[lrInd], F[urInd], Y[lr[1]], Y[ur[1]]); //*comment this!!!
            }
            else if (edge2 == 2){
                //interpolate ul -> ur
                pt2[0] = CalculatePosition(IsoVal, F[ulInd], F[urInd], X[ul[0]], X[ur[0]]);
                pt2[1] = Y[ul[1]];
            }
            else if (edge2 == 3){
                //interpolate ll -> ul
                pt2[0] = X[ll[0]];
                pt2[1] = CalculatePosition(IsoVal, F[llInd], F[ulInd], Y[ll[1]], Y[ul[1]]); ///**** comment this!
            }
            fprintf(stdout, "\tpt1x %f,%f, pt2x %f,%f\n", pt1[0], pt1[1], pt2[0], pt2[1]);
            //sl.AddSegment(pt1[0], pt1[1], pt2[0], pt2[1]);
        }
    }

    vtkPolyData *pd = tl.MakePolyData();

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
