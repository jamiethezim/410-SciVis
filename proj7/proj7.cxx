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
#include "tricase.cxx"

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
    idx[0] = pointId%dims[0];
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
    int v0[3] = {coords[0], coords[1], coords[2]};
    int v1[3] = {coords[0]+1, coords[1], coords[2]};
    int v2[3] = {coords[0], coords[1]+1, coords[2]};
    int v3[3] = {coords[0]+1, coords[1]+1, coords[2]};
    int v4[3] = {coords[0], coords[1], coords[2]+1};
    int v5[3] = {coords[0]+1, coords[1], coords[2]+1};
    int v6[3] = {coords[0], coords[1]+1, coords[2]+1};
    int v7[3] = {coords[0]+1, coords[1]+1, coords[2]+1};

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

// *************************************
// given an edge case number (which edge?), and the ptindices of the endpoints of that egge
// determine what A and B should be, and if they change in X direction, Y direction, and Z direction
// interpolate the mid-distance location
// and store the known other coordinates away too
// return value is float** location, such that locations[0] is the X value, locations[1] is the Y value, and locations[2] Zvalue
void DetermineLocation(float* locations, int edgenum, int pt1, int pt2, float isoval, const float* F, const float* X, const float* Y, const float* Z, const int* dims){
    // need to get logical point indices of the two endpoints so we can get actual XYZ locations
    int pt1_log[3];
    int pt2_log[3];
    GetLogicalPointIndex(pt1_log, pt1, dims);
    GetLogicalPointIndex(pt2_log, pt2, dims);

    float FA, FB, A, B, Xval, Yval, Zval = 0.0;

    if ((edgenum == 0) || (edgenum == 2) || (edgenum == 4) || (edgenum == 6)){
        // we know Y and Z, but don't know X
        A = X[pt1_log[0]];
        B = X[pt2_log[0]];
        Xval = CalculatePosition(isoval, F[pt1], F[pt2], A, B);
        Yval = Y[pt1_log[1]];
        Zval = Z[pt1_log[2]];
    } else if ((edgenum == 1) || (edgenum == 3) || (edgenum == 5) || (edgenum == 7)){
        // we know X and Z, but don't know Y
        Xval = X[pt1_log[0]];
        A = Y[pt1_log[1]];
        B = Y[pt2_log[1]];
        Yval = CalculatePosition(isoval, F[pt1], F[pt2], A, B);
        Zval = Z[pt1_log[2]];
    } else { //edgenum 8, 9, 10, 11 are all delta Z
        // we know X and Y but don't know Z
        Xval = X[pt1_log[0]];
        Yval = Y[pt1_log[1]];
        A = Z[pt1_log[2]];
        B = Z[pt2_log[2]];
        Zval = CalculatePosition(isoval, F[pt1], F[pt2], A, B);
    }
    // now store the results in the return value pointer, and locations contains the actual location of the vertex of the triangle
    locations[0] = Xval;
    locations[1] = Yval;
    locations[2] = Zval;
}


// ******************************************************************************************
// helper function
// tri is an input to the function that will contain the return value
// returns 3 vertices composed of 3 endpoints
// tri[0] is the point of the triang on edge1
    // tri[0][0] is the actual x dimension, tri[0][1] is the actual y-dimension, tri[0][2] is the actual z dimension
// tri[1] is the point of the triang on edge2
// tri[2] is the point of the triang on edge3
// int** tri_by_edges is the three edges described by vertex endpoints
// static int Edges[12][2] = 
//    {{0, 1}, {1, 3}, {2, 3}, {0, 2},
//     {4, 5}, {5, 7}, {6, 7}, {4, 6},
//     {0, 4}, {1, 5}, {2, 6}, {3, 7} };

void InterpolateTriangle(float** tri, int* tri_by_edges, int cellId, float isoval, const int* dims, const float* F, const float* X, const float* Y, const float* Z){
    int pt0[2] = {Edges[tri_by_edges[0]][0], Edges[tri_by_edges[0]][1]}; // edge 0 between vertex 0 and vertex 1
    int pt1[2] = {Edges[tri_by_edges[1]][0], Edges[tri_by_edges[1]][1]}; // edge 3 between vertex 0 and vertex 2
    int pt2[2] = {Edges[tri_by_edges[2]][0], Edges[tri_by_edges[2]][1]}; // edge 8 between vertex 0 and vertex 4

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

    float XYZlocation[3] = {0.0, 0.0, 0.0};
    float x1, y1, z1, x2, y2, z2, x3, y3, z3 = 0.0;

    DetermineLocation(XYZlocation, tri_by_edges[0], edge1pt1_ind, edge1pt2_ind, isoval, F, X, Y, Z, dims);
    x1 = XYZlocation[0];
    y1 = XYZlocation[1];
    z1 = XYZlocation[2];

    DetermineLocation(XYZlocation, tri_by_edges[1], edge2pt1_ind, edge2pt2_ind, isoval, F, X, Y, Z, dims);
    x2 = XYZlocation[0];
    y2 = XYZlocation[1];
    z2 = XYZlocation[2];

    DetermineLocation(XYZlocation, tri_by_edges[2], edge3pt1_ind, edge3pt2_ind, isoval, F, X, Y, Z, dims);
    x3 = XYZlocation[0];
    y3 = XYZlocation[1];
    z3 = XYZlocation[2];

    // set the return value!
    tri[0][0] = x1;
    tri[0][1] = y1;
    tri[0][2] = z1;
    tri[1][0] = x2;
    tri[1][1] = y2;
    tri[1][2] = z2;
    tri[2][0] = x3;
    tri[2][1] = y3;
    tri[2][2] = z3;
}

int main()
{
    int i;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("test_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
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
    float actual_triangle[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    float x1, y1, z1, x2, y2, z2, x3, y3, z3 = 0.0;
    //iterate through cells
    fprintf(stdout, "num cells is %d\n", GetNumberOfCells(dims));
    for (i = 0; i < GetNumberOfCells(dims); i++){
        icase = IdentifyCase(i, IsoVal, dims, F);
        int* edges = triCase[icase];
        while (*edges != -1){
            edge0 = *edges++;
            edge1 = *edges++;
            edge2 = *edges++;
            new_triangle_by_edges[0] = edge0;
            new_triangle_by_edges[1] = edge1;
            new_triangle_by_edges[2] = edge2; //0, 3, 8

            //get real points now
            InterpolateTriangle(actual_triangle, new_triangle_by_edges, i, IsoVal, dims, F, X, Y, Z);
            x1 = actual_triangle[0][0];
            y1 = actual_triangle[0][1];
            z1 = actual_triangle[0][2];
            x2 = actual_triangle[1][0];
            y2 = actual_triangle[1][1];
            z2 = actual_triangle[1][2];
            x3 = actual_triangle[2][0];
            y3 = actual_triangle[2][1];
            z3 = actual_triangle[2][2];

            tl.AddTriangle(x1, y1, z1, x2, y2, z2, x3, y3, z3);
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
