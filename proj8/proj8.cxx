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

#include <vtkLookupTable.h> // for the color look up table - don't need this
#include <vtkPlane.h> // for slicing
#include <vtkCutter.h> //for slicing
#include <vtkHedgeHog.h>
#include <vtkGlyph3D.h> // use either glyph or hedgehog?
#include <vtkArrowSource.h> // for glyphing
#include <vtkMaskPoints.h>
#include <vtkPointSource.h> // for streamlines
#include <vtkLineSource.h>
#include <vtkRungeKutta4.h> // for particle advection
#include <vtkStreamTracer.h> //for streamlines

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


int main()
{

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj8.vtk");
    rdr->Update();
    rdr->GetOutput()->GetPointData()->SetActiveAttribute("grad", vtkDataSetAttributes::VECTORS);

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);
    //rgrid->SetSampleRate(10,10,10);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    //float *F_vectors = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);

    // (xmin, ymin, xmax, ymax)
    double bottomleft[4] = {0.0, 0.0, 0.5, 0.5};
    double topleft[4] = {0.0, 0.5, 0.5, 1.0};
    double bottomright[4] = {0.5, 0.0, 1.0, 0.5};
    double topright[4] = {0.5, 0.5, 1.0, 1.0};

    // create the renderers, window, interactor
    vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren2 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren3 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren4 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    
    //and some source/map/actor stuff
    vtkSmartPointer<vtkSphereSource> sphere_source = vtkSmartPointer<vtkSphereSource>::New();
    vtkSmartPointer<vtkPolyDataMapper> sphere_map = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> sphere_actor1 = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> sphere_actor2 = vtkSmartPointer<vtkActor>::New();

    // but now the actual actors we care about
    vtkSmartPointer<vtkPolyDataMapper> contour_map = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> contour_actor = vtkSmartPointer<vtkActor>::New();
    vtkContourFilter *cf = vtkContourFilter::New();
    cf->SetNumberOfContours(2);
    //Set the ith contour value.
    cf->SetValue(0, 2.5);
    cf->SetValue(1, 5.0);
    cf->SetInputConnection(rdr->GetOutputPort());

    // create color map
    //vtkSmartPointer<vtkLookupTable> colortable = vtkSmartPointer<vtkLookupTable>::New();
    //colortable->SetNumberOfColors(2);
    //colortable->SetHueRange(0.0, 0.8);
    //colortable->SetTableRange(0,2);
    //colortable->Build();
    // map isosurface to actor
    contour_map->SetInputConnection(cf->GetOutputPort());
    //contour_map->SetLookupTable(colortable);
    contour_map->SetScalarRange(0.0, 6.0);
    contour_actor->SetMapper(contour_map);
    //contour_actor->GetProperty()->SetColor(1,1,0);



    // render 2 - slices
    // get the vector data from rdr
    //vtkDataSet *scalarData = rdr->GetOutput();
    //scalarData->GetPointData()->SetActiveScalars("hardyglobal");
    //scalarData->Update();

    vtkSmartPointer<vtkPolyDataMapper> cubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cubeMapper->SetInputConnection(rdr->GetOutputPort());
    cubeMapper->SetScalarRange(2.5, 5.0);


    // Create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
    vtkSmartPointer<vtkPlane> plane1 = vtkSmartPointer<vtkPlane>::New();
    plane1->SetOrigin(rdr->GetOutput()->GetCenter());
    plane1->SetNormal(0,0,1);

    // Create a second plane
    vtkSmartPointer<vtkPlane> plane2 = vtkSmartPointer<vtkPlane>::New();
    plane2->SetOrigin(rdr->GetOutput()->GetCenter());
    plane2->SetNormal(0,1,0);

        // Create a second plane
    vtkSmartPointer<vtkPlane> plane3 = vtkSmartPointer<vtkPlane>::New();
    plane3->SetOrigin(rdr->GetOutput()->GetCenter());
    plane3->SetNormal(1,0,0);
 
    // Create cutter
    vtkSmartPointer<vtkCutter> cutter1 = vtkSmartPointer<vtkCutter>::New();
    cutter1->SetCutFunction(plane1);
    cutter1->SetInputConnection(rdr->GetOutputPort());
    cutter1->Update();

    // Create cutter
    vtkSmartPointer<vtkCutter> cutter2 = vtkSmartPointer<vtkCutter>::New();
    cutter2->SetCutFunction(plane2);
    cutter2->SetInputConnection(rdr->GetOutputPort());
    cutter2->Update();

     // Create cutter
    vtkSmartPointer<vtkCutter> cutter3 = vtkSmartPointer<vtkCutter>::New();
    cutter3->SetCutFunction(plane3);
    cutter3->SetInputConnection(rdr->GetOutputPort());
    cutter3->Update();
 
    vtkSmartPointer<vtkPolyDataMapper> cutterMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    cutterMapper1->SetInputConnection( cutter1->GetOutputPort());
    cutterMapper1->SetScalarRange(1.5, 5.0);
    vtkSmartPointer<vtkPolyDataMapper> cutterMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    cutterMapper2->SetInputConnection( cutter2->GetOutputPort());
    cutterMapper2->SetScalarRange(1.5, 5.0);
    vtkSmartPointer<vtkPolyDataMapper> cutterMapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
    cutterMapper3->SetInputConnection( cutter3->GetOutputPort());
    cutterMapper3->SetScalarRange(1.5, 5.0);

 
   // Create cube actor
    vtkSmartPointer<vtkActor> cubeActor = vtkSmartPointer<vtkActor>::New();
    //cubeActor->GetProperty()->SetColor(0.5,1,0.5);
    cubeActor->GetProperty()->SetOpacity(0.1);
    cubeActor->SetMapper(cubeMapper);
  // Create plane actor
    vtkSmartPointer<vtkActor> plane1Actor = vtkSmartPointer<vtkActor>::New();
    //plane1Actor->GetProperty()->SetColor(1.0,1,0);
    plane1Actor->GetProperty()->SetLineWidth(3);
    plane1Actor->SetMapper(cutterMapper1);

  // Create plane actor
    vtkSmartPointer<vtkActor> plane2Actor = vtkSmartPointer<vtkActor>::New();
    //plane2Actor->GetProperty()->SetColor(1.0,1,0);
    plane2Actor->GetProperty()->SetLineWidth(3);
    plane2Actor->SetMapper(cutterMapper2);

  // Create plane actor
    vtkSmartPointer<vtkActor> plane3Actor = vtkSmartPointer<vtkActor>::New();
    //plane3Actor->GetProperty()->SetColor(1.0,1,0);
    plane3Actor->GetProperty()->SetLineWidth(3);
    plane3Actor->SetMapper(cutterMapper3);


    // Render 3 - hedgehog
   
    //vtkSmartPointer<vtkHedgeHog> hedgehog = vtkSmartPointer<vtkHedgeHog>::New();
    //hedgehog->SetInputConnection(rdr->GetOutputPort()); // or is it SetInputData(sgrid)? 
    //hedgehog->SetScaleFactor(0.8);
    //vtkSmartPointer<vtkPolyDataMapper> hedgeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    //hedgeMapper->SetInputConnection(hedgehog->GetOutputPort());
    //hedgeMapper->SetScalarRange(2.5, 5.0);
    //vtkSmartPointer<vtkActor> HedgeActor = vtkSmartPointer<vtkActor>::New();
    //HedgeActor->SetMapper(hedgeMapper);


    // vtkArrowSource arrow
    // arrow SetTipResolution 6
    // arrow SetTipRadius 0.1
    // arrow SetTipLength 0.35
    // arrow SetShaftResolution 6
    // arrow SetShaftRadius 0.03
    // vtkGlyph3D glyph
    //glyph SetInput [reader GetOutputPort]
    //glyph SetSource [arrow GetOutputPort]
    //glyph SetVectorModeToUseVector
    //glyph SetColorModeToColorByScalar
    //glyph SetScaleModeToDataScalingOff
    //glyph OrientOn
    //glyph SetScaleFactor 0.2

    vtkSmartPointer<vtkArrowSource> arrow = vtkSmartPointer<vtkArrowSource>::New();
    arrow->SetTipResolution(6);
    arrow->SetTipRadius(0.1);
    arrow->SetTipLength(0.35);
    arrow->SetShaftResolution(6);
    arrow->SetShaftRadius(0.03);

    vtkSmartPointer<vtkMaskPoints> ptMask = vtkSmartPointer<vtkMaskPoints>::New();
    ptMask->SetInputConnection(rdr->GetOutputPort());
    ptMask->SetOnRatio(63);
    ptMask->Update();

    vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
    glyph->SetSourceConnection(arrow->GetOutputPort());
    glyph->SetInputData(ptMask->GetOutput()); // or is it set input (rdr get Output port)?
    glyph->SetVectorModeToUseVector();
    glyph->SetColorModeToColorByScalar();
    glyph->SetScaleModeToDataScalingOff();
    glyph->OrientOn();
    glyph->SetScaleFactor(1.5);
    glyph->Update();
    vtkSmartPointer<vtkPolyDataMapper> glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(glyph->GetOutputPort());
    glyphMapper->SetScalarRange(2.5, 5.0); //for now?
    vtkSmartPointer<vtkActor> glyphActor = vtkSmartPointer<vtkActor>::New();
    glyphActor->SetMapper(glyphMapper);

    //HedgeActor->GetProperty()->SetColor(0,0,0);



    // renderer 4 - streamlines
    // vtkPointSource seeds
    //   seeds SetRadius 3.0
    //   eval seeds SetCenter [[reader GetOutput] GetCenter]
    //   seeds SetNumberOfPoints 100    //

    //vtkRungeKutta4 integ    //

    //vtkStreamTracer streamer
    //   streamer SetInputConnection  [reader GetOutputPort]
    //   streamer SetSourceConnection [seeds GetOutputPort]
    //   streamer SetMaximumPropagation 100
    //   streamer SetMaximumPropagationUnitToTimeUnit
    //   streamer SetInitialIntegrationStepUnitToCellLengthUnit
    //   streamer SetInitialIntegrationStep 0.1
    //   streamer SetIntegrationDirectionToBoth
    //   streamer SetIntegrator integ    //

    //vtkPolyDataMapper mapStreamLines
    //   mapStreamLines SetInputConnection [streamer GetOutputPort]
    //   eval mapStreamLines SetScalarRange [[reader GetOutput] GetScalarRange]

    vtkSmartPointer<vtkLineSource> seeds = vtkSmartPointer<vtkLineSource>::New();
    seeds->SetPoint1(-9,0,0);
    seeds->SetPoint2(9,0,0);
    seeds->SetResolution(19);

    // create rk4 for particle advection
    vtkSmartPointer<vtkStreamTracer> streamer = vtkSmartPointer<vtkStreamTracer>::New();
    streamer->SetInputConnection(rdr->GetOutputPort());
    streamer->SetSourceConnection(seeds->GetOutputPort());
    streamer->SetMaximumPropagation(100);
    //streamer->SetMaximumPropagationUnitToTimeUnit();
    //streamer->SetInitialIntegrationStepUnitToCellLengthUnit();
    streamer->SetInitialIntegrationStep(0.1);
    streamer->SetIntegratorTypeToRungeKutta4();

    vtkSmartPointer<vtkPolyDataMapper> streamlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    streamlineMapper->SetInputConnection(streamer->GetOutputPort());
    streamlineMapper->SetScalarRange(2.9, 5.0); //rdr->GetOutput()->GetScalarRange()

    vtkSmartPointer<vtkActor> streamActor = vtkSmartPointer<vtkActor>::New();
    streamActor->SetMapper(streamlineMapper);
    
    // add renderer to window, and adjust placing
    renWin->AddRenderer(ren1);
    ren1->SetViewport(bottomleft);
    renWin->AddRenderer(ren2);
    ren2->SetViewport(topleft);
    renWin->AddRenderer(ren3);
    ren3->SetViewport(bottomright);
    renWin->AddRenderer(ren4);
    ren4->SetViewport(topright);

    // put the isosurface contour filter in ll corner
    ren1->AddActor(contour_actor);
    // put the cuts in render 2
    ren2->AddActor(plane1Actor);
    ren2->AddActor(plane2Actor);
    ren2->AddActor(plane3Actor);
    ren2->AddActor(cubeActor);
    ren3->AddActor(glyphActor);
    ren4->AddActor(streamActor);

    // start the window
    renWin->SetSize(800,800);
    iren->SetRenderWindow(renWin);
    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */
/*
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
*/
    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

}
