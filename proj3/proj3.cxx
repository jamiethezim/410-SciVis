#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
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

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F)
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


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void
ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    RGB[0] = (unsigned char) (F*(255-0));
    RGB[1] = (unsigned char) (F*(255-0));
    RGB[2] = (unsigned char) (F*(255-128)+128);
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128) 
//        F=0.5: (255,255,255) 
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    unsigned char r, g, b;
    if (F == 0){
        r = 0;
        g = 0;
        b = 128;
    } else if ((F>0) && (F<0.5)){
        r = (unsigned char) ((F*2)*(255-0));
        g = (unsigned char) ((F*2)*(255-0));
        b = (unsigned char) ((F*2)*(255-128)+128);
    } else if (F==0.5){
        r = 255;
        g = 255;
        b = 255;
    } else if ((F>0.5) && (F<1)){
        r = (unsigned char) ((1-((F-0.5)/(1-0.5))) * (255-128) + 128);
        g = (unsigned char) ((1-((F-0.5)/(1-0.5))) *(255-0));
        b = (unsigned char) ((1-((F-0.5)/(1-0.5))) *(255-0));
    } else if (F==1){
        r = 128;
        g = 0;
        b = 0;
    }
    RGB[0] = r;
    RGB[1] = g;
    RGB[2] = b;
}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char *RGB)
{
    float r, g, b;
    float saturation = 1.0;
    float value = 1.0;
    float hue = (F*(360-0))/60.0;
    int i = floor(hue);
    float newf = hue - i;
    float p = value * (1.0 - saturation);
    float q = value * (1.0 - saturation*newf);
    float t = value * (1.0 - saturation*(1.0-newf));

    switch(i){
        case 0:
            r = value;
            g = t;
            b = p;
            break;
        case 1:
            r = q;
            g = value;
            b = p;
            break;
        case 2:
            r = p;
            g = value;
            b = t;
            break;
        case 3:
            r = p;
            g = q;
            b = value;
            break;
        case 4:
            r = t;
            g = p;
            b = value;
            break;
        case 5:
            r = value;
            g = p;
            b = q;
            break;
    }

    RGB[0] = r*255;
    RGB[1] = g*255;
    RGB[2] = b*255;
}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    /*fprintf(stdout, "dims are %d x %d\n", dims[0], dims[1]);
    for (int k =0; k<10; k++){
        fprintf(stdout, "coords %f x %f : %f\n", X[k], Y[k], F[k]);
    }*/
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3*nx*ny ; i++)
        for (j = 0 ; j < 3 ; j++)
            buffer[j][i] = 0;

    for (i = 0 ; i < nx ; i++)
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            pt[0] = (i/(nx-1.0)) * 18 + -9;
            pt[1] = (j/(ny-1.0)) * 18 + -9;
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
            float normalizedF = (f-1.2)/(5.02-1.2);
            
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
