#include "FitsReader.h"
#include <ctime>

std::ofstream output_file;
#define LOG(msg)
//#define LOG(msg) vtkWarningMacro(<< __PRETTY_FUNCTION__ << ", l" << __LINE__ << ": " << msg);
//#define LOG(msg) {output_file << __PRETTY_FUNCTION__ << ", l" << __LINE__ << ": " << msg << std::endl;}
template <class T>
void suppress_unused_warning(const T &) {}

#include <time.h>
#include <assert.h>

#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStringArray.h"

#include <casa/Utilities/CountedPtr.h>
#include <images/Images/FITSImage.h>
#include <images/Images/ImageOpener.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>

#include <mpi.h>

vtkStandardNewMacro(FitsReader)

#ifdef _MSC_VER
// Let us get rid of this funny warning on /W4:
// warning C4611: interaction between '_setjmp' and C++ object
// destruction is non-portable
#pragma warning( disable : 4611 )
#endif

FitsReader::FitsReader()
{
    ImageOpener::registerOpenImageFunction(ImageOpener::FITS, &FITSImage::openFITSImage);
}
FitsReader::~FitsReader()
{
}

//----------------------------------------------------------------------------
void FitsReader::ExecuteInformation()
{
    int rank = 0;
    //MPI_Comm_rank (MPI_COMM_WORLD, &rank);    /* get current process id */


    if (!output_file.is_open()) {
        char filename[80];
        sprintf(filename, "/tmp/%d.txt", rank);
        output_file.open(filename);
        LOG("thisptr: " << this);
    }

    using namespace casa;

    this->ComputeInternalFileName(this->DataExtent[4]);
    if (this->InternalFileName == NULL)
    {
        return;
    }

    // Open the image (of any type). Make sure the data type is float.
    // Note that use of CountedPtr takes care of automatic object deletion.
    CountedPtr<LatticeBase> casa_lattice(
                    ImageOpener::openImage(this->InternalFileName));
    ImageInterface<float>* casa_image =
                    dynamic_cast<ImageInterface<float>*>(casa_lattice.operator->());
    if (casa_image == NULL)
    {
        vtkErrorMacro(
                        << "2. Couldn't convert the fits file to floats: " << this->InternalFileName << ".");
        return;
    }

    this->SetNumberOfScalarComponents(1);
    this->SetDataScalarTypeToFloat();
    this->DataExtent[0] = 0;
    this->DataExtent[1] = casa_image->shape()(0) - 1;
    this->DataExtent[2] = 0;
    this->DataExtent[3] = casa_image->shape()(1) - 1;
    this->DataExtent[4] = 0;
    this->DataExtent[5] = 0;
    if (casa_image->shape().size() >= 3)
    {
        this->DataExtent[4] = 0;
        this->DataExtent[5] = casa_image->shape()(2) - 1;
    }
    LOG("Dimensions:      "
                    << this->DataExtent[1] << ", "
                    << this->DataExtent[3] << ", "
                    << this->DataExtent[5]);

    Vector<double> increment = casa_image->coordinates().increment();
    Vector<Double> reference_value = casa_image->coordinates().referenceValue();
    Vector<Double> reference_pixel = casa_image->coordinates().referencePixel();

    increment[2] *= 1.0e-9;
    reference_value[2] = 0.0;

    // The reference pixel in Paraview is always (0,0,0).
    for (int i=0; i<3; ++i)
        reference_value[i] -= reference_pixel[i] * increment[i];

    for (int i=0; i<3; ++i)
        increment[i] = std::abs(increment[i]);

    this->DataSpacing[0] = increment[0];
    this->DataSpacing[1] = increment[1];
    this->DataSpacing[2] = increment[2];

    this->SetDataOrigin(reference_value[0], reference_value[1], reference_value[2]);

    LOG("Spacing:         "
                    << this->DataSpacing[0] << ", "
                    << this->DataSpacing[1] << ", "
                    << this->DataSpacing[2]);
    LOG("Reference value: "
                    << reference_value[0] << ", "
                    << reference_value[1] << ", "
                    << reference_value[2]);
    LOG("Reference pixel: "
                    << reference_pixel[0] << ", "
                    << reference_pixel[1] << ", "
                    << reference_pixel[2]);
}

//----------------------------------------------------------------------------
// This function reads a data from a file.  The datas extent/axes
// are assumed to be the same as the file extent/order.
#if PARAVIEW_VERSION_MINOR <= 20
void FitsReader::ExecuteData(vtkDataObject* output)
#else
void FitsReader::ExecuteDataWithInformation(vtkDataObject* output,
                vtkInformation* outInfo)
#endif
{
    using namespace casa;

#if PARAVIEW_VERSION_MINOR <= 20
    vtkImageData* data = this->AllocateOutputData(output);
#else
    vtkImageData* data = this->AllocateOutputData(output, outInfo);
#endif

    if (this->InternalFileName == NULL)
    {
        vtkErrorMacro(<< "Either a FileName or FilePrefix must be specified.");
        return;
    }

    // Open the image (of any type). Make sure the data type is float.
    // Note that use of CountedPtr takes care of automatic object deletion.
    CountedPtr<LatticeBase> casa_lattice(
                    ImageOpener::openImage(this->InternalFileName));
    ImageInterface<float>* casa_image =
                    dynamic_cast<ImageInterface<float>*>(casa_lattice.operator->());

    if (casa_image == NULL)
    {
        vtkErrorMacro(
                        << "3. Couldn't load data: " << this->InternalFileName << ".");
        return;
    }

    data->GetPointData()->GetScalars()->SetName("intensity");

    this->ComputeDataIncrements();

    IPosition shape = casa_image->shape();
    IPosition pos(shape.size(), 0);

    clock_t start_copy_time = std::clock();
    {
        Vector<double> increment = casa_image->coordinates().increment();

        bool flip[4];
        for (size_t i=0; i<increment.size(); ++i) flip[i] = (increment[i]<0);

        int ext[6];
        data->GetExtent(ext);
        int dx,dy,dz, x0,y0,z0;
        dx = ext[1]-ext[0]; // Careful one off: -1
        dy = ext[3]-ext[2]; // Careful one off: -1
        dz = ext[5]-ext[4]; // Careful one off: -1
        x0 = ext[0];
        y0 = ext[2];
        z0 = ext[4];
        if (flip[0]) x0 = this->DataExtent[1] - ext[1];
        if (flip[1]) y0 = this->DataExtent[3] - ext[3];
        if (flip[2]) z0 = this->DataExtent[5] - ext[5];
        LOG(std::string("Extent: ") << ext[0] << " " << ext[1] << " " << ext[2] << " " << ext[3] << " " << ext[4] << " " << ext[5]);
//        LOG("flip: " << flip[0] << " " << flip[1] << " " << flip[2]);

        switch (increment.size()) {
        case 2: {
            assert(false && "Not yet implemented, no 2D data yet.");
            break;
        }
        case 3: {
            Slicer s(IPosition(3, x0,y0,z0), IPosition(3, dx,dy,dz));
            Cube<float> cube(casa_image->get(false)(s));
            float *output_data = (float*)data->GetScalarPointer();

            for (int k=0; k<=dz; ++k)
            {
                int z = (flip[2]?dz-k:k);
                for (int j=0; j<=dy; ++j)
                {
                    int y = (flip[1]?dy-j:j);
                    for (int i=0; i<=dx; ++i)
                    {
                        int x = (flip[0]?dx-i:i);
                        *output_data = cube(x, y, z);
                        ++output_data;
                    }

                }
            }
            break;
        }
        case 4: {
            Slicer s(IPosition(4, x0,y0,z0,0), IPosition(4, dx,dy,dz,0));
            Array<float> cube(casa_image->get(false)(s));
            float *output_data = (float*)data->GetScalarPointer();

            for (int k=0; k<=dz; ++k)
            {
                int z = (flip[2]?dz-k:k);
                for (int j=0; j<=dy; ++j)
                {
                    int y = (flip[1]?dy-j:j);
                    for (int i=0; i<=dx; ++i)
                    {
                        int x = (flip[0]?dx-i:i);
                        *output_data = cube(IPosition(4, x, y, z, 0));
                        ++output_data;
                    }

                }
            }
            break;
        }
        }
    }
    clock_t end_copy_time = std::clock();
    double dt = ((end_copy_time-start_copy_time)*1.0/CLOCKS_PER_SEC);
    suppress_unused_warning(dt);


    {
        vtkFieldData* fd = vtkFieldData::New();

        const DirectionCoordinate& coordinate =
                        casa_image->coordinates().directionCoordinate(0);
        const char* names[] = { "AxisTitleForX", "AxisTitleForY", "AxisTitleForZ" };
        vtkStringArray* TitleArray;
        for (size_t i = 0; i < coordinate.axisNames(MDirection::DEFAULT).size() && i < 3; ++i)
        {
            TitleArray = vtkStringArray::New();
            TitleArray->SetName(names[i]);
            TitleArray->InsertNextValue((
                            coordinate.axisNames(MDirection::DEFAULT)[i] + " (" +
                            MDirection(coordinate.directionType()).getRefString() + ")").c_str());

            fd->AddArray(TitleArray);
            TitleArray->Delete();
        }

        output->SetFieldData(fd);
    }
}

//----------------------------------------------------------------------------
int FitsReader::CanReadFile(const char* fname)
{
    ImageOpener::ImageTypes image_type = ImageOpener::imageType(fname);
    if (image_type == ImageOpener::FITS || image_type == ImageOpener::HDF5)
        return 3;

    return 0;
}
#ifdef _MSC_VER
// Put the warning back
#pragma warning( default : 4611 )
#endif

//----------------------------------------------------------------------------
void FitsReader::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}
