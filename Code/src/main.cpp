#include <casa/Utilities/CountedPtr.h>
#include <images/Images/FITSImage.h>
#include <images/Images/ImageOpener.h>
#include <images/Images/ImageProxy.h>

int main(int argc, char* argv[])
{
	using namespace casa;

	ImageOpener::registerOpenImageFunction(ImageOpener::FITS, &FITSImage::openFITSImage);

	const char* filename = "~/Projects/LOFAR/Data/h5img/L70181_SBG001_sky.h5";
	if (argc > 1)
		filename = argv[1];


	// Open the image (of any type). Make sure the data type is float.
	// Note that use of CountedPtr takes care of automatic object deletion.
	CountedPtr<LatticeBase> lattice (ImageOpener::openImage (filename));

	ImageInterface<float>* imagePtr = dynamic_cast<ImageInterface<float>*>(lattice.operator->());
	if (imagePtr == NULL)
	{
		std::cout << "2. Couldn't convert the fits file to floats: " << filename << "." << std::endl;
		return 1;
	}

	IPosition shape = imagePtr->shape();
	if (true)
	{
		std::cout << "Dimensions: " << std::endl;
		for (size_t i=0; i<shape.size(); ++i)
			std::cout << " [" << i << "]: " << shape(i) << std::endl;

		std::cout << "Increment: " << std::endl;
		for (size_t i=0; i<shape.size(); ++i)
			std::cout << " [" << i << "]: " << imagePtr->coordinates().increment()[i] << std::endl;

		CoordinateSystem coord_sys = imagePtr->coordinates();
		Vector<Double> reference_value = coord_sys.referenceValue();
		std::cout << "Reference pixel: " << reference_value[0] << ", "
				<< reference_value[1] << ", "
				<< reference_value[2] << std::endl;

		reference_value[2] = 0;
		coord_sys.setReferenceValue(reference_value);
		imagePtr->setCoordinateInfo(coord_sys);

		coord_sys = imagePtr->coordinates();
		reference_value = coord_sys.referenceValue();
		std::cout << "Reference pixel: " << reference_value[0] << ", "
				<< reference_value[1] << ", "
				<< reference_value[2] << std::endl;
	}

	if (false)   // Show data
	{
		IPosition pos(shape.size(), 0);
		for (int j=0; j<shape(1); ++j)
		{
			pos[1] = j;
			for (int i=0; i<shape(0); ++i)
			{
				pos[0] = i;
				std::cout << "i: " << i << ", j: " << j << ", val: " << imagePtr->getAt(pos) << std::endl;
			}
		}
	}
	return 0;
}
