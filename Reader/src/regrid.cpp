#include <casa/Utilities/CountedPtr.h>
#include <images/Images/FITSImage.h>
#include <images/Images/ImageOpener.h>
#include <images/Images/ImageProxy.h>

int main(int argc, char* argv[])
{
	using namespace casa;

	//ImageOpener::registerOpenImageFunction(ImageOpener::FITS, &FITSImage::openFITSImage);

	const char* filename = "~/Projects/LOFAR/Data/h5img/L70181_SBG001_sky.h5";
	if (argc > 1)
		filename = argv[1];

	ImageProxy image_proxy(filename, String(""), vector<ImageProxy>());
	ImageInterface<float>* imagePtr = dynamic_cast<ImageInterface<float>*>(image_proxy.getLattice());

	if (imagePtr == NULL)
	{
		std::cout << "1. Couldn't load the fits file: " << filename << "." << std::endl;
		return 1;
	}

	IPosition shape = imagePtr->shape();
	Vector<Double> increment = imagePtr->coordinates().increment();

	CoordinateSystem coord_sys = imagePtr->coordinates();
	Vector<Double> orig_reference_value = coord_sys.referenceValue();

	std::cout << "shape:           " << shape[0] << ", " << shape[1] << std::endl;
	std::cout << "increment:       " << increment[0] << ", " << increment[1] << std::endl;
	std::cout << "reference_value: " << orig_reference_value[0] << ", " << orig_reference_value[1] << std::endl;

	for (int x=0; x<10; ++x) {
		for (int y=0; y<10; ++y) {
			// Update the reference value
			Vector<Double> new_reference_value(3);
			for (int i=0; i<3; ++i)
				new_reference_value[i] = orig_reference_value[i];
			new_reference_value[0] += x * (shape[0]-1) * increment[0];
			new_reference_value[1] += y * (shape[1]-1) * increment[1];

			// Change the reference value
			coord_sys.setReferenceValue(new_reference_value);
			imagePtr->setCoordinateInfo(coord_sys);

			// Write file
			char new_filename[80];
			sprintf(new_filename, "offset_%d_%d.fits", x, y);
			std::cout << new_filename << std::endl;
			image_proxy.toFits(new_filename);
		}
	}

	return 0;
}
