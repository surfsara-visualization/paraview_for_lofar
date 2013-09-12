/*=========================================================================

  Program:   Visualization Toolkit
  Module:    LofarSliceImage.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME LofarSliceImage - Computes the gradient vector.
// .SECTION Description
// LofarSliceImage computes the gradient vector of an image.  The
// vector results are stored as scalar components. The Dimensionality
// determines whether to perform a 2d or 3d gradient. The default is
// two dimensional XY gradient.  OutputScalarType is always
// double. Gradient is computed using central differences.

#ifndef __LOFARSLICEIMAGE_H
#define __LOFARSLICEIMAGE_H

#include "vtkImagingGeneralModule.h" // For export macro
#include "vtkThreadedImageAlgorithm.h"
#include <vector>

class vtkParametricFunctionSource;
class vtkPolyDataAlgorithm;
class vtkImplicitFunction;
class vtkParametricSpline;

class VTKIMAGINGGENERAL_EXPORT LofarSliceImage : public vtkThreadedImageAlgorithm
{
public:
	static LofarSliceImage *New();
	vtkTypeMacro(LofarSliceImage,vtkThreadedImageAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	virtual void SetCurve(vtkParametricSpline *);
	virtual void SetCurve(vtkPolyDataAlgorithm *);

	virtual unsigned long GetMTime();

	void get(int delta, int &x, int &y) const;

	// Returns whether the point lies in the input domain.
	bool WorldToPixelCoordinates(const double world_pos[3], int pixel_pos[3]);

protected:

	LofarSliceImage();
	~LofarSliceImage() {};

	virtual int RequestInformation (vtkInformation*,
			vtkInformationVector**,
			vtkInformationVector*);
	virtual int RequestUpdateExtent(vtkInformation*,
			vtkInformationVector**,
			vtkInformationVector*);
	virtual int RequestData(vtkInformation*,
			vtkInformationVector**,
			vtkInformationVector*);

	void ThreadedRequestData(vtkInformation*,
			vtkInformationVector**,
			vtkInformationVector*,
			vtkImageData*** inData,
			vtkImageData** outData,
			int outExt[6],
			int threadId);

private:
	LofarSliceImage(const LofarSliceImage&);  // Not implemented.
	void operator=(const LofarSliceImage&);  // Not implemented.

	template <class T>
	void LofarSliceImageExecute(LofarSliceImage *self,
			vtkImageData *inData, T *inPtr,
			vtkImageData *outData, double *outPtr,
			int outExt[6], int id);

	enum ClipType {
		LINE,
		SPLINE,
	} m_clip_type;
	vtkParametricSpline *Spline;
	vtkPolyDataAlgorithm *LineSegment;

	int m_pt0_in_pixels[3];
	int m_pt1_in_pixels[3];
	int m_n_pixels;

	double m_input_origin[3];
	double m_input_spacing[3];
	int m_input_extent[6];

	int m_output_extent[6];

private:
	struct SplinePoint {
		double pt[3];
		double cumm_length;
	};
	std::vector<SplinePoint> LinearizeClippingFunction();

};

#endif



