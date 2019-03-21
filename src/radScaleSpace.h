/***********************************************************************
 **  radScaleSpace.h
 ***********************************************************************/

#if !defined(radScaleSpace_H)
#define radScaleSpace_H 1

#include "radImgFunc.h"
#include <itkRecursiveGaussianImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include <itkNeighborhoodIterator.h>

namespace scalespace 
{
	//supplement function
	template <class TInput, class TOutput>
	TOutput Math_Angle(TInput dx, TInput dy);
	template <class TInput>
	unsigned char Math_Direction(TInput a, TInput b);
	template <class TInput>
	unsigned char Math_Direction(TInput dd);


    // "smooth" is a special case of "gaussian" (order_x=order_y=order_z=0)
    // However for better boundary accuracy, please use smooth 
	template <class TInputImage, class TOutputImage>
	typename TOutputImage::Pointer SmoothImage(typename TInputImage::Pointer img, float scale, 
		const float ratio = 1.0, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
	typename TOutputImage::Pointer SmoothImage(typename TInputImage::Pointer img, float scale_x, 
		float scale_y, float scale_z, bool scale_flag = false);

	// Partial derivatives
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DyImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DzImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);

	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxxImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DyyImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DzzImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxyImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxyImage(typename TInputImage::Pointer img, float x_scale, 
		float y_scale, float z_scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxzImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxzImage(typename TInputImage::Pointer img, float x_scale, 
		float y_scale, float z_scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DyzImage(typename TInputImage::Pointer img, float scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DyzImage(typename TInputImage::Pointer img, float x_scale, 
		float y_scale, float z_scale, bool scale_flag = false);

	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxyzImage(typename TInputImage::Pointer img, float x_scale, 
		float y_scale, float z_scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxxxImage(typename TInputImage::Pointer img, float x_scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DyyyImage(typename TInputImage::Pointer img, float y_scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DzzzImage(typename TInputImage::Pointer img, float z_scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxxyImage(typename TInputImage::Pointer img, float x_scale, 
		float y_scale, bool scale_flag = false);
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxxzImage(typename TInputImage::Pointer img, float x_scale, 
		float z_scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxyyImage(typename TInputImage::Pointer img, float x_scale, 
		float y_scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DyyzImage(typename TInputImage::Pointer img, float y_scale, 
		float z_scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DxzzImage(typename TInputImage::Pointer img, float x_scale, 
		float z_scale, bool scale_flag = false);
    template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer DyzzImage(typename TInputImage::Pointer img, float y_scale, 
		float z_scale, bool scale_flag = false);
    
	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer ComputeGradientMagnitude(typename TInputImage::Pointer img);

	template <class TInputImage, class TOutputImage>
	typename TOutputImage::Pointer ComputeGradientMagnitude(typename TInputImage::Pointer Ix, 
		typename TInputImage::Pointer Iy);

	template <class TInputImage, class TOutputImage>
	typename TOutputImage::Pointer ComputeGradientMagnitude(typename TInputImage::Pointer Ix, 
		typename TInputImage::Pointer Iy, typename TInputImage::Pointer Iz);

	template <class TInputImage, class TOutputImage>
	typename TOutputImage::Pointer ComputeLaplacian(typename TInputImage::Pointer img, float scale, bool scale_flag = false);

	template <class TInputImage, class TOutputImage>
	typename TOutputImage::Pointer ComputeAngleInRadian(typename TInputImage::Pointer Ix, 
		typename TInputImage::Pointer Iy);

    template <class TInputImage, class TOutputImage>
	typename TOutputImage::Pointer ComputeAngle(typename TInputImage::Pointer Ix, 
		typename TInputImage::Pointer Iy);

	template <class TInputImage, class TOutputImage>
	typename TOutputImage::Pointer ComputeDirection(typename TInputImage::Pointer Ix, 
		typename TInputImage::Pointer Iy);

	template <class TInputImage, class TOutputImage>
    typename TOutputImage::Pointer ComputeDirection(typename TInputImage::Pointer Ix, 
		typename TInputImage::Pointer Iy);

} // End of namespace scalespace

#include "radScaleSpace.txx"

#endif // radScaleSpace_H
