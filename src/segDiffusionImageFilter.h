/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkGeodesicActiveContourLevelSetImageFilter.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __segDiffusionImageFilter_h
#define __segDiffusionImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkNeighborhoodIterator.h>
#include "itkImageLinearIteratorWithIndex.h"
#include <itkArray.h>
#include <itkPointSet.h>
#include <itkCastImageFilter.h>

#include "radDefinition.h"

namespace itk {

/** \class GeodesicActiveContourLevelSetImageFilter
 * \brief Segments structures in images based on a user supplied edge potential map.
 *
 * \par IMPORTANT
 * The SegmentationLevelSetImageFilter class and the
 * GeodesicActiveContourLevelSetFunction class contain additional information necessary
 * to gain full understanding of how to use this filter.
 *
 * \par OVERVIEW
 * This class is a level set method segmentation filter. An initial contour
 * is propagated outwards (or inwards) until it ''sticks'' to the shape boundaries.
 * This is done by using a level set speed function based on a user supplied
 * edge potential map.
 *
 * \par INPUTS
 * This filter requires two inputs.  The first input is a initial level set.
 * The initial level set is a real image which contains the initial contour/surface
 * as the zero level set. For example, a signed distance function from the initial
 * contour/surface is typically used. Unlike the simpler ShapeDetectionLevelSetImageFilter
 * the initial contour does not have to lie wholly within the shape to be segmented.
 * The intiial contour is allow to overlap the shape boundary. The extra advection term
 * in the update equation behaves like a doublet and attracts the contour to the boundary.
 * This approach for segmentation follows that of Caselles et al (1997).
 *
 * \par
 * The second input is the feature image.  For this filter, this is the edge
 * potential map. General characteristics of an edge potential map is that
 * it has values close to zero in regions near the edges and values close
 * to one inside the shape itself. Typically, the edge potential map is compute
 * from the image gradient, for example:
 *
 * \f[ g(I) = 1 / ( 1 + | (\nabla * G)(I)| ) \f]
 * \f[ g(I) = \exp^{-|(\nabla * G)(I)|} \f]
 * 
 * where \f$ I \f$ is image intensity and
 * \f$ (\nabla * G) \f$ is the derivative of Gaussian operator. 
 *
 * \par
 * See SegmentationLevelSetImageFilter and SparseFieldLevelSetImageFilter 
 * for more information on Inputs.
 *
 * \par PARAMETERS
 * The PropagationScaling parameter can be used to switch from propagation outwards
 * (POSITIVE scaling parameter) versus propagating inwards (NEGATIVE scaling 
 * parameter). 
 *
 * This implementation allows the user to set the weights between the propagation, advection
 * and curvature term using methods SetPropagationScaling(), SetAdvectionScaling(),
 * SetCurvatureScaling(). In general, the larger the CurvatureScaling, the smoother the
 * resulting contour. To follow the implementation in Caselles et al paper,
 * set the PropagationScaling to \f$ c \f$ (the inflation or ballon force) and
 * AdvectionScaling and CurvatureScaling both to 1.0.
 *
 * \par OUTPUTS
 * The filter outputs a single, scalar, real-valued image.
 * Negative values in the output image represent the inside of the segmented region
 * and positive values in the image represent the outside of the segmented region.  The
 * zero crossings of the image correspond to the position of the propagating
 * front.
 *
 * \par
 * See SparseFieldLevelSetImageFilter and
 * SegmentationLevelSetImageFilter for more information.
 *
 * \par REFERENCES
 * \par  
 *    "Geodesic Active Contours",
 *    V. Caselles, R. Kimmel and G. Sapiro.
 *    International Journal on Computer Vision,
 *    Vol 22, No. 1, pp 61-97, 1997
 *
 * \sa SegmentationLevelSetImageFilter
 * \sa GeodesicActiveContourLevelSetFunction
 * \sa SparseFieldLevelSetImageFilter 
 *
 * \ingroup LevelSetSegmentation
 */
template< class TInputImage, class TOutputImage>
class segDiffusionImageFilter
: public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs */
	typedef segDiffusionImageFilter Self;
	typedef ImageToImageFilter< TInputImage, TOutputImage> Superclass;
	typedef SmartPointer<Self>        Pointer;
	typedef SmartPointer<const Self>  ConstPointer;

	/** Inherited typedef from the superclass. */
	typedef typename Superclass::OutputImageType  OutputImageType;
	typedef typename Superclass::InputImageType InputImageType;
	
	/** Image dimension enumeration. */
	itkStaticConstMacro (ImageDimension, unsigned int, TInputImage::ImageDimension);

	/** Define double type image for level-set propagation. */
	typedef double																	InternalRealPixelType;
	typedef Image<InternalRealPixelType, itkGetStaticConstMacro(ImageDimension)>	InternalRealImageType;
	typedef typename OutputImageType::PixelType										OutputPixelType;

	typedef typename InputImageType::Pointer										InputImagePointer;
	typedef typename OutputImageType::Pointer										OutputImagePointer;
	typedef typename InternalRealImageType::Pointer									InternalRealImagePointer;

	typedef itk::ImageRegionIteratorWithIndex<InputImageType>						InputImageIteratorType;
	typedef itk::ImageRegionIteratorWithIndex<OutputImageType>						OutputImageIteratorType;
	typedef itk::ImageRegionIteratorWithIndex<InternalRealImageType>				InternalRealImageIteratorType;

	typedef itk::ImageLinearIteratorWithIndex< InternalRealImageType >				InternalRealLinearIteratorType;

	typedef typename TInputImage::IndexType											InputImageIndexType;
	typedef typename TOutputImage::IndexType										OutputImageIndexType;
	typedef typename InternalRealImageType::IndexType								InternalRealImageIndexType;

	 /** Type definition for the input image region type. */
	typedef typename TInputImage::RegionType			InputImageRegionType;
	typedef typename TOutputImage::RegionType			OutputImageRegionType;
	typedef typename InternalRealImageType::RegionType	InternalRealImageRegionType;

	/** Used for AOS computation*/
	typedef itk::Array< InternalRealPixelType > AOSVectorType;

	/** Run-time type information (and related methods). */
	itkTypeMacro(segDiffusionImageFilter, ImageToImageFilter);

	/** Method for creation through the object factory */
	itkNewMacro(Self);

	/** Set the time step. */
	itkSetClampMacro( TimeStep, double, 0.1, NumericTraits<unsigned int>::max() );

	/** Set the number of iterations. */
	itkSetClampMacro( NumberOfIterations, unsigned int, 1, NumericTraits<unsigned int>::max() );

	/** Set the weight of region terms. */
	itkSetClampMacro( DiffusionWeight, double, 0, NumericTraits<unsigned int>::max() );

protected:
	~segDiffusionImageFilter() {}
	segDiffusionImageFilter();

	virtual void PrintSelf(std::ostream &os, Indent indent) const override;

	segDiffusionImageFilter(const Self &); // purposely not implemented
	void operator=(const Self&); //purposely not implemented

	/** Overridden from Superclass to handle the case when PropagationScaling is zero.*/
	void GenerateData() override;

private:

	double m_TimeStep; // time step 
	unsigned int m_NumberOfIterations; //number of iterations for level-set propagation
	double m_DiffusionWeight; //Diffusion type

	//aided functions - frequently used by other functions
	template <class TInput, class TOutput>
	typename TOutput::Pointer AllocateImage(typename TInput::Pointer in, typename TOutput::PixelType value);

	template <class TInput, class TOutput>
	typename TOutput::Pointer CopyImage(typename TInput::Pointer input, typename TOutput::Pointer output);

	template <class TInput, class TOutput>
	typename TOutput::Pointer CastImage(typename TInput::Pointer in);

	//functions
	void DiffuseImage();
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "segDiffusionImageFilter.txx"
#endif

#endif
