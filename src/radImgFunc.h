/*
 *  radImgFunc.h
 *  
 *
 *  Created by Jianfei Liu on 3/22/12.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef RADIMGFUNC_H
#define RADIMGFUNC_H

#include "radDefinition.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDirectory.h>
#include <itkAddImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkCastImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkImageMaskSpatialObject.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "itkChangeInformationImageFilter.h"

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAbsImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkSubtractImageFilter.h>

#include <locale>
#include <algorithm>

//IO functions

template <class TImageType>
typename TImageType::Pointer ReadImage(const char *fileName)
{
	typedef itk::ImageFileReader< TImageType >  ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(fileName);
	reader->SetReleaseDataFlag(true);
	
	try
	{
		reader->UpdateLargestPossibleRegion();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << err << std::endl;
	}
	return reader->GetOutput();
}

template <class TImageType>
void WriteImage(const char *fileName, typename TImageType::Pointer img)
{
	typedef itk::ImageFileWriter< TImageType > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fileName);
	writer->SetInput(img);
	writer->Update();
}

template <class TImage>
typename TImage::Pointer CreateTiffImage(typename TImage::Pointer img)
{
	typedef itk::ChangeInformationImageFilter< TImage >  FilterType;
	typename FilterType::Pointer filter = FilterType::New();
	typename TImage::SpacingType img_spacing;
	img_spacing.Fill(1);
	typename TImage::PointType img_origin;
	img_origin.Fill(0);

	filter->SetOutputSpacing(img_spacing);
	filter->ChangeSpacingOn();
	filter->SetOutputOrigin(img_origin);
	filter->ChangeOriginOn();
	filter->SetInput(img);
	filter->Update();
	return filter->GetOutput();
}

//=====================================================================================
//=====================================================================================
//=====================================================================================

template <typename TImageType>
void ComputeIntensityRange(typename TImageType::Pointer img, 
                           typename TImageType::PixelType & low, 
                           typename TImageType::PixelType & high)
{
	typedef itk::MinimumMaximumImageCalculator<TImageType> MinimumMaximumImageCalculatorType;
	typename MinimumMaximumImageCalculatorType::Pointer intensityRangeCalc = MinimumMaximumImageCalculatorType::New();
	intensityRangeCalc->SetImage(img);
	intensityRangeCalc->ComputeMinimum();
	low = intensityRangeCalc->GetMinimum();
	intensityRangeCalc->ComputeMaximum();
	high = intensityRangeCalc->GetMaximum();
}

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer AllocateImage(typename TInputImage::Pointer input_img)
{
	typename TOutputImage::Pointer output_img = TOutputImage::New();
	output_img->SetRegions( input_img->GetLargestPossibleRegion() );
	output_img->SetSpacing( input_img->GetSpacing() );
	output_img->CopyInformation( input_img );
	output_img->Allocate();
    return output_img;
}

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer AllocateImage(typename TInputImage::RegionType region)
{
	typename TOutputImage::Pointer output_img = TOutputImage::New();
	output_img->SetRegions( region );
	output_img->Allocate();
    return output_img;
}

template <class TImageType>
typename TImageType::Pointer AllocateImage(typename TImageType::IndexType start, 
                   typename TImageType::SizeType size)
{
	typename TImageType::Pointer output_img = TImageType::New();
	typename TImageType::RegionType region;
    
	region.SetSize(size);
    region.SetIndex(start);
    output_img->SetRegions(region);
    output_img->Allocate();
    return output_img;
}

template<class TImage>
void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output)
{
	output->SetRegions(input->GetLargestPossibleRegion());
	output->Allocate();
 
	itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());
 
	while(!inputIterator.IsAtEnd())
    {
		outputIterator.Set(inputIterator.Get());
		++inputIterator;
		++outputIterator;
    }
}
 
template <class TImageType>
typename TImageType::Pointer ExtractImage(typename TImageType::Pointer input_img,
                  typename TImageType::IndexType &lower_index, 
                  typename TImageType::IndexType &upper_index)
{
	typename TImageType::IndexType desiredStart;
	typename TImageType::SizeType desiredSize;

	int i;
	for (i=0; i<input_img->GetImageDimension(); i++)
	{
		desiredStart[i] = lower_index[i];
		desiredSize[i] = upper_index[i]-lower_index[i];
	}
	
	typename TImageType::RegionType desiredRegion(desiredStart, desiredSize);
	typedef itk::ExtractImageFilter<TImageType, TImageType> ExtractImageFilterType;
	typename ExtractImageFilterType::Pointer extractFilter = ExtractImageFilterType::New();
	extractFilter->SetExtractionRegion(desiredRegion);
	extractFilter->SetInput(input_img);
	extractFilter->Update();
	return extractFilter->GetOutput();
}

template <class TImageType>
typename TImageType::Pointer ExtractImage(typename TImageType::Pointer input_img,
                  typename TImageType::RegionType & img_region)
{
	typedef itk::ExtractImageFilter<TImageType, TImageType> ExtractImageFilterType;
	typename ExtractImageFilterType::Pointer extractFilter = ExtractImageFilterType::New();
	extractFilter->SetExtractionRegion(img_region);
	extractFilter->SetInput(input_img);
	extractFilter->Update();
	return extractFilter->GetOutput();
}

template <typename TImageType>
typename TImageType::Pointer DuplicateImage(typename TImageType::Pointer input_img)
{
	typedef itk::ImageDuplicator< TImageType > DuplicatorType;
	typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(input_img);
	duplicator->Update();
	return duplicator->GetOutput();
}

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer CastImage(typename TInputImage::Pointer input_img)
{
	typedef itk::CastImageFilter< TInputImage, TOutputImage > CastFilterType;
	typename CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(input_img);
	castFilter->SetReleaseDataFlag(true);
	castFilter->Update();
	return castFilter->GetOutput();
}

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer MultipleImageWithConstant(typename TInputImage::Pointer input_img, float val)
{
	typedef itk::MultiplyImageFilter<TInputImage, TOutputImage, TOutputImage> MultiplyImageFilterType;
	typename MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
    multiplyImageFilter->SetInput(input_img);
	multiplyImageFilter->SetConstant(val);
	//multiplyImageFilter->SetReleaseDataFlag(true);
	multiplyImageFilter->Update();
	return multiplyImageFilter->GetOutput();
}

template <class TInputImage1, class TInputImage2, class TOutputImage>
typename TOutputImage::Pointer MultipleImages(typename TInputImage1::Pointer input_img1, 
                    typename TInputImage2::Pointer input_img2)
{
	typedef itk::MultiplyImageFilter <TInputImage1, TInputImage2, TOutputImage > MultiplyImageFilterType;
	typename MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New();
	multiplyFilter->SetInput1(input_img1);
	multiplyFilter->SetInput2(input_img2);
	multiplyFilter->SetReleaseDataFlag(true);
	multiplyFilter->Update();
	return multiplyFilter->GetOutput();
}

template <class TInputImage1, class TInputImage2, class TOutputImage>
typename TOutputImage::Pointer SubtractImages(typename TInputImage1::Pointer input_img1, 
                    typename TInputImage2::Pointer input_img2)
{
	typedef itk::SubtractImageFilter <TInputImage1, TInputImage2, TOutputImage > SubtractImageFilterType;
	typename SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New();
	subtractFilter->SetInput1(input_img1);
	subtractFilter->SetInput2(input_img2);
	//subtractFilter->SetReleaseDataFlag(true);
	subtractFilter->Update();
	return subtractFilter->GetOutput();
}

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer AbsoluteImage(typename TInputImage::Pointer input_img)
{
	typedef itk::AbsImageFilter<TInputImage, TOutputImage> AbsImageFilterType;
	typename AbsImageFilterType::Pointer abs_filter = AbsImageFilterType::New();
    abs_filter->SetInput(input_img);
	abs_filter->Update();
	return abs_filter->GetOutput();
}

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer RescaleImageIntensity(typename TInputImage::Pointer input_img,
                           typename TOutputImage::PixelType low,
                           typename TOutputImage::PixelType high)
{
	typedef itk::RescaleIntensityImageFilter< TInputImage, TOutputImage > RescaleFilterType;
	typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(input_img);
	rescaleFilter->SetOutputMinimum(low);
	rescaleFilter->SetOutputMaximum(high);
	rescaleFilter->Update();
	return rescaleFilter->GetOutput();
}

template <class TInputImage1, class TInputImage2, class TOutputImage, class TTransform>
typename TOutputImage::Pointer ResampleImage(typename TInputImage1::Pointer fixed_img,
				   typename TInputImage2::Pointer moving_img, 
				   typename TTransform::Pointer transform)
{
    typedef itk::ResampleImageFilter< TInputImage2, TOutputImage >    ResampleFilterType;
	typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
	resample->SetTransform( transform );
	resample->SetInput( moving_img );

	resample->SetSize(    fixed_img->GetLargestPossibleRegion().GetSize() );
	resample->SetOutputOrigin(  fixed_img->GetOrigin() );
	resample->SetOutputSpacing( fixed_img->GetSpacing() );
	resample->SetOutputDirection( fixed_img->GetDirection() );
	resample->SetDefaultPixelValue( 0 );
	//resample->SetReleaseDataFlag(true);
	resample->UpdateLargestPossibleRegion();
	return resample->GetOutput();
}


template <typename TImageType, unsigned int VImageDimension>
void DetermineBoundingBox(typename TImageType::Pointer input_img, typename TImageType::RegionType & output_region)	
{
	typedef itk::ImageMaskSpatialObject<VImageDimension> ImageMaskSpatialObject;
	typename ImageMaskSpatialObject::Pointer maskSO = ImageMaskSpatialObject::New();
	maskSO->SetReleaseDataFlag(true);
	maskSO->SetImage ( input_img );
	output_region = maskSO->GetAxisAlignedBoundingBoxRegion();
	std::cout << "Bounding Box Region: " << output_region << std::endl;
}


template <class TITKImage, class TType>
void ConvertITKImageToVTKImage(typename TITKImage::Pointer input_img, int number_of_components, 
                           vtkSmartPointer<vtkImageData> output_img, int vtk_TType)
{
	if (input_img->GetImageDimension() == 3)
	{
		output_img->SetDimensions(input_img->GetLargestPossibleRegion().GetSize()[0], 
			input_img->GetLargestPossibleRegion().GetSize()[1], 
			input_img->GetLargestPossibleRegion().GetSize()[2]);
		output_img->SetSpacing(input_img->GetSpacing()[0], input_img->GetSpacing()[1], input_img->GetSpacing()[2]);
		output_img->SetOrigin(input_img->GetOrigin()[0], input_img->GetOrigin()[1], input_img->GetOrigin()[2]);
	}
	else
	{
		output_img->SetDimensions(input_img->GetLargestPossibleRegion().GetSize()[0], 
			input_img->GetLargestPossibleRegion().GetSize()[1], 1);
		output_img->SetSpacing(input_img->GetSpacing()[0], input_img->GetSpacing()[1], 1.0);
		output_img->SetOrigin(input_img->GetOrigin()[0], input_img->GetOrigin()[1], 0);
	}

	output_img->AllocateScalars(vtk_TType, number_of_components);
    
    typedef itk::ImageRegionIteratorWithIndex<TITKImage>   TIteratorType;
    TIteratorType it(input_img, input_img->GetLargestPossibleRegion());
    //int index = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
    {
		TType * pixel;
		if (input_img->GetImageDimension() == 3)
			pixel = static_cast<TType *>(output_img->GetScalarPointer(it.GetIndex()[0],it.GetIndex()[1],it.GetIndex()[2]));
		else
			pixel = static_cast<TType *>(output_img->GetScalarPointer(it.GetIndex()[0],it.GetIndex()[1],0));

		/*if (number_of_components == 1)
			pixel[0] = it.Get();
		else*/
		{
			for (int i=0; i<number_of_components; i++)
				pixel[i] = it.Get()[i];
		}
    }
}

template <class TITKImage, class TType>
typename TITKImage::Pointer convertVtkToItk(vtkImageData* vtkInput)
{
	typename TITKImage::Pointer itkoutput = TITKImage::New();
	typename TITKImage::RegionType region;
	typename TITKImage::SizeType size;
	typename TITKImage::IndexType start;
	typename TITKImage::SpacingType spacing;
	typename TITKImage::PointType origin;

	for (int i=0; i<3; i++)
	{
		size[i] = vtkInput->GetDimensions()[i];
		start[i] = 0;
		spacing[i] = vtkInput->GetSpacing()[i];
		origin[i] = vtkInput->GetOrigin()[i];
	}
    
	region.SetSize(size);
    region.SetIndex(start);
    itkoutput->SetRegions(region);
	itkoutput->SetSpacing( spacing );
	itkoutput->SetOrigin(origin);
    itkoutput->Allocate();
    
    typedef itk::ImageRegionIteratorWithIndex<TITKImage>   TIteratorType;
    TIteratorType it(itkoutput, itkoutput->GetLargestPossibleRegion());
    //int index = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) 
    {
		TType * pixel = static_cast<TType *>(vtkInput->GetScalarPointer(it.GetIndex()[0],it.GetIndex()[1],it.GetIndex()[2]));
		it.Set(pixel[0]);
    }

	return itkoutput;
}

template <typename TImageType>
typename TImageType::Pointer BinaryDilate(typename TImageType::Pointer input_img, unsigned int radius)
{
    typedef itk::BinaryBallStructuringElement<typename TImageType::PixelType, TImageType::ImageDimension>  StructuringElementType;
	StructuringElementType se;

	/*typename TImageType::SizeType rad;
	rad.Fill(radius);*/

	se.SetRadius(radius);
	se.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter<TImageType, TImageType, StructuringElementType> BinaryDilateImageFilterType;
	typename BinaryDilateImageFilterType::Pointer dilater = BinaryDilateImageFilterType::New();
	dilater->SetInput(input_img);
	dilater->SetKernel(se);
	dilater->SetForegroundValue(255);
	dilater->SetBackgroundValue(0);
	dilater->Update();

	return dilater->GetOutput();
}

template <typename TInputImage, typename TOutputImage>
typename TOutputImage::Pointer DistanceTransform(typename TInputImage::Pointer input_img, bool inside_flag)
{
	typedef itk::SignedMaurerDistanceMapImageFilter<TInputImage, TOutputImage> DistanceTransformType;
	typename DistanceTransformType::Pointer distance_filter = DistanceTransformType::New();
	distance_filter->SetInput(input_img);
	distance_filter->SetInsideIsPositive(inside_flag);
	distance_filter->SetUseImageSpacing(1);

	//distance_filter->SetReleaseDataFlag(true);
	distance_filter->Update();
	return distance_filter->GetOutput();
}

template <typename TImageType>
typename TImageType::Pointer BinaryOpening(typename TImageType::Pointer input_img, unsigned int radius)
{
    typedef itk::BinaryBallStructuringElement<typename TImageType::PixelType, TImageType::ImageDimension>  StructuringElementType;
	StructuringElementType se;

	/*typename TImageType::SizeType rad;
	rad.Fill(radius);*/

	se.SetRadius(radius);
	se.CreateStructuringElement();

	typedef itk::BinaryMorphologicalOpeningImageFilter<TImageType, TImageType, StructuringElementType> BinaryOpeningImageFilterType;
	typename BinaryOpeningImageFilterType::Pointer opening = BinaryOpeningImageFilterType::New();
	opening->SetInput(input_img);
	opening->SetKernel(se);
	opening->SetForegroundValue(255);
	opening->Update();

	return opening->GetOutput();
}

#endif // RADIMGFUNC_H
