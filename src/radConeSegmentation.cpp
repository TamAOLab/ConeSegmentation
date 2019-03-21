/*
 *  radConeSegmentation.cpp
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "radconesegmentation.h"
#include "segDiffusionImageFilter.h"
#include "itkVoronoiDiagram2DGenerator.h"
#include <QDebug>

#include "QuickView.h"
#include "itkListSample.h"
#include "itkMeanSampleFilter.h"
#include "itkCovarianceSampleFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkZeroCrossingImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"

#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "vnl/vnl_math.h"
#include <vnl/vnl_matrix.h>
#include "vnl/algo/vnl_determinant.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include <vnl/vnl_vector.h>
#include <QMessageBox>

using namespace scalespace;

//used for scale selection
const unsigned int Number_Of_Scale_Levels = 15;
const double  Scale_Interval = 1.2;
const double  Initial_Scale = 0.5;

template <class T>
static inline T EnforceRange(const T & x,const int & MaxValue) {return std::min<T>(std::max<T>(x, 0), MaxValue-1);};

template <typename TVal>
void ComputeMeanAndVariance(vector<TVal> & vecs, TVal & mean, TVal & deviation)
{
	const unsigned int MeasurementVectorLength = 1;
	typedef itk::Vector< TVal, MeasurementVectorLength > MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
    typename SampleType::Pointer sample = SampleType::New();
	sample->SetMeasurementVectorSize( MeasurementVectorLength );

    for (typename vector<TVal>::iterator it = vecs.begin(); it != vecs.end(); ++it)
	{
		MeasurementVectorType mv;
		mv[0] = *it;
		sample->PushBack( mv );
	}

	typedef itk::Statistics::MeanSampleFilter< SampleType > MeanAlgorithmType;
    typename MeanAlgorithmType::Pointer meanAlgorithm = MeanAlgorithmType::New();
	meanAlgorithm->SetInput( sample );
	meanAlgorithm->Update();
	mean = meanAlgorithm->GetMean()[0];

	typedef itk::Statistics::CovarianceSampleFilter< SampleType > CovarianceAlgorithmType;
    typename CovarianceAlgorithmType::Pointer covarianceAlgorithm = CovarianceAlgorithmType::New();
	covarianceAlgorithm->SetInput( sample );
	covarianceAlgorithm->Update();
	deviation = sqrt(covarianceAlgorithm->GetCovarianceMatrix()(0, 0));
}

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer ComputeZeroCrossing(typename TInputImage::Pointer input_img)
{
	typedef itk::ZeroCrossingImageFilter< TInputImage, TOutputImage >  filterType;
    typename filterType::Pointer filter = filterType::New();
	filter->SetInput( input_img );
	filter->SetBackgroundValue(0);
	filter->SetForegroundValue(255);
	filter->Update();
	return filter->GetOutput();
}

//==============================================================================
//======================used for level set======================================
//==============================================================================
template<class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
	typedef CommandIterationUpdate   Self;
	typedef itk::Command             Superclass;
	typedef itk::SmartPointer<Self>  Pointer;
	itkNewMacro( Self );
protected:
	CommandIterationUpdate() {};

public:
	void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
		Execute( (const itk::Object *) caller, event);
    }
  
	void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
		const TFilter * filter = static_cast< const TFilter * >( object );
		if( typeid( event ) != typeid( itk::IterationEvent ) )
		{ return; }
    }
};

radConeSegmentation::radConeSegmentation()
{
	PairDistance = 30;
	HessianThreshold = 100;
	DarkConeThresholdRatio = 1.5;
	PairGACIterationNumber = 10;
	SingleGACIterationNumber = 10;

	BlobSizeThreshold = 0;
	ConeRadiusThreshold = ConeRoundnessThreshold = std::numeric_limits<double>::max();
}

radConeSegmentation::~radConeSegmentation()
{ 
}

void radConeSegmentation::SetParameters(double time_step, unsigned int iteration_step, double weight)
{
	DiffusionTimeStep = time_step;
	NumberofSteps = iteration_step;
	DiffusionWeight = weight;
}

void radConeSegmentation::SetParameters(ConeSegmentationParameters & params)
{
	SetScale(ConeSegmentation::scale);
	DiffusionTimeStep = ConeSegmentation::diffusionTimeStep;
	NumberofSteps = ConeSegmentation::diffusionNumberOfSteps;
	DiffusionWeight = ConeSegmentation::diffusionWeight;
	SetHessianThreshold(params.HessianThreshold);
	SetDarkConeThresholdRatio(ConeSegmentation::DarkConeThresholdRatio);
	SetPairGACIterationNumber(params.GACIterationNumber);
	SetSingleGACIterationNumber(params.GACIterationNumber);
	SetPairDistance(ConeSegmentation::PairDistance);
}

void radConeSegmentation::LevelsetSegmentation(FloatImageType::Pointer input_img, BinaryImageType::Pointer hessian_img, 
											   int val, BinaryImageType::Pointer result_img)
{
	typedef itk::segDiffusionImageFilter<FloatImageType, FloatImageType> DiffusionFilter;
	DiffusionFilter::Pointer diffusion_filter = DiffusionFilter::New();
	diffusion_filter->SetTimeStep(DiffusionTimeStep);
	diffusion_filter->SetNumberOfIterations(NumberofSteps);
	diffusion_filter->SetDiffusionWeight(DiffusionWeight);
	diffusion_filter->SetInput(input_img);
	diffusion_filter->Update();
	FloatImageType::Pointer diffused_img = diffusion_filter->GetOutput();

	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< FloatImageType, FloatImageType >  GradientFilterType;
	typedef itk::SigmoidImageFilter< FloatImageType, FloatImageType >  SigmoidFilterType;
	typedef itk::BinaryThresholdImageFilter< FloatImageType, BinaryImageType > ThresholdingFilterType;

	GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
	SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();

	//gradientMagnitude->SetInput( smoothing->GetOutput() );
	gradientMagnitude->SetInput( diffused_img );
	sigmoid->SetInput( gradientMagnitude->GetOutput() );
	gradientMagnitude->SetSigma(  m_Scale  );

	gradientMagnitude->Update();

	FloatImageType::Pointer gradient_magn = gradientMagnitude->GetOutput();
	float min_magn, max_magn;
	ComputeIntensityRange<FloatImageType>(gradient_magn, min_magn, max_magn);
	sigmoid->SetOutputMinimum(  0.0  );
	sigmoid->SetOutputMaximum(  1.0  );
	const double alpha =  -0.5;
	const double beta  =  2.0;
	// Software Guide : BeginCodeSnippet
	sigmoid->SetAlpha( alpha );
	sigmoid->SetBeta(  beta  );

	BinaryImageType::Pointer dilate_img = BinaryDilate<BinaryImageType>(hessian_img, 2);
	FloatImageType::Pointer dist_img = DistanceTransform<BinaryImageType, FloatImageType>(dilate_img, false);
	
	double propagationScaling = 2.0;
	typedef  itk::GeodesicActiveContourLevelSetImageFilter< FloatImageType, FloatImageType >  GeodesicActiveContourFilterType;
	GeodesicActiveContourFilterType::Pointer geodesicActiveContour = GeodesicActiveContourFilterType::New();
	geodesicActiveContour->SetPropagationScaling( propagationScaling );
	geodesicActiveContour->SetCurvatureScaling( 1.0 );
	geodesicActiveContour->SetAdvectionScaling( 1.0 );
	geodesicActiveContour->SetMaximumRMSError( 0.005 );
	geodesicActiveContour->SetNumberOfIterations( PairGACIterationNumber );
	geodesicActiveContour->SetInput( dist_img );
	geodesicActiveContour->SetFeatureImage( sigmoid->GetOutput() );
	geodesicActiveContour->Update();
	FloatImageType::Pointer gac_img = geodesicActiveContour->GetOutput();

	thresholder->SetInput( geodesicActiveContour->GetOutput() );
	thresholder->SetLowerThreshold( -1000.0 );
	thresholder->SetUpperThreshold( 0.0  );
	thresholder->SetOutsideValue(  itk::NumericTraits< BinaryImageType::PixelType >::min()  );
	thresholder->SetInsideValue(  itk::NumericTraits< BinaryImageType::PixelType >::max() );
	thresholder->Update();

	BinaryIteratorType seg_it(ConvexHullImage, ConvexHullImage->GetLargestPossibleRegion());
	BinaryIteratorType threshold_it(thresholder->GetOutput(), thresholder->GetOutput()->GetLargestPossibleRegion());
	BinaryIteratorType result_it(result_img, result_img->GetLargestPossibleRegion());
	for (seg_it.GoToBegin(), threshold_it.GoToBegin(), result_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it, ++threshold_it, ++result_it)
	{
		if (threshold_it.Get() != 0)
		{
			seg_it.Set(val);
			result_it.Set(255);
		}
	}
}

FloatImageType::Pointer radConeSegmentation::ComputeImageMagn(FloatImageType::Pointer img, float sigma)
{
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< FloatImageType, FloatImageType >  filterType;
	filterType::Pointer gradientFilter = filterType::New();
	gradientFilter->SetInput(img);
	gradientFilter->SetSigma(sigma);
	gradientFilter->SetNormalizeAcrossScale(true);
	gradientFilter->Update();
	return gradientFilter->GetOutput();
}

itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer radConeSegmentation::LabelObjectSegmentation(BinaryImageType::Pointer img)
{
	typedef itk::BinaryImageToShapeLabelMapFilter<BinaryImageType> BinaryImageToShapeLabelMapFilterType;
	BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToShapeLabelMapFilter->SetInput(img);
	binaryImageToShapeLabelMapFilter->Update();
	return binaryImageToShapeLabelMapFilter->GetOutput();
}

void radConeSegmentation::LabelBinaryImage(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType *labelObject, 
										   short val, ShortImageType::Pointer label_img)
{
	for (unsigned int j=0; j < labelObject->GetNumberOfPixels(); j++)
		label_img->SetPixel(labelObject->GetIndex(j), val);
}

void radConeSegmentation::GetRegionLabelsBelongingToCurrentObject(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType *labelObject, 
																  vector<unsigned int> & label_list)
{
	label_list.clear();
	for (unsigned int i=0; i<labelObject->GetNumberOfPixels(); i++)
	{
		if (std::find(label_list.begin(), label_list.end(), ObjectLabelImage->GetPixel(labelObject->GetIndex(i))) 
			== label_list.end() && ObjectLabelImage->GetPixel(labelObject->GetIndex(i)) != 0)
		{
			label_list.push_back(ObjectLabelImage->GetPixel(labelObject->GetIndex(i)));
		}
	}
}

void radConeSegmentation::ComputeSegmentationCenters(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer labeling, 
													 DoublePointArray & pts, vector<float> & radius)
{
	pts.clear();
	radius.clear();

	DoublePointType pt;
	for(unsigned int i = 0; i < labeling->GetNumberOfLabelObjects(); i++)
    {
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= labeling->GetNthLabelObject(i);
		// Output the bounding box (an example of one possible property) of the ith region
		//std::cout << "Object " << i << " has bounding box " << labelObject->GetBoundingBox() << std::endl;
		pt[0] = (labelObject->GetCentroid()[0]-InputImages.first->GetOrigin()[0])/InputImages.first->GetSpacing()[0];
		pt[1] = (labelObject->GetCentroid()[1]-InputImages.first->GetOrigin()[1])/InputImages.first->GetSpacing()[1];

		pts.push_back(pt);
		radius.push_back(labelObject->GetEquivalentSphericalRadius());
    }
}

int radConeSegmentation::SelectLeftBlob(int l_pt_id, int r_pt_id, vector<bool> &left_selection_flags)
{
	DoublePointType r_pt = RightSegmentationCenters.at(r_pt_id);
	DoublePointType l_pt = LeftSegmentationCenters.at(l_pt_id);
	DoublePointType tmp_pt;
	double dist = r_pt.EuclideanDistanceTo(l_pt);
	int selection_id = -1;
	// double min_dist = std::numeric_limits<double>::max();

	for (int i=12; i<=24; i++)
	{
		double angle = (double)i*10.0/180*PI;
		tmp_pt[0] = dist*cos(angle)+r_pt[0];
		tmp_pt[1] = dist*sin(angle)+r_pt[1];

		for (int j=0; j<left_selection_flags.size(); j++)
		{
			if (!left_selection_flags[j] && tmp_pt.EuclideanDistanceTo(LeftSegmentationCenters[j]) < LeftBlobRadius[j]
			&& LeftBlobRadius[l_pt_id] < LeftBlobRadius[j])
			{
				selection_id = j;
			}
		}
	}

	return selection_id;
}

void radConeSegmentation::DecidePairsAccordingToCenters()
{
	SegmentationPairs.clear();

	DoublePointType l_pt, r_pt;
	int dist, min_dist;
	unsigned int r_id;
	vector<bool> right_selection_flags(RightSegmentationCenters.size(), false);
	vector<int> dist_array(RightSegmentationCenters.size(), 10000);
	pair<unsigned int, unsigned int> point_link;
	//current radius
	double angle, min_angle;

	for (unsigned int i=0; i<LeftSegmentationCenters.size(); i++)
	{
		l_pt = LeftSegmentationCenters.at(i);
		r_id = std::numeric_limits<unsigned int>::max();
		min_dist = PairDistance*PairDistance+1;
		for (unsigned int j=0; j<RightSegmentationCenters.size(); j++)
		{
			r_pt = RightSegmentationCenters.at(j);
			angle = atan2(r_pt[1]-l_pt[1], r_pt[0]-l_pt[0]);
			if (angle > 75.0/180*PI || angle < -75.0/180*PI)
				continue;

			dist = (r_pt[0] -l_pt[0])*(r_pt[0] -l_pt[0]) + (r_pt[1] -l_pt[1])*(r_pt[1] -l_pt[1]);
			if (r_pt[0] > l_pt[0] && dist < PairDistance*PairDistance && min_dist > dist)
			{
				min_dist = dist;
				r_id = j;
				min_angle = angle;
			}
		}

		if (r_id != std::numeric_limits<unsigned int>::max())
		{
			//selection relaxation 
			for (unsigned int j=0; j<RightSegmentationCenters.size(); j++)
			{
				r_pt = RightSegmentationCenters.at(j);
				angle = atan2(r_pt[1]-l_pt[1], r_pt[0]-l_pt[0]);
			
				dist = (r_pt[0] -l_pt[0])*(r_pt[0] -l_pt[0]) + (r_pt[1] -l_pt[1])*(r_pt[1] -l_pt[1]);
				if (r_pt[0] > l_pt[0] && dist/min_dist < 1.2 && fabs(min_angle) > fabs(angle))
				{
					min_dist = dist;
					r_id = j;
				}
			}

			//next, slight expand min_dist and based on it to search for largest blob regions
			if (!right_selection_flags[r_id])
			{
				point_link.first = i;
				point_link.second = r_id;
				SegmentationPairs.push_back(point_link);
				dist_array[r_id] = min_dist;
				//cout << i << ": " << LeftSegmentationCenters.at(i) << ", " << RightSegmentationCenters.at(r_id) << std::endl;
			}
			else if (right_selection_flags[r_id] && dist_array[r_id] > min_dist)
			{
				dist_array[r_id] = min_dist;
				for (int k=0; k<SegmentationPairs.size(); k++)
				{
					if (SegmentationPairs[k].second == r_id)
					{
						point_link.first = SegmentationPairs[k].first;
						point_link.second = std::numeric_limits<unsigned int>::max();

						SegmentationPairs[k].first = i;
						SegmentationPairs.push_back(point_link);
					}
				}
			}
			else
			{
				point_link.first = i;
				point_link.second = std::numeric_limits<unsigned int>::max();
				SegmentationPairs.push_back(point_link);
			}
		
			right_selection_flags[r_id] = true;
		}
		else
		{
			point_link.first = i;
			point_link.second = r_id;
			SegmentationPairs.push_back(point_link);
		}
	}

	//locally adjust left point selection
	vector<bool> left_selection_flags(LeftSegmentationCenters.size(), false);
	for (vector< pair<unsigned int, unsigned int> >::iterator it = SegmentationPairs.begin(); it != SegmentationPairs.end(); ++it)
	{
		if (it->first != std::numeric_limits<unsigned int>::max() && it->second != std::numeric_limits<unsigned int>::max())
		{
			left_selection_flags[it->first] = true;
		}
	}

	for (vector< pair<unsigned int, unsigned int> >::iterator it = SegmentationPairs.begin(); it != SegmentationPairs.end(); ++it)
	{
		if (it->first != std::numeric_limits<unsigned int>::max() && it->second != std::numeric_limits<unsigned int>::max())
		{
			// double current_left_radius = LeftBlobRadius.at(it->first);
			r_pt = RightSegmentationCenters.at(it->second);
			l_pt = LeftSegmentationCenters.at(it->first);

			double angle = atan2(r_pt[1]-l_pt[1], r_pt[0]-l_pt[0]);
			if (angle > PI/18 || angle < - PI/18)
			{
				int selection_id = SelectLeftBlob(it->first, it->second, left_selection_flags);
				if (selection_id != -1)
				{
					left_selection_flags[it->first] = false;
					left_selection_flags[selection_id] = true;
					it->first = selection_id;
				}
			}
		}
	}

	for (vector<bool>::iterator it = right_selection_flags.begin(); it != right_selection_flags.end(); ++it)
	{
		if (!(*it))
		{
			point_link.first = std::numeric_limits<unsigned int>::max();
			point_link.second = it - right_selection_flags.begin();
			SegmentationPairs.push_back(point_link);
		}
	}

	//compute average dist and std deviation
	vector<pair<unsigned int, unsigned int> > tmp_pair_array;
	std::copy(SegmentationPairs.begin(), SegmentationPairs.end(), back_inserter(tmp_pair_array));

	vector<double> valid_radius_array;
	double radius_mean, radius_deviation;
	for (unsigned int i=0; i<dist_array.size(); i++)
	{
		if (dist_array[i] != 10000)
			valid_radius_array.push_back(dist_array[i]);
	}
	ComputeMeanAndVariance<double>(valid_radius_array, radius_mean, radius_deviation);
	SegmentationPairs.clear();
	for (unsigned int i=0; i<tmp_pair_array.size(); i++)
	{
		if (tmp_pair_array[i].first == std::numeric_limits<unsigned int>::max() || 
			tmp_pair_array[i].second == std::numeric_limits<unsigned int>::max())
		{
			SegmentationPairs.push_back(tmp_pair_array[i]);
		}
		else
		{
			r_id = tmp_pair_array[i].second;
			if (dist_array[r_id] > radius_mean + 3*radius_deviation)
			{
				tmp_pair_array[i].first = std::numeric_limits<unsigned int>::max();
				right_selection_flags[r_id] = false;
			}
			SegmentationPairs.push_back(tmp_pair_array[i]);
		}
	}

	//compute average pair distance
	AveragePairDistance = 0;
	unsigned int number_of_pairs = 0;
	for (vector< pair<unsigned int, unsigned int> >::iterator it = SegmentationPairs.begin(); it != SegmentationPairs.end(); ++it)
	{
		if (it->first != std::numeric_limits<unsigned int>::max() 
			&& it->second != std::numeric_limits<unsigned int>::max())
		{
			number_of_pairs++;
			AveragePairDistance += sqrt((LeftSegmentationCenters[it->first][0]-RightSegmentationCenters[it->second][0])
				*(LeftSegmentationCenters[it->first][0]-RightSegmentationCenters[it->second][0])
				+(LeftSegmentationCenters[it->first][1]-RightSegmentationCenters[it->second][1])
				*(LeftSegmentationCenters[it->first][1]-RightSegmentationCenters[it->second][1]));
		}
	}

	if (number_of_pairs != 0) AveragePairDistance /= number_of_pairs;
}

float radConeSegmentation::ComputeDistanceBetweenTwoBlobs(BinaryImageType::Pointer binary_img, FloatImageType::Pointer dist_img)
{
	float min_dist = 10000;
	BinaryIteratorType binary_it(binary_img, binary_img->GetLargestPossibleRegion());
	FloatIteratorType dist_it(dist_img, dist_img->GetLargestPossibleRegion());

	for (binary_it.GoToBegin(), dist_it.GoToBegin(); !binary_it.IsAtEnd(); ++binary_it, ++dist_it)
	{
		if (binary_it.Get() != 0 && min_dist > dist_it.Get())
		{
			min_dist = dist_it.Get();
		}
	}

	return min_dist;
}

void radConeSegmentation::FillGapBetweenBlobs(FloatImageType::Pointer left_dist, FloatImageType::Pointer right_dist, 
											  BinaryImageType::Pointer segmented_img, float dist_val, DoublePointType & left_ct, 
											  DoublePointType & right_ct, float left_rd, float right_rd)
{
	FloatIteratorType left_it(left_dist, left_dist->GetLargestPossibleRegion());
	FloatIteratorType right_it(right_dist, right_dist->GetLargestPossibleRegion());
	BinaryIteratorType seg_it(segmented_img, segmented_img->GetLargestPossibleRegion());
	float left_dist_val, right_dist_val;

	for (left_it.GoToBegin(), right_it.GoToBegin(), seg_it.GoToBegin(); !left_it.IsAtEnd(); ++left_it, ++right_it, ++seg_it)
	{
		if (left_it.Get() > 0 && right_it.Get() > 0 && fabs(left_it.Get()-right_it.Get()) < 1.5*dist_val)
		{
			left_dist_val = (left_it.GetIndex()[0]-left_ct[0])*(left_it.GetIndex()[0]-left_ct[0])
				+ (left_it.GetIndex()[1]-left_ct[1])*(left_it.GetIndex()[1]-left_ct[1]);
			right_dist_val = (right_it.GetIndex()[0]-right_ct[0])*(right_it.GetIndex()[0]-right_ct[0])
				+ (right_it.GetIndex()[1]-right_ct[1])*(right_it.GetIndex()[1]-right_ct[1]);

			if (left_dist_val < (left_rd+dist_val)*(left_rd+dist_val) && right_dist_val < (right_rd+dist_val)*(right_rd+dist_val))
				seg_it.Set(255);
		}
	}
}

void radConeSegmentation::FillBlobsBasedOnDistance(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer left_labeling, 
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer right_labeling)
{
	for (unsigned int i=0; i<SegmentationPairs.size(); i++)
	{
		if (SegmentationPairs[i].first != std::numeric_limits<unsigned int>::max()
			&& SegmentationPairs[i].second != std::numeric_limits<unsigned int>::max())
		{
			itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* left_labelObject 
				= left_labeling->GetNthLabelObject(SegmentationPairs[i].first);
			itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* right_labelObject 
				= right_labeling->GetNthLabelObject(SegmentationPairs[i].second);

			//assign label image in case of objects touch with each other
			LabelBinaryImage(left_labelObject, i+1, ObjectLabelImage);
			LabelBinaryImage(right_labelObject, i+1, ObjectLabelImage);

			BinaryImageType::RegionType paired_region;
			BinaryImageType::IndexType start_index, end_index;
			BinaryImageType::SizeType region_size;

			for (unsigned int j=0; j<2; j++)
			{
				start_index[j] = std::min<unsigned int>(left_labelObject->GetBoundingBox().GetIndex()[j], 
					right_labelObject->GetBoundingBox().GetIndex()[j]);
				end_index[j] = std::max<unsigned int>(left_labelObject->GetBoundingBox().GetIndex()[j]
					+left_labelObject->GetBoundingBox().GetSize()[j], right_labelObject->GetBoundingBox().GetIndex()[j]
					+right_labelObject->GetBoundingBox().GetSize()[j]);
				region_size[j] = end_index[j] - start_index[j];
			}
			paired_region.SetIndex(start_index);
			paired_region.SetSize(region_size);

			//create two binary images
			BinaryImageType::Pointer left_binary_img = AllocateImage<BinaryImageType, BinaryImageType>(paired_region);
			left_binary_img->FillBuffer(0);
			BinaryImageType::Pointer right_binary_img = AllocateImage<BinaryImageType, BinaryImageType>(left_binary_img);
			right_binary_img->FillBuffer(0);

			//assign two binary images
			BinaryIteratorType left_binary_it(left_binary_img, left_labelObject->GetBoundingBox());
			BinaryIteratorType right_binary_it(right_binary_img, right_labelObject->GetBoundingBox());
			BinaryIteratorType left_labeling_it(DarkSegmentationBlobs.first, left_labelObject->GetBoundingBox());
			BinaryIteratorType right_labeling_it(DarkSegmentationBlobs.second, right_labelObject->GetBoundingBox());

			for (left_binary_it.GoToBegin(), left_labeling_it.GoToBegin(); !left_binary_it.IsAtEnd(); 
				++left_binary_it, ++left_labeling_it)
				left_binary_it.Set(left_labeling_it.Get());

			for (right_binary_it.GoToBegin(), right_labeling_it.GoToBegin(); !right_binary_it.IsAtEnd(); 
				++right_binary_it, ++right_labeling_it)
				right_binary_it.Set(right_labeling_it.Get());

			FloatImageType::Pointer left_dist_img = DistanceTransform<BinaryImageType, FloatImageType>(left_binary_img, false);
			FloatImageType::Pointer right_dist_img = DistanceTransform<BinaryImageType, FloatImageType>(right_binary_img, false);

			float left_min_dist = ComputeDistanceBetweenTwoBlobs(left_binary_img, right_dist_img);
			float right_min_dist = ComputeDistanceBetweenTwoBlobs(right_binary_img, left_dist_img);

			BinaryImageType::Pointer local_seg = AllocateImage<BinaryImageType, BinaryImageType>(left_binary_img);
			local_seg->FillBuffer(0);
			FillGapBetweenBlobs(left_dist_img, right_dist_img, local_seg, std::max<float>(left_min_dist, right_min_dist), 
				LeftSegmentationCenters[SegmentationPairs[i].first], RightSegmentationCenters[SegmentationPairs[i].second], 
				LeftBlobRadius[SegmentationPairs[i].first], RightBlobRadius[SegmentationPairs[i].second]);

			BinaryIteratorType local_it(local_seg, local_seg->GetLargestPossibleRegion());
			BinaryIteratorType seg_it(ConvexHullImage, local_seg->GetLargestPossibleRegion());

			//bool, do extra checking
			bool flag_3 = false;
			for (local_it.GoToBegin(), seg_it.GoToBegin(); !local_it.IsAtEnd(); ++local_it, ++seg_it)
			{
				if (local_it.Get() != 0)
				{
					seg_it.Set(3);
					ObjectLabelImage->SetPixel(seg_it.GetIndex(), i+1);
					flag_3 = true;
				}
			}

			//check a neighborhood point is from other side, in the case of left and right regions touching with each other
			if (!flag_3)
			{
				BinaryImageType::SizeType radius;
				radius[0] = 1;
				radius[1] = 1;
				itk::NeighborhoodIterator<BinaryImageType> neighbor_it(radius, left_binary_img, left_binary_img->GetLargestPossibleRegion());
				for (neighbor_it.GoToBegin(); !neighbor_it.IsAtEnd(); ++neighbor_it)
				{
					if (neighbor_it.GetCenterPixel() != 0)
					{
						for (unsigned int k=0; k<neighbor_it.Size(); k++)
						{
							bool IsInBounds;
							neighbor_it.GetPixel(k, IsInBounds);
							if (IsInBounds)
							{
								if (right_binary_img->GetPixel(neighbor_it.GetIndex(k)))
								{
									ConvexHullImage->SetPixel(neighbor_it.GetIndex(k), 3);
									ObjectLabelImage->SetPixel(neighbor_it.GetIndex(k), i+1);
								}
							}
						}
					}
				}
			}

		}
	} 

}

double radConeSegmentation::ComputeDarkBlobHistogram()
{
	const unsigned int MeasurementVectorLength = 1;
	typedef itk::Vector< unsigned char, MeasurementVectorLength > MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > ListSampleType;
	ListSampleType::Pointer listSample = ListSampleType::New();
	listSample->SetMeasurementVectorSize( MeasurementVectorLength );
	
	MeasurementVectorType mv;
	BinaryIteratorType dark_it(DarkSegmentationBlobs.first, DarkSegmentationBlobs.first->GetLargestPossibleRegion());
	FloatIteratorType img_it(InputImages.first, InputImages.first->GetLargestPossibleRegion());
	for (dark_it.GoToBegin(), img_it.GoToBegin(); !dark_it.IsAtEnd(); ++dark_it, ++img_it)
	{
		if (dark_it.Get() != 0)
		{
			mv[0] = img_it.Get();
			listSample->PushBack(mv);
		}
	}

	dark_it = BinaryIteratorType(DarkSegmentationBlobs.second, DarkSegmentationBlobs.second->GetLargestPossibleRegion());
	img_it = FloatIteratorType(InputImages.second, InputImages.second->GetLargestPossibleRegion());
	for (dark_it.GoToBegin(), img_it.GoToBegin(); !dark_it.IsAtEnd(); ++dark_it, ++img_it)
	{
		if (dark_it.Get() != 0)
		{
			mv[0] = img_it.Get();
			listSample->PushBack(mv);
		}
	}

	HistogramType::SizeType size(NumberOfHistogramComponents);
	size.Fill(50);
	HistogramType::MeasurementVectorType lowerBound( NumberOfHistogramComponents );
	HistogramType::MeasurementVectorType upperBound( NumberOfHistogramComponents );
	lowerBound[0] = 0;
	upperBound[0] = 255;

	typedef itk::Statistics::SampleToHistogramFilter< ListSampleType, HistogramType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( listSample );
	filter->SetHistogramSize( size );
	filter->SetHistogramBinMinimum( lowerBound );
	filter->SetHistogramBinMaximum( upperBound );
	filter->Update();

	const HistogramType* histogram = filter->GetOutput();
	HistogramType::ConstIterator iter = histogram->Begin();
	double result_val = 0.;
	double accumulate_frequency = 0;
	while ( iter != histogram->End() )
    {
//		std::cout << "Measurement vectors = " << iter.GetMeasurementVector()
//			<< " frequency = " << iter.GetFrequency() << std::endl;
		accumulate_frequency += (double)iter.GetFrequency()/(double)histogram->GetTotalFrequency();
		DarkConeThreshold1 = iter.GetMeasurementVector()[0];
		if (accumulate_frequency >= 0.8)
		{
			result_val = iter.GetMeasurementVector()[0];
			break;
		}
		++iter;
    }
//	std::cout << "Size = " << histogram->Size() << std::endl;
//	std::cout << "Total frequency = " << histogram->GetTotalFrequency() << std::endl;

	return result_val*DarkConeThresholdRatio;
}

void radConeSegmentation::ExpandSegmentationPairs(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer left_labeling, 
												  itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer right_labeling)
{
	ComputeSegmentationCenters(left_labeling, LeftSegmentationCenters, LeftBlobRadius);
	ComputeSegmentationCenters(right_labeling, RightSegmentationCenters, RightBlobRadius);
	DecidePairsAccordingToCenters();
	FillBlobsBasedOnDistance(left_labeling, right_labeling);

	DarkConeThreshold = ComputeDarkBlobHistogram();
}

void radConeSegmentation::FillConvexHullImage(BinaryImageType::Pointer convex_hull_img, 
											  vtkSmartPointer<vtkPoints> contour_pts)
{
	//first transform contour_pts
	// BinaryImageType::PointType img_origin = convex_hull_img->GetOrigin();
	// BinaryImageType::SpacingType img_spacing = convex_hull_img->GetSpacing();
	BinaryImageType::IndexType pixelIndex;
	vector<float> tmp_contour_pts;
	DoublePointType point;

	tmp_contour_pts.clear();
	for (unsigned int j=0; j<contour_pts->GetNumberOfPoints(); j++)
	{
		for (unsigned int k=0; k<2; k++)
		{
			point[k] = contour_pts->GetPoint(j)[k];
		}
		convex_hull_img->TransformPhysicalPointToIndex(point, pixelIndex);

		for (unsigned int k=0; k<2; k++)
			tmp_contour_pts.push_back(pixelIndex[k]);
	}

	int  nodes, i1, j1, swap;
	int polyCorners = tmp_contour_pts.size()/2;
	int nodeX[1500];
	int IMAGE_RIGHT = convex_hull_img->GetLargestPossibleRegion().GetSize()[0], IMAGE_LEFT = 0;
	//  Loop through the rows of the image.
	for (unsigned int y=0; y<convex_hull_img->GetLargestPossibleRegion().GetSize()[1]; y++) 
	{
		//  Build a list of nodes.
		nodes=0; j1=polyCorners-1;
		for (i1=0; i1<polyCorners; i1++) 
		{
			if ((tmp_contour_pts[2*i1+1] < (double)y && tmp_contour_pts[2*j1+1] >= (double)y)
				|| (tmp_contour_pts[2*j1+1] < (double)y && tmp_contour_pts[2*i1+1] >= (double)y)) 
			{
				nodeX[nodes++]=(int)(tmp_contour_pts[2*i1]+(y-tmp_contour_pts[2*i1+1])/(tmp_contour_pts[2*j1+1]-tmp_contour_pts[2*i1+1])
				*(tmp_contour_pts[2*j1]-tmp_contour_pts[2*i1])); 
			}
    
			j1=i1; 
		}

		//  Sort the nodes, via a simple ?Bubble? sort.
		i1=0;
		while (i1<nodes-1) 
		{
			if (nodeX[i1]>nodeX[i1+1]) 
			{
				swap=nodeX[i1]; 
				nodeX[i1]=nodeX[i1+1]; 
				nodeX[i1+1]=swap; 
				if (i1) i1--; 
			}
			else 
			{
				i1++; 
			}
		}

		//  Fill the pixels between node pairs.
		for (i1=0; i1<nodes; i1+=2) 
		{
			if (nodeX[i1]>=IMAGE_RIGHT) 
				break;
			
			if (nodeX[i1+1]> IMAGE_LEFT ) 
			{
				if (nodeX[i1] < IMAGE_LEFT ) 
					nodeX[i1]=IMAGE_LEFT ;
      
				if (nodeX[i1+1]> IMAGE_RIGHT) 
					nodeX[i1+1]=IMAGE_RIGHT;
      
				for (j1=nodeX[i1]; j1<nodeX[i1+1]; j1++) 
				{
					pixelIndex[0] = j1;
					pixelIndex[1] = y;

					convex_hull_img->SetPixel(pixelIndex, 255);
				}
			}
		}
	}
}

void radConeSegmentation::FillConvexHullImage(BinaryImageType::Pointer convex_hull_img)
{
	//first transform contour_pts
	// BinaryImageType::PointType img_origin = convex_hull_img->GetOrigin();
	// BinaryImageType::SpacingType img_spacing = convex_hull_img->GetSpacing();
	BinaryImageType::IndexType pixelIndex;
	vector<float> tmp_contour_pts;
	DoublePointType point;

	short min_label, max_label;
	ComputeIntensityRange<ShortImageType>(ObjectLabelImage, min_label, max_label);

	short current_label;
	for (unsigned int i=0; i<ConvexHullPoints.size(); i++)
	{
		current_label = i + max_label + 1;
		tmp_contour_pts.clear();
		for (unsigned int j=0; j<ConvexHullPoints[i]->GetNumberOfPoints(); j++)
		{
			for (unsigned int k=0; k<2; k++)
			{
				point[k] = ConvexHullPoints[i]->GetPoint(j)[k];
			}
			convex_hull_img->TransformPhysicalPointToIndex(point, pixelIndex);

			for (unsigned int k=0; k<2; k++)
				tmp_contour_pts.push_back(pixelIndex[k]);
		}

		int  nodes, i1, j1, swap;
		int polyCorners = tmp_contour_pts.size()/2;
		int nodeX[1500];
		int IMAGE_RIGHT = convex_hull_img->GetLargestPossibleRegion().GetSize()[0], IMAGE_LEFT = 0;
		//  Loop through the rows of the image.
		for (unsigned int y=0; y<convex_hull_img->GetLargestPossibleRegion().GetSize()[1]; y++) 
		{
			//  Build a list of nodes.
			nodes=0; j1=polyCorners-1;
			for (i1=0; i1<polyCorners; i1++) 
			{
				if ((tmp_contour_pts[2*i1+1] < (double)y && tmp_contour_pts[2*j1+1] >= (double)y)
					|| (tmp_contour_pts[2*j1+1] < (double)y && tmp_contour_pts[2*i1+1] >= (double)y)) 
				{
					nodeX[nodes++]=(int)(tmp_contour_pts[2*i1]+(y-tmp_contour_pts[2*i1+1])/(tmp_contour_pts[2*j1+1]-tmp_contour_pts[2*i1+1])
					*(tmp_contour_pts[2*j1]-tmp_contour_pts[2*i1])); 
				}
    
				j1=i1; 
			}

			//  Sort the nodes, via a simple ?Bubble? sort.
			i1=0;
			while (i1<nodes-1) 
			{
				if (nodeX[i1]>nodeX[i1+1]) 
				{
					swap=nodeX[i1]; 
					nodeX[i1]=nodeX[i1+1]; 
					nodeX[i1+1]=swap; 
					if (i1) i1--; 
				}
				else 
				{
					i1++; 
				}
			}

			//  Fill the pixels between node pairs.
			for (i1=0; i1<nodes; i1+=2) 
			{
				if (nodeX[i1]>=IMAGE_RIGHT) 
					break;
			
				if (nodeX[i1+1]> IMAGE_LEFT ) 
				{
					if (nodeX[i1] < IMAGE_LEFT ) 
						nodeX[i1]=IMAGE_LEFT ;
      
					if (nodeX[i1+1]> IMAGE_RIGHT) 
						nodeX[i1+1]=IMAGE_RIGHT;
      
					for (j1=nodeX[i1]; j1<nodeX[i1+1]; j1++) 
					{
						pixelIndex[0] = j1;
						pixelIndex[1] = y;

						convex_hull_img->SetPixel(pixelIndex, 255);
						ObjectLabelImage->SetPixel(pixelIndex, current_label);
					}
				}
			}
		}
	}
}

void radConeSegmentation::FillContourImage(DoublePointArray & pts, BinaryImageType::Pointer sub_img)
{
	if (!sub_img) return;

	//first transform contour_pts
	BinaryImageType::IndexType img_start = sub_img->GetLargestPossibleRegion().GetIndex();
	BinaryImageType::IndexType pixelIndex;
	DoublePointType point;
	vector<float> tmp_contour_pts;

	for (unsigned int i=0; i<pts.size(); i++)
	{
		sub_img->TransformPhysicalPointToIndex(pts[i], pixelIndex);
		for (unsigned int k=0; k<2; k++)
			tmp_contour_pts.push_back(pixelIndex[k]-img_start[k]);
	}

	int  nodes, i1, j1, swap;
	int polyCorners = tmp_contour_pts.size()/2;
	int nodeX[1500];
	int IMAGE_RIGHT = sub_img->GetLargestPossibleRegion().GetSize()[0], IMAGE_LEFT = 0;
	
	//  Loop through the rows of the image.
	for (unsigned int y=0; y<sub_img->GetLargestPossibleRegion().GetSize()[1]; y++) 
	{
		//  Build a list of nodes.
		nodes=0; j1=polyCorners-1;
		for (i1=0; i1<polyCorners; i1++) 
		{
			if ((tmp_contour_pts[2*i1+1] < (double)y && tmp_contour_pts[2*j1+1] >= (double)y)
				|| (tmp_contour_pts[2*j1+1] < (double)y && tmp_contour_pts[2*i1+1] >= (double)y)) 
			{
				nodeX[nodes++]=(int)(tmp_contour_pts[2*i1]+(y-tmp_contour_pts[2*i1+1])/(tmp_contour_pts[2*j1+1]-tmp_contour_pts[2*i1+1])
				*(tmp_contour_pts[2*j1]-tmp_contour_pts[2*i1])); 
			}
    
			j1=i1; 
		}

		//  Sort the nodes, via a simple ?Bubble? sort.
		i1=0;
		while (i1<nodes-1) 
		{
			if (nodeX[i1]>nodeX[i1+1]) 
			{
				swap=nodeX[i1]; 
				nodeX[i1]=nodeX[i1+1]; 
				nodeX[i1+1]=swap; 
				if (i1) i1--; 
			}
			else 
			{
				i1++; 
			}
		}

		//  Fill the pixels between node pairs.
		for (i1=0; i1<nodes; i1+=2) 
		{
			if (nodeX[i1]>=IMAGE_RIGHT) 
				break;
			
			if (nodeX[i1+1]> IMAGE_LEFT ) 
			{
				if (nodeX[i1] < IMAGE_LEFT ) 
					nodeX[i1]=IMAGE_LEFT ;
      
				if (nodeX[i1+1]> IMAGE_RIGHT) 
					nodeX[i1+1]=IMAGE_RIGHT;
      
				for (j1=nodeX[i1]; j1<nodeX[i1+1]; j1++) 
				{
					pixelIndex[0] = j1+img_start[0];
					pixelIndex[1] = y+img_start[1];

					sub_img->SetPixel(pixelIndex, 255);
				}
			}
		}
	}
}

void radConeSegmentation::ExpandConvexHull(BinaryImageType::Pointer sub_cone_img, 
										   itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer sub_cone_list)
{
	unsigned int max_id, max_num = 0;
	for(unsigned int i = 0; i < sub_cone_list->GetNumberOfLabelObjects(); i++)
    {
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= sub_cone_list->GetNthLabelObject(i);

		if (labelObject->GetNumberOfPixels() > max_num)
		{
			max_num = labelObject->GetNumberOfPixels();
			max_id = i;
		}
    }

	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= sub_cone_list->GetNthLabelObject(max_id);
	int radius = (labelObject->GetEquivalentSphericalRadius()+labelObject->GetEquivalentEllipsoidDiameter()[1]/2)/2;
	sub_cone_img->FillBuffer(0);
	
	BinaryIteratorType sub_cone_it(sub_cone_img, sub_cone_img->GetLargestPossibleRegion());
	for (sub_cone_it.GoToBegin(); !sub_cone_it.IsAtEnd(); ++sub_cone_it)
	{
		if ((sub_cone_it.GetIndex()[0]-labelObject->GetCentroid()[0])*(sub_cone_it.GetIndex()[0]-labelObject->GetCentroid()[0])
			+(sub_cone_it.GetIndex()[1]-labelObject->GetCentroid()[1])*(sub_cone_it.GetIndex()[1]-labelObject->GetCentroid()[1])
			<= radius*radius)
		{
			sub_cone_it.Set(255);
		}
	}
}

float radConeSegmentation::BilinearInterpolate(FloatImageType::Pointer img, double x, double y)
{
	int xx, yy, m, n;
	xx = x;
	yy = y;
	double dx, dy, s;
	dx = std::max<double>(std::min<double>(x-xx, 1), 0);
	dy = std::max<double>(std::min<double>(y-yy, 1), 0);
	FloatImageType::SizeType img_size = img->GetLargestPossibleRegion().GetSize();
	FloatImageType::IndexType pixelIndex;

	float result = 0;
	for (m=0; m<=1; m++)
		for (n=0; n<=1; n++)
		{
			pixelIndex[0] = EnforceRange(xx+m, img_size[0]);
			pixelIndex[1] = EnforceRange(yy+n, img_size[1]);
			s = fabs(1 - m - dx) * fabs(1 - n - dy);
			result += img->GetPixel(pixelIndex) * s;
		}
	return result;
}

float radConeSegmentation::BicubicInterpolate(FloatImageType::Pointer img, FloatImageType::Pointer img_dx, 
											  FloatImageType::Pointer img_dy, FloatImageType::Pointer img_dxy, 
											  double x, double y)
{
	int x0 = x;
	int y0 = y;
	int x1 = x0+1;
	int y1 = y0+1;

	FloatImageType::SizeType img_size = img->GetLargestPossibleRegion().GetSize();

	x0 = std::min<int>(std::max<int>(x0, 0), img_size[0]-1);
	x1 = std::min<int>(std::max<int>(x1, 0), img_size[0]-1);
	y0 = std::min<int>(std::max<int>(y0, 0), img_size[1]-1);
	y1 = std::min<int>(std::max<int>(y1, 0), img_size[1]-1);

	double dx = x - x0;
	double dy = y- y0;
	double dx2 = dx*dx;
	double dy2 = dy*dy;
	double dx3 = dx*dx2;
	double dy3 = dy*dy2;

	FloatImageType::IndexType pixelIndex[4];
	pixelIndex[0][0] = x0;
	pixelIndex[0][1] = y0;  //offsets[0][0]
	pixelIndex[1][0] = x1;
	pixelIndex[1][1] = y0;  //offsets[1][0]
	pixelIndex[2][0] = x0;
	pixelIndex[2][1] = y1;  //offsets[0][1]
	pixelIndex[3][0] = x1;
	pixelIndex[3][1] = y1;  //offsets[1][1]

	double a[4][4];
	a[0][0] = img->GetPixel(pixelIndex[0]);
	a[1][0] = img_dx->GetPixel(pixelIndex[0]);
	a[2][0] = -3*img->GetPixel(pixelIndex[0]) + 3*img->GetPixel(pixelIndex[1]) -2*img_dx->GetPixel(pixelIndex[0]) - img_dx->GetPixel(pixelIndex[1]);
	a[3][0] =   2*img->GetPixel(pixelIndex[0]) -  2*img->GetPixel(pixelIndex[1]) +   img_dx->GetPixel(pixelIndex[0]) +img_dx->GetPixel(pixelIndex[1]);

	a[0][1] = img_dy->GetPixel(pixelIndex[0]);
	a[1][1] = img_dxy->GetPixel(pixelIndex[0]);
	a[2][1] = -3*img_dy->GetPixel(pixelIndex[0]) + 3*img_dy->GetPixel(pixelIndex[1]) - 2*img_dxy->GetPixel(pixelIndex[0]) - img_dxy->GetPixel(pixelIndex[1]);
	a[3][1] = 2*img_dy->GetPixel(pixelIndex[0]) - 2*img_dy->GetPixel(pixelIndex[1]) + img_dxy->GetPixel(pixelIndex[0]) + img_dxy->GetPixel(pixelIndex[1]);

	a[0][2] =      -3*img->GetPixel(pixelIndex[0])      + 3*img->GetPixel(pixelIndex[2])       -2*img_dy->GetPixel(pixelIndex[0])        - img_dy->GetPixel(pixelIndex[2]);
	a[1][2] = -3*img_dx->GetPixel(pixelIndex[0]) + 3*img_dx->GetPixel(pixelIndex[2]) -2*img_dxy->GetPixel(pixelIndex[0]) - img_dxy->GetPixel(pixelIndex[2]);
	a[2][2] =		     9*img->GetPixel(pixelIndex[0])      -        9*img->GetPixel(pixelIndex[1])     -        9*img->GetPixel(pixelIndex[2])     +    9*img->GetPixel(pixelIndex[3]) +
								6*img_dx->GetPixel(pixelIndex[0])   +    3*img_dx->GetPixel(pixelIndex[1])   -     6*img_dx->GetPixel(pixelIndex[2]) -    3*img_dx->GetPixel(pixelIndex[3]) +
								6*img_dy->GetPixel(pixelIndex[0])   -     6*img_dy->GetPixel(pixelIndex[1]) +      3*img_dy->GetPixel(pixelIndex[2]) -    3*img_dy->GetPixel(pixelIndex[3]) +
							4*img_dxy->GetPixel(pixelIndex[0]) + 2*img_dxy->GetPixel(pixelIndex[1]) + 2*img_dxy->GetPixel(pixelIndex[2]) + img_dxy->GetPixel(pixelIndex[3]);
	a[3][2] =		    -6*img->GetPixel(pixelIndex[0])      +      6*img->GetPixel(pixelIndex[1])     +       6*img->GetPixel(pixelIndex[2])     -     6*img->GetPixel(pixelIndex[3]) +
							(-3)*img_dx->GetPixel(pixelIndex[0])   -     3*img_dx->GetPixel(pixelIndex[1])   +    3*img_dx->GetPixel(pixelIndex[2]) +   3*img_dx->GetPixel(pixelIndex[3]) +
							(-4)*img_dy->GetPixel(pixelIndex[0])   +    4*img_dy->GetPixel(pixelIndex[1])    -    2*img_dy->GetPixel(pixelIndex[2]) +   2*img_dy->GetPixel(pixelIndex[3]) +
						(-2)*img_dxy->GetPixel(pixelIndex[0])  - 2*img_dxy->GetPixel(pixelIndex[1])   -    img_dxy->GetPixel(pixelIndex[2])   -  img_dxy->GetPixel(pixelIndex[3]);

	a[0][3] =      2*img->GetPixel(pixelIndex[0])        - 2*img->GetPixel(pixelIndex[2])       + img_dy->GetPixel(pixelIndex[0])        + img_dy->GetPixel(pixelIndex[2]);
	a[1][3] = 2*img_dx->GetPixel(pixelIndex[0])  - 2*img_dx->GetPixel(pixelIndex[2])  + img_dxy->GetPixel(pixelIndex[0]) + img_dxy->GetPixel(pixelIndex[2]);
	a[2][3] =		    -6*img->GetPixel(pixelIndex[0])      +      6*img->GetPixel(pixelIndex[1])     +       6*img->GetPixel(pixelIndex[2])     -     6*img->GetPixel(pixelIndex[3]) +
							(-4)*img_dx->GetPixel(pixelIndex[0])   -     2*img_dx->GetPixel(pixelIndex[1])   +    4*img_dx->GetPixel(pixelIndex[2]) +   2*img_dx->GetPixel(pixelIndex[3]) +
							(-3)*img_dy->GetPixel(pixelIndex[0])   +    3*img_dy->GetPixel(pixelIndex[1])    -    3*img_dy->GetPixel(pixelIndex[2]) +   3*img_dy->GetPixel(pixelIndex[3]) +
						(-2)*img_dxy->GetPixel(pixelIndex[0])  -     img_dxy->GetPixel(pixelIndex[1]) -  2*img_dxy->GetPixel(pixelIndex[2])   -  img_dxy->GetPixel(pixelIndex[3]);
	a[3][3] =		     4*img->GetPixel(pixelIndex[0])      -        4*img->GetPixel(pixelIndex[1])     -        4*img->GetPixel(pixelIndex[2])     +    4*img->GetPixel(pixelIndex[3]) +
								2*img_dx->GetPixel(pixelIndex[0])   +    2*img_dx->GetPixel(pixelIndex[1])   -     2*img_dx->GetPixel(pixelIndex[2]) -    2*img_dx->GetPixel(pixelIndex[3]) +
								2*img_dy->GetPixel(pixelIndex[0])   -     2*img_dy->GetPixel(pixelIndex[1]) +      2*img_dy->GetPixel(pixelIndex[2]) -    2*img_dy->GetPixel(pixelIndex[3]) +
								img_dxy->GetPixel(pixelIndex[0]) +     img_dxy->GetPixel(pixelIndex[1]) +      img_dxy->GetPixel(pixelIndex[2]) + img_dxy->GetPixel(pixelIndex[3]);

	return (a[0][0] + a[0][1]*dy + a[0][2]*dy2 + a[0][3]*dy3 
		+ a[1][0]*dx +   a[1][1]*dx*dy   + a[1][2]*dx*dy2   + a[1][3]*dx*dy3 
		+ a[2][0]*dx2 + a[2][1]*dx2*dy + a[2][2]*dx2*dy2 + a[2][3]*dx2*dy3
		+ a[3][0]*dx3 + a[3][1]*dx3*dy + a[3][2]*dx3*dy2 + a[3][3]*dx3*dy3);
}

bool cmp(const pair<DoublePointType, double> d1, const pair<DoublePointType, double> d2)
{
	return d1.second < d2.second;
}

bool radConeSegmentation::SearchForSubPixelLocation(BinaryImageType::IndexType & region_center, 
													BinaryImageType::IndexType & boundary_pixel, 
													FloatImageType::Pointer laplacian_img, 
													FloatImageType::Pointer laplacian_dx, 
													FloatImageType::Pointer laplacian_dy,
													FloatImageType::Pointer laplacian_dxy,
													int start_id, int end_id, double step, 
													DoublePointType & result_pt)
{
	if (start_id > end_id) return false;

	DoublePointType direction;
	double magn = sqrt((double)(region_center[0]-boundary_pixel[0])*(region_center[0]-boundary_pixel[0])
		+(double)(region_center[1]-boundary_pixel[1])*(region_center[1]-boundary_pixel[1]));
	for (unsigned int i=0; i<2; i++)
		direction[i] = (boundary_pixel[i]-region_center[i])/magn;

	vector< pair<DoublePointType, double> > zero_crossing_pts;
	for (int id = start_id; id<=end_id; id++)
	{
		pair<DoublePointType, double> tmp_pair;
		for (unsigned int i=0; i<2; i++)
			tmp_pair.first[i] = id*direction[i]*step + boundary_pixel[i];

		//tmp_pair.second = fabs(BilinearInterpolate(laplacian_img, tmp_pair.first[0], tmp_pair.first[1]));
		tmp_pair.second = fabs(BicubicInterpolate(laplacian_img, laplacian_dx, laplacian_dy, laplacian_dxy, tmp_pair.first[0], tmp_pair.first[1]));
		zero_crossing_pts.push_back(tmp_pair);
	}

	//search for local minimum
	vector< pair<DoublePointType, double> >::iterator min_it = std::min_element(zero_crossing_pts.begin(), 
		zero_crossing_pts.end(), cmp);

	result_pt = min_it->first;
	return true;
}

void radConeSegmentation::ExtractBoundaryPixels(BinaryImageType::Pointer seg_img, DoublePointArray & contour_pts, 
												int foreground)
{
	contour_pts.clear();
	FloatImageType::SizeType radius;
	radius.Fill(1);
	
	vector< pair<BinaryImageType::IndexType, bool> > tmp_boundary_index_list;
	itk::NeighborhoodIterator<BinaryImageType> iterator(radius, seg_img, seg_img->GetLargestPossibleRegion());
	while(!iterator.IsAtEnd())
    {
		// Set the current pixel to white
		if (iterator.GetCenterPixel() == foreground) //object pixel
		{
			for(unsigned int i = 0; i < 9; i++)
			{
				if (i == 0 || i == 2 || i == 4 || i == 6 || i == 9)
					continue;

				// BinaryImageType::IndexType index = iterator.GetIndex(i);
				
				bool IsInBounds;
				int val = iterator.GetPixel(i, IsInBounds);
				if(!IsInBounds || val != foreground)
				{
					tmp_boundary_index_list.push_back(std::make_pair(iterator.GetIndex(), false));
					break;
				}
			}
		}
		++iterator;
    }

	//try to connect all boundary pixels within a loop
	if (tmp_boundary_index_list.empty())
		return;

	unsigned int number_of_boundary_pixels = tmp_boundary_index_list.size();
	BinaryImageType::IndexType start_index = tmp_boundary_index_list[0].first;
	unsigned int id = 0;
	unsigned int current_dist, min_dist;
	while (contour_pts.size() != number_of_boundary_pixels)
	{
		//select a seed point
		DoublePointType tmp_pt;
		tmp_pt[0] = start_index[0];
		tmp_pt[1] = start_index[1];

		contour_pts.push_back(tmp_pt);
		tmp_boundary_index_list[id].second = true;

		min_dist = std::numeric_limits<int>::max();
		id = number_of_boundary_pixels;
		//search for the closest point to the seed point
		for (unsigned int i=0; i<tmp_boundary_index_list.size(); i++)
		{
			if (tmp_boundary_index_list[i].second)
				continue;

			current_dist = (start_index[0]-tmp_boundary_index_list[i].first[0])
				* (start_index[0]-tmp_boundary_index_list[i].first[0])
				+ (start_index[1]-tmp_boundary_index_list[i].first[1])
				* (start_index[1]-tmp_boundary_index_list[i].first[1]);

			if (current_dist < min_dist && current_dist <= 2)
			{
				min_dist = current_dist;
				id = i;
			}
		}

		if (id == number_of_boundary_pixels)
			break;
		else
			start_index = tmp_boundary_index_list[id].first;
	}

	if (contour_pts.size() != number_of_boundary_pixels)
	{
		unsigned int current_dist1;
		for (unsigned int i=0; i<tmp_boundary_index_list.size(); i++)
		{
			if (!tmp_boundary_index_list[i].second)
			{
				for (unsigned int j=0; j<contour_pts.size(); j++)
				{
					if (j == 0)
					{
						current_dist = (contour_pts[contour_pts.size()-1][0]-tmp_boundary_index_list[i].first[0])
							* (contour_pts[contour_pts.size()-1][0]-tmp_boundary_index_list[i].first[0])
							+ (contour_pts[contour_pts.size()-1][1]-tmp_boundary_index_list[i].first[1])
							* (contour_pts[contour_pts.size()-1][1]-tmp_boundary_index_list[i].first[1]);
					}
					else
					{
						current_dist = (contour_pts[j-1][0]-tmp_boundary_index_list[i].first[0])
							* (contour_pts[j-1][0]-tmp_boundary_index_list[i].first[0])
							+ (contour_pts[j-1][1]-tmp_boundary_index_list[i].first[1])
							* (contour_pts[j-1][1]-tmp_boundary_index_list[i].first[1]);
					}

					current_dist1 = (contour_pts[j][0]-tmp_boundary_index_list[i].first[0])
							* (contour_pts[j][0]-tmp_boundary_index_list[i].first[0])
							+ (contour_pts[j][1]-tmp_boundary_index_list[i].first[1])
							* (contour_pts[j][1]-tmp_boundary_index_list[i].first[1]);

					if (current_dist <=2 && current_dist1 <= 2)
					{
						DoublePointType tmp_pt;
						tmp_pt[0] = tmp_boundary_index_list[i].first[0];
						tmp_pt[1] = tmp_boundary_index_list[i].first[1];

						contour_pts.insert(contour_pts.begin()+j, tmp_pt);
						break;
					}
				}
			}
		}
	}
}

vnl_matrix<double> radConeSegmentation::ComputeP(double alpha, double beta, double gamma, double N) throw (int)
{
	double a = gamma*(2*alpha+6*beta)+1;
	double b = gamma*(-alpha-4*beta);
	double c = gamma*beta;
 
	vnl_matrix<double> P(N,N);
	P.fill(0);
 
	//fill diagonal
	P.fill_diagonal(a);
 
	//fill next two diagonals
	for (int i=0; i<(N-1); i++)
    {
		P(i+1,i) = b;
		P(i,i+1) = b;
    }
  
	//Moreover
	P(0, N-1)=b;
	P(N-1, 0)=b;
 
	//fill next two diagonals
	for (int i=0; i<(N-2); i++)
    {
		P(i+2,i) = c;
		P(i,i+2) = c;
    }
  
	//Moreover
	P(0, N-2)=c;
	P(1, N-1)=c;
	P(N-2, 0)=c;
	P(N-1, 1)=c;
 
	if ( vnl_determinant(P) == 0.0 )
    {
		std::cerr << "Singular matrix. Determinant is 0." << std::endl;
		throw 2;
    }
 
	//Compute the inverse of the matrix P
	vnl_matrix< double > Pinv;
	Pinv = vnl_matrix_inverse< double >(P);

	return Pinv.transpose();
}

vnl_vector<double> radConeSegmentation::SampleImage(vnl_vector<double> x, vnl_vector<double> y, 
													VectorImageType::Pointer gradient, int position)
{
	int size;
	size = x.size();
	vnl_vector<double> ans(size);
 
	VectorImageType::IndexType index;
	for (int i=0; i<size; i++)
    {
		index[0] = x[i];
		index[1] = y[i];

		if( gradient->GetLargestPossibleRegion().IsInside(index) )
		{
			ans[i] = gradient->GetPixel(index)[position];
		}
		else
		{
			ans.clear();
			return ans;
		}
    }
	return ans;
}


bool radConeSegmentation::SnakeRefinement(FloatImageType::Pointer sub_img, DoublePointArray &initinal_contour)
{
	/*double alpha = 0.05;
	double beta = 2.0;
	double gamma = 0.001;
	double sigma = 1;
	double iterations = 100;*/

	double alpha = 0.01;
	double beta = 10;
	double gamma = 0.002;
	double sigma = 2;
	double iterations = 300;
	
	double N;
	int n = initinal_contour.size();
	vnl_matrix<double> P;
	vnl_vector<double> v(2*(n+1));

	for (int i=0; i<n; i++)
    {
		v[2*i] = initinal_contour[i][0];
		v[2*i+1] = initinal_contour[i][1];
    }
	v[2*n]=v[0];
	v[2*n+1]=v[1];

	//Computes P matrix.
	N = v.size()/2;
	P = ComputeP(alpha, beta, gamma, N);

	typedef itk::GradientMagnitudeImageFilter<FloatImageType, FloatImageType >   GradMagfilterType;
	//Computes the magnitude gradient
	GradMagfilterType::Pointer gradientMagnitudeFilter = GradMagfilterType::New();
	gradientMagnitudeFilter->SetInput( sub_img );
	gradientMagnitudeFilter->Update();

	typedef itk::GradientRecursiveGaussianImageFilter<FloatImageType, VectorImageType>  FilterType;
	FilterType::Pointer gradientFilter = FilterType::New();
	gradientFilter->SetInput( gradientMagnitudeFilter->GetOutput() );
	gradientFilter->SetSigma( sigma );
	gradientFilter->Update();

	//Loop
	vnl_vector<double> x(N);
	vnl_vector<double> y(N);
 
	//std::cout << "Initial snake" << std::endl;
	for (int i = 0; i < N; i++)
    {
		x[i] = v[2*i];
		y[i] = v[2*i+1];
		//std::cout << "(" << x[i] << ", " << y[i] << ")" << std::endl;
	}

	for (int i = 0; i < iterations; i++)
    {
		vnl_vector<double> fex;
		vnl_vector<double> fey;
		fex = SampleImage(x, y, gradientFilter->GetOutput(), 0);
		fey = SampleImage(x, y, gradientFilter->GetOutput(), 1);
 
		if (!fex.empty() && !fey.empty())
		{
			x = (x+gamma*fex).post_multiply(P);
			y = (y+gamma*fey).post_multiply(P);
		}
		else
		{
			//QMessageBox msgBox;
			//msgBox.setText("A contour is outside the gradient field!");
			//msgBox.exec();
			return false;
		}
    }

	//std::cout << "Final snake after " << iterations << " iterations" << std::endl;
	initinal_contour.clear();
	initinal_contour.resize(N);
	for (int i=0; i<N; i++)
    {
		initinal_contour[i][0] = x[i];
		initinal_contour[i][1] = y[i];
		//std::cout << "(" << x[i] << ", " << y[i] << ")" << std::endl;
    }
	return true;
}

bool radConeSegmentation::IsCurrentIndexInRightPosition(BinaryImageType::IndexType index, vector< BinaryImageType::IndexType > &indice_array)
{
	if (indice_array.empty())
		return true;
	else
	{
		int last_pos = indice_array.size()-1;
		BinaryImageType::IndexType last_index = indice_array[last_pos];
		if (abs(index[0]-last_index[0])+abs(index[1]-last_index[1]) < 5)
			return true;
		else
			return false;
	}
}

void radConeSegmentation::SubPixelImprovement(FloatImageType::Pointer sub_img, FloatImageType::Pointer laplacian_img, 
											  BinaryImageType::Pointer sub_seg_img, DoublePointArray & cell_contour, 
											  FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy, 
											  FloatImageType::Pointer laplacian_dxy)
{
	FloatImageType::SizeType radius;
	radius.Fill(1);
	vector< BinaryImageType::IndexType > boundary_index_list;

	typedef itk::ShapedNeighborhoodIterator<BinaryImageType> ShapedNeighborhoodIteratorType;
	
	ShapedNeighborhoodIteratorType::OffsetType top = {{0,-1}};
	ShapedNeighborhoodIteratorType::OffsetType bottom = {{0,1}};
	ShapedNeighborhoodIteratorType::OffsetType left = {{-1,0}};
	ShapedNeighborhoodIteratorType::OffsetType right = {{1,0}};

	typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< BinaryImageType > FaceCalculatorType;
	FaceCalculatorType faceCalculator;
	FaceCalculatorType::FaceListType faceList;
	FaceCalculatorType::FaceListType::iterator faceListIterator;
 
	faceList = faceCalculator( sub_seg_img, sub_seg_img->GetRequestedRegion(), radius );

	bool boundary_flag;
	vector< pair<BinaryImageType::IndexType, bool> > tmp_boundary_index_list;
	for ( faceListIterator=faceList.begin(); faceListIterator != faceList.end(); ++faceListIterator)
	{
		ShapedNeighborhoodIteratorType it( radius, sub_seg_img, *faceListIterator );
		it.ActivateOffset(top);
		it.ActivateOffset(bottom);
		it.ActivateOffset(left);
		it.ActivateOffset(right);

		//main region
		if (faceListIterator->GetIndex()[0] - sub_seg_img->GetLargestPossibleRegion().GetIndex()[0] == 1
			&& faceListIterator->GetIndex()[1] - sub_seg_img->GetLargestPossibleRegion().GetIndex()[1] == 1)
		{
			for (it.GoToBegin(); !it.IsAtEnd(); ++it)
			{	
				if (sub_seg_img->GetPixel(it.GetIndex()) != 0)
				{
					boundary_flag = false;
					ShapedNeighborhoodIteratorType::ConstIterator ci;
					for (ci = it.Begin(); ci != it.End(); ci++)
					{
						if ((int)ci.Get() == 0)
						{
							boundary_flag = true;
							break;
						}
					}

					if (boundary_flag)
					{
						tmp_boundary_index_list.push_back(std::make_pair(it.GetIndex(), false));
					}
				}
			}
		}
		else
		{
			for (it.GoToBegin(); !it.IsAtEnd(); ++it)
			{
				if (sub_seg_img->GetPixel(it.GetIndex()) != 0)
				{
					tmp_boundary_index_list.push_back(std::make_pair(it.GetIndex(), false));
				}
			}
		}
	}

	//try to connect all boundary pixels within a loop
	if (tmp_boundary_index_list.empty())
		return;

	unsigned int number_of_boundary_pixels = tmp_boundary_index_list.size();
	BinaryImageType::IndexType start_index = tmp_boundary_index_list[0].first;
	unsigned int id = 0;
	unsigned int current_dist, min_dist;
	while (boundary_index_list.size() != number_of_boundary_pixels)
	{
		//select a seed point
		if (IsCurrentIndexInRightPosition(start_index, boundary_index_list))
		{
			boundary_index_list.push_back(start_index);
			tmp_boundary_index_list[id].second = true;
		}
		else
			break;

		min_dist = std::numeric_limits<int>::max();
		id = number_of_boundary_pixels;
		//search for the closest point to the seed point
		for (unsigned int i=0; i<tmp_boundary_index_list.size(); i++)
		{
			if (tmp_boundary_index_list[i].second)
				continue;

			current_dist = (start_index[0]-tmp_boundary_index_list[i].first[0])
				* (start_index[0]-tmp_boundary_index_list[i].first[0])
				+ (start_index[1]-tmp_boundary_index_list[i].first[1])
				* (start_index[1]-tmp_boundary_index_list[i].first[1]);

			if (current_dist < min_dist)
			{
				min_dist = current_dist;
				id = i;
			}
		}

		if (id == number_of_boundary_pixels)
			break;
		else
		{
			start_index = tmp_boundary_index_list[id].first;
		}
	}

	if (boundary_index_list.empty())
		return;

	//next for each boundary point, calculate sub-pixel in terms of zero-crossing
	//create an status image
	BinaryImageType::IndexType region_center;
	region_center.Fill(0);

	for (unsigned int i=0; i<boundary_index_list.size(); i++)
	{
		for (unsigned int j=0; j<2; j++)
			region_center[j] += boundary_index_list[i][j];
	}

	for (unsigned int j=0; j<2; j++)
		region_center[j] /= boundary_index_list.size();

	DoublePointType tmp_pt;
	unsigned int sampling_rate = 1;
	for (unsigned int i=0; i<boundary_index_list.size(); i=i+sampling_rate)
	{
		//temporally remove 
		if (i >= boundary_index_list.size()) break;
		if (boundary_index_list[i][0] == 0 || boundary_index_list[i][0] == laplacian_img->GetLargestPossibleRegion().GetSize()[0]-1
			|| boundary_index_list[i][1] == 0 || boundary_index_list[i][1] == laplacian_img->GetLargestPossibleRegion().GetSize()[1]-1)
		{
			for (unsigned int j=0; j<2; j++)
				tmp_pt[j] = boundary_index_list[i][j];
		}
		else
			SearchForSubPixelLocation(region_center, boundary_index_list[i], laplacian_img, laplacian_dx, laplacian_dy, 
				laplacian_dxy, -8, 8, 0.1, tmp_pt);

		cell_contour.push_back(tmp_pt);
	}
	
	//last step use spline snake
	if (cell_contour.size() < 3)
		return;

	SnakeRefinement(sub_img, cell_contour);
}

bool radConeSegmentation::IsSegmentationValid(BinaryImageType::Pointer tmp_seg, BinaryImageType::Pointer seg)
{
	BinaryIteratorType tmp_seg_it(tmp_seg, tmp_seg->GetLargestPossibleRegion());
	BinaryIteratorType seg_it(seg, tmp_seg->GetLargestPossibleRegion());

	float num_of_overlap_seg = 0, num_of_tmp_seg = 0;
	for (tmp_seg_it.GoToBegin(), seg_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it, ++tmp_seg_it)
	{
		if (tmp_seg_it.Get() != 0)
		{
			num_of_tmp_seg += 1.0;
			if (seg_it.Get() != 0)
				num_of_overlap_seg += 1.0;
		}
	}

	if (num_of_tmp_seg == 0) return false;
	else if (num_of_overlap_seg/num_of_tmp_seg > 0.5) return false;
	else 
		return true;

}

bool radConeSegmentation::IsSegmentationValid(BinaryImageType::Pointer img, double radius_threshold, 
											  double roundness_threshold)
{
	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer label_img = LabelObjectSegmentation(img);

	if (label_img->GetNumberOfLabelObjects() != 1)
		return false;

	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= label_img->GetNthLabelObject(0);
	
	if (labelObject->GetEquivalentSphericalRadius() < radius_threshold
		|| labelObject->GetRoundness() < roundness_threshold)
		return false;

	return true;
}

bool radConeSegmentation::IsSegmentationValid(BinaryImageType::Pointer img1, 
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType *img2_label)
{
	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer label_img = LabelObjectSegmentation(img1);
	if (label_img->GetNumberOfLabelObjects() != 1)
		return false;

	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* img1_label 
			= label_img->GetNthLabelObject(0);

	if (img1_label->GetEquivalentSphericalRadius() >= BlobSizeThreshold
		|| img2_label->GetEquivalentSphericalRadius() >= BlobSizeThreshold)
		return true;
	else
		return false;
}

void radConeSegmentation::AssignLabelImage(BinaryImageType::Pointer sub_seg_img, int label, 
										   ShortImageType::Pointer label_img)
{
	if (!label_img || label == -1)
		return;

	//label images for the segmenation from region pairs, for checking segmenation touching in the remaining steps
	BinaryIteratorType seg_it(sub_seg_img, sub_seg_img->GetLargestPossibleRegion());
	ShortIteratorType label_it(label_img, sub_seg_img->GetLargestPossibleRegion());
				
	for (seg_it.GoToBegin(), label_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it, ++label_it) 
	{
		if (seg_it.Get() != 0)
			label_it.Set(label);
	}
}

void radConeSegmentation::ReassignBinaryImageFromLabelImage(ShortImageType::Pointer label_img, int label, 
															BinaryImageType::Pointer res_img)
{
	itk::ImageRegionConstIterator<ShortImageType> label_it(label_img, label_img->GetLargestPossibleRegion());
	itk::ImageRegionIterator<BinaryImageType> res_it(res_img, res_img->GetLargestPossibleRegion());

	for (label_it.GoToBegin(), res_it.GoToBegin(); !label_it.IsAtEnd(); ++label_it, ++res_it)
	{
		if (label_it.Get() == label)
			res_it.Set(255);
	}
}

bool radConeSegmentation::ImproveConeSegmentation(FloatImageType::Pointer sub_img, BinaryImageType::Pointer sub_cone_img, 
												  FloatImageType::Pointer laplacian_img, DoublePointArray & cell_contour,
												  FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
												  FloatImageType::Pointer laplacian_dxy, ShortImageType::Pointer label_img,
												  int label, bool determination_flag)
{
	if (sub_img->GetLargestPossibleRegion().GetSize()[0] < 8 
		|| sub_img->GetLargestPossibleRegion().GetSize()[1] < 8 )
		return false;

	typedef  itk::CurvatureAnisotropicDiffusionImageFilter< FloatImageType, FloatImageType > SmoothingFilterType;
	SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
	smoothing->SetTimeStep( 0.125 );
	smoothing->SetNumberOfIterations( 5 );
	smoothing->SetConductanceParameter( 9.0 );
	smoothing->SetInput( sub_img );

	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< FloatImageType, FloatImageType >  GradientFilterType;
	typedef itk::SigmoidImageFilter< FloatImageType, FloatImageType >  SigmoidFilterType;
	typedef itk::BinaryThresholdImageFilter< FloatImageType, BinaryImageType > ThresholdingFilterType;

	GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
	SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();

	gradientMagnitude->SetInput( smoothing->GetOutput() );
	sigmoid->SetInput( gradientMagnitude->GetOutput() );
	gradientMagnitude->SetSigma(  m_Scale  );

	gradientMagnitude->Update();
	FloatImageType::Pointer gradient_magn = gradientMagnitude->GetOutput();
	float min_magn, max_magn;
	ComputeIntensityRange<FloatImageType>(gradient_magn, min_magn, max_magn);

	sigmoid->SetOutputMinimum(  0.0  );
	sigmoid->SetOutputMaximum(  1.0  );
	const double alpha =  (min_magn-max_magn)/9;
	const double beta  =  (min_magn+max_magn)/6;
	// Software Guide : BeginCodeSnippet
	sigmoid->SetAlpha( alpha );
	sigmoid->SetBeta(  beta  );

	//perform symmetry
	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer sub_cone_labeling = LabelObjectSegmentation(sub_cone_img);
	ExpandConvexHull(sub_cone_img, sub_cone_labeling);


	FloatImageType::Pointer dist_img = DistanceTransform<BinaryImageType, FloatImageType>(sub_cone_img, false);
	FloatImageType::Pointer abs_img = AbsoluteImage<FloatImageType, FloatImageType>(dist_img);
	FloatImageType::PixelType min_val, max_val;
	ComputeIntensityRange<FloatImageType>(abs_img, min_val, max_val);
	FloatImageType::Pointer shape_img = MultipleImageWithConstant<FloatImageType, FloatImageType>(abs_img, 1.0/max_val);

	sigmoid->Update();
	FloatImageType::Pointer sigmoid_img = sigmoid->GetOutput();

	float thresh = 0.9;
	FloatImageType::Pointer speed_img = AllocateImage<FloatImageType, FloatImageType>(sigmoid_img);
	FloatIteratorType speed_it(speed_img, speed_img->GetLargestPossibleRegion());
	FloatIteratorType sigmoid_it(sigmoid_img, sigmoid_img->GetLargestPossibleRegion());
	FloatIteratorType shape_it(shape_img, shape_img->GetLargestPossibleRegion());
	for (speed_it.GoToBegin(), sigmoid_it.GoToBegin(), shape_it.GoToBegin(); !speed_it.IsAtEnd(); 
		++speed_it, ++shape_it, ++sigmoid_it)
	{
		if (shape_it.Get() < 0.3)
			speed_it.Set(sigmoid_it.Get()*(1-thresh));
			//speed_it.Set(0);
		else
			speed_it.Set(thresh+sigmoid_it.Get()*(1-thresh));
			//speed_it.Set(1);
	}
	
	double propagationScaling = 2.0;
	typedef  itk::GeodesicActiveContourLevelSetImageFilter< FloatImageType, FloatImageType >  GeodesicActiveContourFilterType;
	GeodesicActiveContourFilterType::Pointer geodesicActiveContour = GeodesicActiveContourFilterType::New();
	geodesicActiveContour->SetPropagationScaling( propagationScaling );
	geodesicActiveContour->SetCurvatureScaling( 1.0 );
	geodesicActiveContour->SetAdvectionScaling( 1.0 );
	geodesicActiveContour->SetMaximumRMSError( 0.005 );
	geodesicActiveContour->SetInput( dist_img );
	geodesicActiveContour->SetFeatureImage( speed_img );

	geodesicActiveContour->SetNumberOfIterations( SingleGACIterationNumber );
	geodesicActiveContour->Update();
	
	thresholder->SetInput( geodesicActiveContour->GetOutput() );
	thresholder->SetLowerThreshold( -1000.0 );
	thresholder->SetUpperThreshold( 0.0  );
	thresholder->SetOutsideValue(  itk::NumericTraits< BinaryImageType::PixelType >::min()  );
	thresholder->SetInsideValue(  itk::NumericTraits< BinaryImageType::PixelType >::max() );
	thresholder->Update();

	if (!IsSegmentationValid(thresholder->GetOutput(), FinalSegmentedImage))
	{
		AssignLabelImage(thresholder->GetOutput(), label, label_img);
		return false;
	}

	//here expand threshold image to sub-pixel level segmentation & regenerate binary segmentation
	SubPixelImprovement(sub_img, laplacian_img, thresholder->GetOutput(), cell_contour, laplacian_dx, laplacian_dy, laplacian_dxy);
	BinaryImageType::Pointer sub_tmp_img = AllocateImage<BinaryImageType, BinaryImageType>(thresholder->GetOutput());
	sub_tmp_img->FillBuffer(0);
	FillContourImage(cell_contour, sub_tmp_img);
	
	BinaryIteratorType sub_seg_it(sub_tmp_img, sub_tmp_img->GetLargestPossibleRegion());
	BinaryIteratorType seg_it(FinalSegmentedImage, sub_tmp_img->GetLargestPossibleRegion());

	if (determination_flag)
	{
		if (!IsSegmentationValid(sub_tmp_img, ConeRadiusThreshold, ConeRoundnessThreshold))
		{
			AssignLabelImage(sub_tmp_img, label, label_img);
			return false;
		}

		BinaryIteratorType sub_cone_it(sub_cone_img, sub_cone_img->GetLargestPossibleRegion());
		double post_seg_num = 0, prior_seg_num = 0;
		for (sub_cone_it.GoToBegin(), sub_seg_it.GoToBegin(); !sub_cone_it.IsAtEnd(); ++sub_cone_it, ++sub_seg_it)
		{
			if (sub_cone_it.Get() != 0) 
				prior_seg_num += 1;
			if (sub_seg_it.Get() != 0)
				post_seg_num += 1;
		}

		if (post_seg_num/prior_seg_num > 1.6)
		{
			for (seg_it.GoToBegin(), sub_cone_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it, ++sub_cone_it)
			{
				if (sub_cone_it.Get() != 0)
				{
					seg_it.Set(255);
				}
			}

			AssignLabelImage(sub_cone_img, label, label_img);
			return true;
		}
	}

	for (seg_it.GoToBegin(), sub_seg_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it, ++sub_seg_it)
	{
		if (sub_seg_it.Get() != 0)
		{
			seg_it.Set(255);
		}
	}

	AssignLabelImage(sub_tmp_img, label, label_img);
	return true;
}


bool radConeSegmentation::ImproveSingleConeSegmentation(BinaryImageType::Pointer seg_img, FloatImageType::Pointer laplacian_img, 
														FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
														FloatImageType::Pointer laplacian_dxy)
{
	BinaryImageType::IndexType lower_index, upper_index;
	for (unsigned int i=0; i<seg_img->GetImageDimension(); i++)
	{
		lower_index[i] = seg_img->GetLargestPossibleRegion().GetSize()[i];
		upper_index[i] = 0;
	}

	BinaryIteratorType seg_it(seg_img, seg_img->GetLargestPossibleRegion());
	for (seg_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it)
	{
		if (seg_it.Get() != 0)
		{
			for (unsigned int i=0; i<seg_img->GetImageDimension(); i++)
			{
				if (lower_index[i] > seg_it.GetIndex()[i])
					lower_index[i] = seg_it.GetIndex()[i];
				if (upper_index[i] < seg_it.GetIndex()[i])
					upper_index[i] = seg_it.GetIndex()[i];
			}
		}
	}

	double ratio = 1;
	//Here, segment each cone
	BinaryImageType::SizeType region_size;
	BinaryImageType::RegionType cone_region;
	for (unsigned int j=0; j<2; j++)
	{
		lower_index[j] = std::max<int>(0, (int)lower_index[j]-(int)(upper_index[j]-lower_index[j])*ratio);
		upper_index[j] = std::min<int>(seg_img->GetLargestPossibleRegion().GetSize()[j], 
				(int)lower_index[j]+(int)(upper_index[j]-lower_index[j])*(1.0+ratio));
		region_size[j] = upper_index[j] - lower_index[j];
	}
	cone_region.SetIndex(lower_index);
	cone_region.SetSize(region_size);

	FloatImageType::Pointer sub_img = ExtractImage<FloatImageType>(InputImages.first, cone_region);
	BinaryImageType::Pointer sub_cone_img = ExtractImage<BinaryImageType>(seg_img, cone_region);

	DoublePointArray cell_contour;
	BinaryImageType::Pointer tmp_res_img; //used for single cone segmentation
	if (ImproveConeSegmentation(sub_img, sub_cone_img, laplacian_img, cell_contour, laplacian_dx, 
		laplacian_dy, laplacian_dxy, NULL, -1, true))
	{
		//assign cone template
		BinaryIteratorType sub_cone_it(sub_cone_img, sub_cone_img->GetLargestPossibleRegion());
		BinaryIteratorType shape_it(ShapePriorImage, sub_cone_img->GetLargestPossibleRegion());
		for (sub_cone_it.GoToBegin(), shape_it.GoToBegin(); !sub_cone_it.IsAtEnd(); ++shape_it, ++sub_cone_it)
		{
			if (sub_cone_it.Get() != 0)
				shape_it.Set(1);
		}

		SegmentationContours.push_back(cell_contour);
		return true;
	}
	else
		return false;
}

bool radConeSegmentation::FindNearestPointPair(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType *labelObject,
											   BinaryImageType::Pointer convex_hull_img, int label, BinaryImageType::IndexType &start_index, 
											   BinaryImageType::IndexType &end_index)
{
	BinaryImageType::IndexType region_start_index;
	BinaryImageType::SizeType region_size;
	BinaryImageType::RegionType label_region;
	if (label == 1) //left side
	{
		region_start_index[0] = labelObject->GetBoundingBox().GetIndex()[0] - labelObject->GetBoundingBox().GetSize()[0]/2 > 0
			? labelObject->GetBoundingBox().GetIndex()[0] - labelObject->GetBoundingBox().GetSize()[0]/2 : 0;
		region_start_index[1] = labelObject->GetBoundingBox().GetIndex()[1];

		region_size[0] = labelObject->GetBoundingBox().GetIndex()[0]-region_start_index[0]
			+ labelObject->GetBoundingBox().GetSize()[0];
		region_size[1] = labelObject->GetBoundingBox().GetSize()[1];
	}
	else if (label == 2)
	{
		region_start_index[0] = labelObject->GetBoundingBox().GetIndex()[0];
		region_start_index[1] = labelObject->GetBoundingBox().GetIndex()[1];

		region_size[0] = region_start_index[0] + 3*labelObject->GetBoundingBox().GetSize()[0]/2 
				< convex_hull_img->GetLargestPossibleRegion().GetSize()[0]
				? 3*labelObject->GetBoundingBox().GetSize()[0]/2 
				: convex_hull_img->GetLargestPossibleRegion().GetSize()[0] - region_start_index[0];
		region_size[1] = labelObject->GetBoundingBox().GetSize()[1];
	}

	double imgW = (double)convex_hull_img->GetLargestPossibleRegion().GetSize()[0];
	double imgH = (double)convex_hull_img->GetLargestPossibleRegion().GetSize()[1];
	if (region_start_index[0] < 0 || region_start_index[0] >= imgW
	|| region_start_index[0] < 1 || region_start_index[1] >= imgH)
		return false;

	start_index.Fill(0);
	end_index.Fill(0);

	label_region.SetIndex(region_start_index);
	label_region.SetSize(region_size);

	unsigned int min_dist = 100000;
	unsigned int cur_dist;
	BinaryIteratorType label_it(convex_hull_img, label_region);
	for (unsigned int i=0; i<labelObject->GetNumberOfPixels(); i++)
	{
		for (label_it.GoToBegin(); !label_it.IsAtEnd(); ++label_it)
		{
			//check the nearest pixels
			if (label_it.Get() == label)
			{
				cur_dist = (labelObject->GetIndex(i)[0] - label_it.GetIndex()[0])
					* (labelObject->GetIndex(i)[0] - label_it.GetIndex()[0])
					+ (labelObject->GetIndex(i)[1] - label_it.GetIndex()[1])
					* (labelObject->GetIndex(i)[1] - label_it.GetIndex()[1]);

				if (cur_dist < min_dist)
				{
					min_dist = cur_dist;
					start_index = labelObject->GetIndex(i);
					end_index = label_it.GetIndex();
				}
			}
		}
	}

	if (start_index[0] == start_index[1] && start_index[1] == end_index[0] 
	&& end_index[0] == end_index[1] && end_index[1] == 0)
		return false;
	else
		return true;
}

void radConeSegmentation::FillImageRegion(BinaryImageType::IndexType & start_index, BinaryImageType::IndexType &end_index, 
										  BinaryImageType::Pointer img)
{
	BinaryImageType::IndexType region_start_index, region_end_index;
	BinaryImageType::SizeType region_size;
	BinaryImageType::RegionType requested_region;

	for (unsigned int i=0; i<2; i++)
	{
		region_start_index[i] = start_index[i] < end_index[i] ? start_index[i] : end_index[i];
		region_end_index[i] = start_index[i] >= end_index[i] ? start_index[i] : end_index[i];

		region_size[i] = region_end_index[i]-region_start_index[i]+1;
	}

	requested_region.SetIndex(region_start_index);
	requested_region.SetSize(region_size);

	itk::ImageRegionIterator<BinaryImageType> it(img, requested_region);
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		it.Set(255);
}

void radConeSegmentation::RefineConePairSegmentation(FloatImageType::Pointer img, FloatImageType::Pointer laplacian_img, 
													 FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
													 FloatImageType::Pointer laplacian_dxy)
{
	BinaryImageType::Pointer tmp_seg_img = AllocateImage<FloatImageType, BinaryImageType>(img);
	tmp_seg_img->FillBuffer(0);
	BinaryIteratorType seg_it(ConvexHullImage, ConvexHullImage->GetLargestPossibleRegion());
	BinaryIteratorType tmp_seg_it(tmp_seg_img, tmp_seg_img->GetLargestPossibleRegion());
	for (seg_it.GoToBegin(), tmp_seg_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it, ++tmp_seg_it)
	{
		if (seg_it.Get() != 0)
			tmp_seg_it.Set(255);
	}

	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer seg_labeling 
		= LabelObjectSegmentation(tmp_seg_img);

	tmp_seg_img->FillBuffer(0);
	bool pair_flag1, pair_flag2, pair_flag3;
	for(unsigned int i = 0; i < seg_labeling->GetNumberOfLabelObjects(); i++)
    {
		pair_flag1 = false;
		pair_flag2 = false;
		pair_flag3 = false;

		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= seg_labeling->GetNthLabelObject(i);
		for (unsigned int j=0; j<labelObject->Size(); j++)
		{
			//make sure the connected region contains left blob, right blob and filled distance gap
			if (ConvexHullImage->GetPixel(labelObject->GetIndex(j)) == 3)
			{
				pair_flag3 = true;
			}

			if (ConvexHullImage->GetPixel(labelObject->GetIndex(j)) == 2)
			{
				pair_flag2 = true;
			}

			if (ConvexHullImage->GetPixel(labelObject->GetIndex(j)) == 1)
			{
				pair_flag1 = true;
			}
		}

		if (!pair_flag3) //no filled distance gap, not region pair connection, skip now
			continue;

		for (unsigned int j=0; j<labelObject->Size(); j++)
		{
			tmp_seg_img->SetPixel(labelObject->GetIndex(j), 255);
		}
	
		BinaryImageType::IndexType start_index, end_index;
		if (!pair_flag1)
		{
			//we are missing the region 1, check the right side to supplement the current tmp_seg_img
			if (!FindNearestPointPair(labelObject, ConvexHullImage, 1, start_index, end_index))
				continue;

			ReassignBinaryImageFromLabelImage(ObjectLabelImage, ObjectLabelImage->GetPixel(end_index), tmp_seg_img);
			FillImageRegion(start_index, end_index, tmp_seg_img);
		}

		if (!pair_flag2)
		{
			//we are missing the region 1, check the right side to supplement the current tmp_seg_img
			if (!FindNearestPointPair(labelObject, ConvexHullImage, 2, start_index, end_index))
				continue;

			ReassignBinaryImageFromLabelImage(ObjectLabelImage, ObjectLabelImage->GetPixel(end_index), tmp_seg_img);
			FillImageRegion(start_index, end_index, tmp_seg_img);
		}

    }

	//fill holes
	typedef itk::BinaryFillholeImageFilter< BinaryImageType > HoleFilterType;
	HoleFilterType::Pointer hole_filter = HoleFilterType::New();
	hole_filter->SetInput( tmp_seg_img );
	hole_filter->SetForegroundValue( itk::NumericTraits< BinaryImageType::PixelType >::max() );
	hole_filter->Update();
	tmp_seg_img = hole_filter->GetOutput();

	ConvexHullPoints.clear();
	seg_labeling = LabelObjectSegmentation(tmp_seg_img);
	vtkSmartPointer<vtkPoints> inPoints = vtkSmartPointer<vtkPoints>::New();
	vector<unsigned int> label_list;
	for(unsigned int i = 0; i < seg_labeling->GetNumberOfLabelObjects(); i++)
    {
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= seg_labeling->GetNthLabelObject(i);

		GetRegionLabelsBelongingToCurrentObject(labelObject, label_list);

		for (unsigned k=0; k<label_list.size(); k++)
		{
			inPoints->Initialize();
			for (unsigned int j=0; j<labelObject->Size(); j++)
			{
				if (ObjectLabelImage->GetPixel(labelObject->GetIndex(j)) == label_list[k])
					inPoints->InsertNextPoint(labelObject->GetIndex(j)[0], labelObject->GetIndex(j)[1], 0.0);
			}

			vtkSmartPointer<vtkPoints> ch_pts = vtkSmartPointer<vtkPoints>::New();
			vtkConvexHull2D::CalculateConvexHull(inPoints, ch_pts, 2.0);
			ConvexHullPoints.push_back(ch_pts);
		}
	}

	BinaryImageType::Pointer convex_hull_img = AllocateImage<FloatImageType, BinaryImageType>(img);
	convex_hull_img->FillBuffer(0);
	ObjectLabelImage->FillBuffer(0);
	FillConvexHullImage(convex_hull_img);
	seg_labeling = LabelObjectSegmentation(convex_hull_img);

	//used for the final segmentation
	FinalSegmentedImage = AllocateImage<BinaryImageType, BinaryImageType>(ConvexHullImage);
	FinalSegmentedImage->FillBuffer(0);
	ShortImageType::Pointer tmp_label_img = AllocateImage<ShortImageType, ShortImageType>(ObjectLabelImage);
	tmp_label_img->FillBuffer(0);

	for (unsigned int i=0; i<seg_labeling->GetNumberOfLabelObjects(); i++)
	{
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= seg_labeling->GetNthLabelObject(i);

		GetRegionLabelsBelongingToCurrentObject(labelObject, label_list);

		for (unsigned int k=0; k<label_list.size(); k++)
		{
			double ratio = 1.0;
			//Here, segment each cone
			BinaryImageType::IndexType low_index, high_index;
			BinaryImageType::SizeType region_size;
			BinaryImageType::RegionType cone_region;

			//determine low index, high index, and region size more precisely
			for (unsigned int j=0; j<2; j++)
			{
				low_index[j] = std::max<int>(0, (int)labelObject->GetBoundingBox().GetIndex()[j]
					-(int)labelObject->GetBoundingBox().GetSize()[j]*ratio);
				high_index[j] = std::min<int>(img->GetLargestPossibleRegion().GetSize()[j], 
					(int)labelObject->GetBoundingBox().GetIndex()[j]+(int)labelObject->GetBoundingBox().GetSize()[j]*(1.0+ratio));
				region_size[j] = high_index[j] - low_index[j];
			}
			cone_region.SetIndex(low_index);
			cone_region.SetSize(region_size);

			FloatImageType::Pointer sub_img = ExtractImage<FloatImageType>(img, cone_region);
			BinaryImageType::Pointer sub_cone_img = ExtractImage<BinaryImageType>(convex_hull_img, cone_region);

			//update sub_cone_img based on the label_list
			BinaryIteratorType sub_cone_img_it(sub_cone_img, sub_cone_img->GetLargestPossibleRegion());
			for (sub_cone_img_it.GoToBegin(); !sub_cone_img_it.IsAtEnd(); ++sub_cone_img_it)
			{
				if (ObjectLabelImage->GetPixel(sub_cone_img_it.GetIndex()) != label_list[k])
					sub_cone_img_it.Set(0);
			}

			DoublePointArray cell_contour;
			if (ImproveConeSegmentation(sub_img, sub_cone_img, laplacian_img, cell_contour, 
				laplacian_dx, laplacian_dy, laplacian_dxy, tmp_label_img, SegmentationContours.size()+1))
			{
				//assign cone template
				BinaryIteratorType sub_cone_it(sub_cone_img, sub_cone_img->GetLargestPossibleRegion());
				BinaryIteratorType shape_it(ShapePriorImage, sub_cone_img->GetLargestPossibleRegion());
				for (sub_cone_it.GoToBegin(), shape_it.GoToBegin(); !sub_cone_it.IsAtEnd(); ++shape_it, ++sub_cone_it)
				{
					if (sub_cone_it.Get() != 0)
						shape_it.Set(1);
				}

				SegmentationContours.push_back(cell_contour);
			}

		}
	}

	DeepCopy<ShortImageType>(tmp_label_img, ObjectLabelImage);
	FinalSegmentedImage = BinaryOpening<BinaryImageType>(FinalSegmentedImage, 2);
}

float radConeSegmentation::SearchLocalDarkPoint(DoublePointType & center_pt, FloatImageType::Pointer split_img, 
												BinaryImageType::Pointer dark_img, bool split_type, 
												ShortPointType & min_pt)
{
	//double degree_interval = 2;
	double angle;
	int start_degree_line, end_degree_line;
	double number_of_points = 0, number_of_total_points = 0;
	double number_of_already_segmented_points = 0;
	float split_val, split_min_val = std::numeric_limits<float>::max();

	if (split_type)
	{
		start_degree_line = -30;
		end_degree_line = 30;
	}
	else
	{
		start_degree_line = 150;
		end_degree_line = 210;
	}

	for (int i=start_degree_line; i<=end_degree_line; i++)
	{
		int search_min = AveragePairDistance*0.5, search_max = AveragePairDistance*1.5;
		ShortPointType pt;
		FloatImageType::IndexType pixelIndex;

		angle = i;
		for (int r=search_min; r<=search_max; r++)
		{
			pt[0] = center_pt[0] + r*cos(angle*PI/180.0);
			pt[1] = center_pt[1] + r*sin(angle*PI/180.0);

			if ((int)(pt[0]+0.5) < 0 || (int)(pt[0]+0.5) >= (int) split_img->GetLargestPossibleRegion().GetSize()[0]
			|| (int)(pt[1]+0.5) < 0 || (int)(pt[1]+0.5) >= (int) split_img->GetLargestPossibleRegion().GetSize()[1])
			{
				number_of_points += 1;
				number_of_total_points += 1;
				continue;
			}

			pixelIndex[0] = (int)(pt[0]+0.5);
			pixelIndex[1] = (int)(pt[1]+0.5);
			
			if (dark_img->GetPixel(pixelIndex) == 0)
			{
				number_of_total_points += 1;
				split_val = split_img->GetPixel(pixelIndex);
				if (split_val < split_min_val)
				{
					split_min_val = split_val;
					min_pt[0] = pixelIndex[0];
					min_pt[1] = pixelIndex[1];
				}
			}
			else
				number_of_already_segmented_points += 1;
		}
	}

	if (number_of_points/number_of_total_points > 0.4 || number_of_already_segmented_points > number_of_total_points)
		return -1;
	else
		return split_min_val;
}

void radConeSegmentation::RegionGrowingSegmentation(ShortPointType & pt, FloatImageType::Pointer split_img, 
													BinaryImageType::Pointer tmp_seg_img)
{
	float timeThreshold = 50;
	typedef  itk::CurvatureAnisotropicDiffusionImageFilter< FloatImageType, FloatImageType > SmoothingFilterType;
	SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
	smoothing->SetTimeStep( 0.125 );
	smoothing->SetNumberOfIterations( 5 );
	smoothing->SetConductanceParameter( 9.0 );
	smoothing->SetInput( split_img );

	typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter< FloatImageType, FloatImageType >  GradientFilterType;
	GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
	gradientMagnitude->SetInput( smoothing->GetOutput() );
	const double sigma = 2;
	gradientMagnitude->SetSigma(  sigma  );
  
	typedef   itk::SigmoidImageFilter< FloatImageType, FloatImageType >  SigmoidFilterType;
	SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	sigmoid->SetOutputMinimum(  0.0  );
	sigmoid->SetOutputMaximum(  1.0  );
	sigmoid->SetInput( gradientMagnitude->GetOutput() );
	sigmoid->SetAlpha( -1 );
	sigmoid->SetBeta(  2  );

	typedef  itk::FastMarchingImageFilter< FloatImageType, FloatImageType >    FastMarchingFilterType;
	FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
	fastMarching->SetInput( sigmoid->GetOutput() );
	typedef FastMarchingFilterType::NodeContainer           NodeContainer;
	typedef FastMarchingFilterType::NodeType                NodeType;
	NodeContainer::Pointer seeds = NodeContainer::New();
	//  Software Guide : EndCodeSnippet
	FloatImageType::IndexType  seedPosition;
	seedPosition[0] = pt[0];
	seedPosition[1] = pt[1];
	NodeType node;
	const double seedValue = 0.0;
	node.SetValue( seedValue );
	node.SetIndex( seedPosition );
	seeds->Initialize();
	seeds->InsertElement( 0, node );
	fastMarching->SetTrialPoints(  seeds  );
	fastMarching->SetOutputSize(tmp_seg_img->GetBufferedRegion().GetSize() );
	const double stoppingTime = std::max(1, SingleGACIterationNumber/2);
	// Software Guide : BeginCodeSnippet
	fastMarching->SetStoppingValue(  stoppingTime  );

	typedef itk::BinaryThresholdImageFilter< FloatImageType, BinaryImageType > ThresholdingFilterType;
	ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
	thresholder->SetLowerThreshold(           0.0  );
	thresholder->SetUpperThreshold( timeThreshold  );
	thresholder->SetOutsideValue(  0  );
	thresholder->SetInsideValue(  255 );
	thresholder->SetInput( fastMarching->GetOutput() );
	thresholder->Update();

	BinaryImageType::Pointer thresholder_img = thresholder->GetOutput();
	BinaryIteratorType thresholder_it(thresholder_img, thresholder_img->GetLargestPossibleRegion());
	BinaryIteratorType tmp_seg_it(tmp_seg_img, tmp_seg_img->GetLargestPossibleRegion());
	for (thresholder_it.GoToBegin(), tmp_seg_it.GoToBegin(); !tmp_seg_it.IsAtEnd(); ++tmp_seg_it, ++thresholder_it)
		tmp_seg_it.Set(thresholder_it.Get());
}

void radConeSegmentation::ComputeExtraConvexHullPoints(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* label_object, 
													   BinaryImageType::Pointer tmp_seg_img, vtkSmartPointer<vtkPoints> convex_hull_points)
{
	vtkSmartPointer<vtkPoints> inPoints = vtkSmartPointer<vtkPoints>::New();
	inPoints->Initialize();
	for(unsigned int i = 0; i < label_object->Size(); i++)
    {
		inPoints->InsertNextPoint(label_object->GetIndex(i)[0], label_object->GetIndex(i)[1], 0);
	}

	BinaryIteratorType tmp_seg_it(tmp_seg_img, tmp_seg_img->GetLargestPossibleRegion());
	for (tmp_seg_it.GoToBegin(); !tmp_seg_it.IsAtEnd(); ++tmp_seg_it)
	{
		if (tmp_seg_it.Get() != 0)
			inPoints->InsertNextPoint(tmp_seg_it.GetIndex()[0], tmp_seg_it.GetIndex()[1], 0);
	}

	vtkConvexHull2D::CalculateConvexHull(inPoints, convex_hull_points, 2.0);
}

void radConeSegmentation::RefineSingleConeSegmentation(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer left_labeling, 
													   itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer right_labeling, 
													   FloatImageType::Pointer laplacian_img, FloatImageType::Pointer laplacian_dx, 
													   FloatImageType::Pointer laplacian_dy, FloatImageType::Pointer laplacian_dxy)
{
	int num_ticks = 5;
	unsigned int blksz = SegmentationPairs.size() / num_ticks;
	unsigned int next_blk = blksz;
	
	ShortPointType min_pt;
	float split_val;
	BinaryImageType::Pointer tmp_seg_img = AllocateImage<FloatImageType, BinaryImageType>(InputImages.first);
	for (unsigned int i=0; i<SegmentationPairs.size(); i++)
	{
		if (i >= next_blk) {
			next_blk += blksz;
			if (num_ticks > 0) {
				emit tick();
				--num_ticks;
			}
		}
		if (SegmentationPairs[i].first == std::numeric_limits<unsigned int>::max())
		{
			split_val = SearchLocalDarkPoint(RightSegmentationCenters[SegmentationPairs[i].second], InputImages.first, 
				FinalSegmentedImage, false, min_pt);
			tmp_seg_img->FillBuffer(0);
			
			if (split_val == -1) //add blob near the boundary
			{
				itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
					= right_labeling->GetNthLabelObject(SegmentationPairs[i].second);
				for (unsigned int j=0; j<labelObject->Size(); j++)
				{
					tmp_seg_img->SetPixel(labelObject->GetIndex(j), 255);
				}

				//if (IsSegmentationValid(tmp_seg_img, BlobSizeThreshold, 0.0))
				if (IsSegmentationValid(tmp_seg_img, ConeRadiusThreshold, ConeRoundnessThreshold))
				//if (IsSegmentationValid(tmp_seg_img, BlobSizeThreshold, ConeRoundnessThreshold))
				{
					for (unsigned int j=0; j<labelObject->Size(); j++)
					{
						FinalSegmentedImage->SetPixel(labelObject->GetIndex(j), 255);
					}

					//get boundary pixels
					DoublePointArray cell_contour;
					ExtractBoundaryPixels(tmp_seg_img, cell_contour, 255);
					//ReportLocations(cell_contour);
					if (!cell_contour.empty())
						SegmentationContours.push_back(cell_contour);
				}
			}
			else if (split_val <= DarkConeThreshold)
			{
				RegionGrowingSegmentation(min_pt, InputImages.first, tmp_seg_img);

				BinaryIteratorType seg_it(ConvexHullImage, ConvexHullImage->GetLargestPossibleRegion());
				BinaryIteratorType tmp_it(tmp_seg_img, tmp_seg_img->GetLargestPossibleRegion());
				for (seg_it.GoToBegin(), tmp_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it, ++tmp_it)
				{
					if (tmp_it.Get() != 0)
						seg_it.Set(1);
				}

				vtkSmartPointer<vtkPoints> convex_hull_points = vtkSmartPointer<vtkPoints>::New();
				itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
					= right_labeling->GetNthLabelObject(SegmentationPairs[i].second);

				if (!IsSegmentationValid(tmp_seg_img, labelObject))
					continue;

				ComputeExtraConvexHullPoints(labelObject, tmp_seg_img, convex_hull_points);
				FillConvexHullImage(tmp_seg_img, convex_hull_points);

				//intentially put here
				if (ImproveSingleConeSegmentation(tmp_seg_img, laplacian_img, laplacian_dx, laplacian_dy, laplacian_dxy))
					ConvexHullPoints.push_back(convex_hull_points);
			}
		}
		if (SegmentationPairs[i].second == std::numeric_limits<unsigned int>::max())
		{
			split_val = SearchLocalDarkPoint(LeftSegmentationCenters[SegmentationPairs[i].first], InputImages.second, 
				FinalSegmentedImage, true, min_pt);
			tmp_seg_img->FillBuffer(0);

			if (split_val == -1) //add blob near the boundary
			{
				itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
					= left_labeling->GetNthLabelObject(SegmentationPairs[i].first);
				for (unsigned int j=0; j<labelObject->Size(); j++)
				{
					tmp_seg_img->SetPixel(labelObject->GetIndex(j), 255);
				}

				if (IsSegmentationValid(tmp_seg_img, ConeRadiusThreshold, ConeRoundnessThreshold))
				{
					for (unsigned int j=0; j<labelObject->Size(); j++)
					{
						FinalSegmentedImage->SetPixel(labelObject->GetIndex(j), 255);
					}

					//get boundary pixels
					DoublePointArray cell_contour;
					ExtractBoundaryPixels(tmp_seg_img, cell_contour, 255);
					if (!cell_contour.empty())
						SegmentationContours.push_back(cell_contour);
				}
			}
			else if (split_val <= DarkConeThreshold)
			{
				RegionGrowingSegmentation(min_pt, InputImages.second, tmp_seg_img);

				BinaryIteratorType seg_it(ConvexHullImage, ConvexHullImage->GetLargestPossibleRegion());
				BinaryIteratorType tmp_it(tmp_seg_img, tmp_seg_img->GetLargestPossibleRegion());
				for (seg_it.GoToBegin(), tmp_it.GoToBegin(); !seg_it.IsAtEnd(); ++seg_it, ++tmp_it)
				{
					if (tmp_it.Get() != 0)
						seg_it.Set(2);
				}

				vtkSmartPointer<vtkPoints> convex_hull_points = vtkSmartPointer<vtkPoints>::New();
				itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
					= left_labeling->GetNthLabelObject(SegmentationPairs[i].first);

				if (!IsSegmentationValid(tmp_seg_img, labelObject))
					continue;

				ComputeExtraConvexHullPoints(labelObject, tmp_seg_img, convex_hull_points);
				FillConvexHullImage(tmp_seg_img, convex_hull_points);

				if (ImproveSingleConeSegmentation(tmp_seg_img, laplacian_img, laplacian_dx, laplacian_dy, laplacian_dxy))
					ConvexHullPoints.push_back(convex_hull_points);
			}
		}
	}

	while (num_ticks > 0) {
		emit tick();
		--num_ticks;
	}

}

FloatImageType::Pointer radConeSegmentation::ComputeLOGImage(FloatImageType::Pointer img, double sigma, bool normalized)
{
	typedef itk::LaplacianRecursiveGaussianImageFilter<FloatImageType, FloatImageType >  filterType;
	filterType::Pointer laplacianFilter = filterType::New();
	laplacianFilter->SetInput( img ); // NOTE: input image type must be double or float
	laplacianFilter->SetSigma(sigma);
	laplacianFilter->SetNormalizeAcrossScale(normalized);
	laplacianFilter->Update();

	return laplacianFilter->GetOutput();
}

void radConeSegmentation::DetectShadingBlobs(FloatImageType::Pointer img, BinaryImageType::Pointer hessian_binary_img)
{
	//step 1: build scale space
	typedef itk::JoinSeriesImageFilter<FloatImageType, FloatImageType3D> JoinSeriesFilterType;
	JoinSeriesFilterType::Pointer joinFilter = JoinSeriesFilterType::New();
	joinFilter->SetOrigin(0.0);
	joinFilter->SetSpacing(1.0);

	JoinSeriesFilterType::Pointer joinFilter1 = JoinSeriesFilterType::New();
	joinFilter1->SetOrigin(0.0);
	joinFilter1->SetSpacing(1.0);

	//build LOG scale space
	for (unsigned int i=0; i< Number_Of_Scale_Levels; i++)
	{
		double sigma = Initial_Scale*pow(Scale_Interval, (double)i);
		FloatImageType::Pointer Ixx = scalespace::DxxImage<FloatImageType, FloatImageType>(img, sigma, true);
		FloatImageType::Pointer Iyy = scalespace::DyyImage<FloatImageType, FloatImageType>(img, sigma, true);
		FloatImageType::Pointer Ixy = scalespace::DxyImage<FloatImageType, FloatImageType>(img, sigma, true);
		FloatImageType::Pointer mult1 = MultipleImages<FloatImageType, FloatImageType, FloatImageType>(Ixx, Iyy);
		FloatImageType::Pointer mult2 = MultipleImages<FloatImageType, FloatImageType, FloatImageType>(Ixy, Ixy);
		FloatImageType::Pointer subtract = SubtractImages<FloatImageType, FloatImageType, FloatImageType>(mult1, mult2);
		joinFilter->SetInput(i, subtract);

		FloatImageType::Pointer tmp_img = ComputeLOGImage(img, sigma);
		joinFilter1->SetInput(i, tmp_img);
	}
	joinFilter->Update();
	FloatImageType3D::Pointer hessian_img = joinFilter->GetOutput();
	joinFilter1->Update();
	FloatImageType3D::Pointer scale_img = joinFilter1->GetOutput();

	//step 2: dark shanding blob selection
	FloatImageType3D::SizeType radius;
	radius.Fill(1);

	itk::NeighborhoodIterator<FloatImageType3D> hessian_it(radius, hessian_img, hessian_img->GetLargestPossibleRegion());
	FloatImageType3D::IndexType maxIndex;
	float max_val, optimized_scale, val;
	bool IsInBounds;
	
	//detect only local maximum
	for (hessian_it.GoToBegin(); !hessian_it.IsAtEnd(); ++hessian_it)
	{
		max_val = hessian_it.GetCenterPixel();
		maxIndex = hessian_it.GetIndex();

		//new added sentence
		for (unsigned int k=0; k<hessian_it.Size(); k++)
		{
			val = hessian_it.GetPixel(k, IsInBounds);

			if (!IsInBounds)
				continue;

			if (val > max_val)
			{
				max_val = val;
				maxIndex = hessian_it.GetIndex(k);
				optimized_scale = Initial_Scale*pow(Scale_Interval, (double)maxIndex[2]);
			}
		}

		if (max_val > HessianThreshold && maxIndex[2] != hessian_img->GetLargestPossibleRegion().GetSize()[2]-1)
		{
			if (scale_img->GetPixel(maxIndex) > 0) 
			{
				BinaryImageType::IndexType tmp_index;
				tmp_index[0] = maxIndex[0];
				tmp_index[1] = maxIndex[1];
				hessian_binary_img->SetPixel(tmp_index, 255);
			}
		}
	}
}

void radConeSegmentation::RefineExtraSegmentations()
{
	ObjectLabelImage->FillBuffer(0);
	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer seg_labeling = LabelObjectSegmentation(FinalSegmentedImage);
	vector<double> radius_arr;
	vector<unsigned int> pixel_num_arr; 

	for(unsigned int i = 0; i < seg_labeling->GetNumberOfLabelObjects(); i++)
    {
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= seg_labeling->GetNthLabelObject(i);

		radius_arr.push_back(labelObject->GetEquivalentSphericalRadius());
		pixel_num_arr.push_back(labelObject->Size());
    }

	std::sort(radius_arr.begin(), radius_arr.end());
	double median_radius = radius_arr.at(radius_arr.size()/2);

	FinalSegmentedImage->FillBuffer(0);

	double radius_selection_thresh1 = 2.0, radius_selection_thresh2 = 2.0;
	unsigned int label_val = 0;
	for(unsigned int i = 0; i < seg_labeling->GetNumberOfLabelObjects(); i++)
    {
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= seg_labeling->GetNthLabelObject(i);

		if (labelObject->GetEquivalentSphericalRadius() > median_radius/radius_selection_thresh1
			&& labelObject->GetEquivalentSphericalRadius() < median_radius*radius_selection_thresh2)
		{
			label_val ++;
			for (unsigned int j=0; j<labelObject->Size(); j++)
			{
				FinalSegmentedImage->SetPixel(labelObject->GetIndex(j), 255);
				ObjectLabelImage->SetPixel(labelObject->GetIndex(j), label_val);
			}
		}
    }

	//check for final segmentation contours
	vector<unsigned int> contour_removal_list;
	DoublePointType center_pt;
	BinaryImageType::IndexType center_index;
	for (unsigned int i=0; i<SegmentationContours.size(); i++)
	{
		center_pt.Fill(0);
		for (unsigned int j=0; j<SegmentationContours[i].size(); j++)
			for (unsigned int k=0; k<2; k++)
				center_pt[k] += SegmentationContours[i][j][k];

		for (unsigned int k=0; k<2; k++)
		{
			center_pt[k] /= SegmentationContours[i].size();
			center_index[k] = (int)(center_pt[k]+0.5);
		}

		if (center_index[0] >= 0 && center_index[0] < (int) FinalSegmentedImage->GetLargestPossibleRegion().GetSize()[0]
			&& center_index[1] >= 0 && center_index[1] < (int) FinalSegmentedImage->GetLargestPossibleRegion().GetSize()[1]
			&& FinalSegmentedImage->GetPixel(center_index) == 0)
			contour_removal_list.push_back(i);
	}

	std::sort(contour_removal_list.begin(), contour_removal_list.end(), std::greater<int>());
	for (vector<unsigned int>::iterator it=contour_removal_list.begin(); it!=contour_removal_list.end(); ++it)
	{
		SegmentationContours.erase(SegmentationContours.begin()+(*it));
	}

	//check connected regions
	std::sort(pixel_num_arr.begin(), pixel_num_arr.end());
	unsigned int median_pixel_num = pixel_num_arr.at(pixel_num_arr.size()/2); //the median number of pixels of a segmented cell region

	//build connection between region and contours
	//update the labeling of the final segmentaiton image since we remove some small regions
	seg_labeling = LabelObjectSegmentation(FinalSegmentedImage);

	vector< vector<unsigned int> > contour_region_indices(seg_labeling->GetNumberOfLabelObjects());
	for (unsigned int i=0; i<SegmentationContours.size(); i++)
	{
		center_pt.Fill(0);
		for (unsigned int j=0; j<SegmentationContours[i].size(); j++)
			for (unsigned int k=0; k<2; k++)
				center_pt[k] += SegmentationContours[i][j][k];

		for (unsigned int k=0; k<2; k++)
		{
			center_pt[k] /= SegmentationContours[i].size();
			center_index[k] = (int)(center_pt[k]+0.5);
		}

		unsigned int val = ObjectLabelImage->GetPixel(center_index);
		if (val != 0)
		{
			contour_region_indices[val-1].push_back(i);
		}
	}

	//search for regions with multiple contours
	//need to improve here..
	vector<DoublePointArray> new_contours;
	vector<unsigned int> deleted_contour_indices;
	//second tmp segmentation image for checking
	BinaryImageType::Pointer tmp_seg_img = AllocateImage<BinaryImageType, BinaryImageType>(FinalSegmentedImage);
	for (unsigned int i=0; i<contour_region_indices.size(); i++)
	{
		if (contour_region_indices[i].size() > 1 && pixel_num_arr[i] < median_pixel_num * 1.2)
		{
			//add addition condition here
			tmp_seg_img->FillBuffer(0);
			ReassignBinaryImageFromLabelImage(ObjectLabelImage, i+1, tmp_seg_img);
			itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer tmp_labeling = LabelObjectSegmentation(tmp_seg_img);

			if (tmp_labeling->GetNumberOfLabelObjects() != 1 || tmp_labeling->GetNthLabelObject(0)->GetRoundness() < 0.8)
				continue;

			//need to merge contours belonging to the current region
			for (unsigned int j=0; j<contour_region_indices[i].size(); j++)
				deleted_contour_indices.push_back(contour_region_indices[i][j]);

			DoublePointArray cell_contour;
			ExtractBoundaryPixels(tmp_seg_img, cell_contour, 255);
			new_contours.push_back(cell_contour);
		}
	}

	std::sort(deleted_contour_indices.begin(), deleted_contour_indices.end(), std::greater<int>());
	for (vector<unsigned int>::iterator it=deleted_contour_indices.begin(); it!=deleted_contour_indices.end(); ++it)
	{
		SegmentationContours.erase(SegmentationContours.begin()+(*it)); //delete old contours
	}
	for (unsigned int i=0; i<new_contours.size(); i++)
		SegmentationContours.push_back(new_contours[i]); //add new contours
}

void radConeSegmentation::ComputeConeThresholds()
{
	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer shape_label_img 
		= LabelObjectSegmentation(FinalSegmentedImage);

	double average_radius = 0, average_roundness = 0;
	unsigned int number_of_objects = 0;
	vector<unsigned int> object_indices;
	for(unsigned int i = 0; i < shape_label_img->GetNumberOfLabelObjects(); i++)
    {
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType* labelObject 
			= shape_label_img->GetNthLabelObject(i);

		if (labelObject->GetBoundingBox().GetIndex()[0] <= AveragePairDistance/2 
			|| labelObject->GetBoundingBox().GetIndex()[1] <= AveragePairDistance/2 
			|| labelObject->GetBoundingBox().GetIndex()[0] + labelObject->GetBoundingBox().GetSize()[0] 
			>= FinalSegmentedImage->GetLargestPossibleRegion().GetSize()[0] - AveragePairDistance/2
			|| labelObject->GetBoundingBox().GetIndex()[1] + labelObject->GetBoundingBox().GetSize()[1] 
			>= FinalSegmentedImage->GetLargestPossibleRegion().GetSize()[1] - AveragePairDistance/2)
			continue;

		//check the number of objects in the current labeling, just in case that objects are touching with each other
		GetRegionLabelsBelongingToCurrentObject(labelObject, object_indices);

		if (object_indices.size() == 1)
		{
			average_radius += labelObject->GetEquivalentSphericalRadius();
			average_roundness += labelObject->GetFlatness();
			number_of_objects ++;
		}
		else
		{
			//it means that the current region potentially contains multiple objects
			BinaryImageType::IndexType start_index = labelObject->GetBoundingBox().GetIndex();
			BinaryImageType::SizeType region_size = labelObject->GetBoundingBox().GetSize();
			BinaryImageType::Pointer tmp_img = AllocateImage<BinaryImageType>(start_index, region_size);
			
			for (unsigned int j=0; j<object_indices.size(); j++)
			{
				tmp_img->FillBuffer(0);

				for (unsigned int k=0; k<labelObject->GetNumberOfPixels(); k++)
				{
					if (ObjectLabelImage->GetPixel(labelObject->GetIndex(k)) == object_indices[j])
					{
						tmp_img->SetPixel(labelObject->GetIndex(k), 255);
					}
				}

				itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer tmp_labels 
					= LabelObjectSegmentation(tmp_img);
				if (tmp_labels->GetNthLabelObject(0)->GetNumberOfPixels() < 8)
					continue;

				average_radius += tmp_labels->GetNthLabelObject(0)->GetEquivalentSphericalRadius();
				average_roundness += tmp_labels->GetNthLabelObject(0)->GetFlatness();
				number_of_objects ++;
			}
		}
    }

	if (number_of_objects == 0)
		return;

	ConeRadiusThreshold = 0.7*average_radius/number_of_objects;
	ConeRoundnessThreshold = 0.7*average_roundness/number_of_objects;

	vector<float>::iterator it;
	BlobSizeThreshold = 0;
	for (it=LeftBlobRadius.begin(); it!=LeftBlobRadius.end(); ++it)
		BlobSizeThreshold += *it;

	for (it=RightBlobRadius.begin(); it!=RightBlobRadius.end(); ++it)
		BlobSizeThreshold += *it;

	BlobSizeThreshold /= (LeftBlobRadius.size()+RightBlobRadius.size());
	BlobSizeThreshold = BlobSizeThreshold * 0.7;
}

void radConeSegmentation::SegmentRgbImage(RGBImageType::Pointer rgbImg)
{
	FloatImageType::IndexType img_start;
	img_start[0] = img_start[1] = 0;
	FloatImageType::SizeType img_size = rgbImg->GetLargestPossibleRegion().GetSize();

	FloatImageType::Pointer first_img = AllocateImage<FloatImageType>(img_start, img_size);
	FloatImageType::Pointer second_img = AllocateImage<FloatImageType>(img_start, img_size);

	RGBConstIteratorType rgb_it(rgbImg, rgbImg->GetLargestPossibleRegion());
	FloatIteratorType first_it(first_img, first_img->GetLargestPossibleRegion());
	FloatIteratorType second_it(second_img, second_img->GetLargestPossibleRegion());

	for (rgb_it.GoToBegin(), first_it.GoToBegin(), second_it.GoToBegin(); !rgb_it.IsAtEnd(); ++rgb_it, ++first_it, ++second_it)
	{
		RGBPixelType pix = rgb_it.Get();
		float gray = floor(pix.GetRed()*0.30 + pix.GetGreen()*0.59 + pix.GetBlue()*0.11 + 0.5);
		if (gray < 0.) gray = 0.;
		else if (gray > 255.) gray = 255.;
		first_it.Set(gray);
		gray = 255. - gray;
		second_it.Set(gray);
	}

	DoublePointArray pts;
	std::pair<FloatImageType::Pointer, FloatImageType::Pointer> imgs;
	imgs.first = first_img;
	imgs.second = second_img;

	SegmentImage(imgs, pts);
}

void radConeSegmentation::SegmentImage(std::pair<FloatImageType::Pointer, FloatImageType::Pointer> & imgs, 
									   DoublePointArray & pts)
{
	InputImages = imgs;
	SeedPoints = pts;
	SegmentationContours.clear();

	HessianImages.first = AllocateImage<FloatImageType, BinaryImageType>(imgs.first);
	HessianImages.second = AllocateImage<FloatImageType, BinaryImageType>(imgs.second);
	HessianImages.first->FillBuffer(0);
	HessianImages.second->FillBuffer(0);

	DetectShadingBlobs(imgs.first, HessianImages.first);
	DetectShadingBlobs(imgs.second, HessianImages.second);

	emit tick();

	ConvexHullImage = AllocateImage<FloatImageType, BinaryImageType>(InputImages.first);
	ConvexHullImage->FillBuffer(0);
	ObjectLabelImage = AllocateImage<FloatImageType, ShortImageType>(InputImages.first);
	ObjectLabelImage->FillBuffer(0);
	ShapePriorImage = AllocateImage<FloatImageType, BinaryImageType>(InputImages.first);
	ShapePriorImage->FillBuffer(0);

	DarkSegmentationBlobs.first = AllocateImage<FloatImageType, BinaryImageType>(InputImages.first);
	DarkSegmentationBlobs.first->FillBuffer(0);
	DarkSegmentationBlobs.second = AllocateImage<FloatImageType, BinaryImageType>(InputImages.second);
	DarkSegmentationBlobs.second->FillBuffer(0);
	LevelsetSegmentation(InputImages.first, HessianImages.first, 1, DarkSegmentationBlobs.first);
	LevelsetSegmentation(InputImages.second, HessianImages.second, 2, DarkSegmentationBlobs.second);

	emit tick();

	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer left_labeling = LabelObjectSegmentation(DarkSegmentationBlobs.first);
	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer right_labeling = LabelObjectSegmentation(DarkSegmentationBlobs.second);
	
	FloatImageType::Pointer laplacian_img = ComputeLOGImage(InputImages.first, 3.0, false);
	FloatImageType::Pointer laplacian_img_dx = scalespace::DxImage<FloatImageType, FloatImageType>(laplacian_img, 1, false);
	FloatImageType::Pointer laplacian_img_dy = scalespace::DyImage<FloatImageType, FloatImageType>(laplacian_img, 1, false);
	FloatImageType::Pointer laplacian_img_dxy = scalespace::DxyImage<FloatImageType, FloatImageType>(laplacian_img, 1, false);

	ExpandSegmentationPairs(left_labeling, right_labeling);
	RefineConePairSegmentation(InputImages.first, laplacian_img, laplacian_img_dx, laplacian_img_dy, laplacian_img_dxy);
	ComputeConeThresholds();
	
	emit tick();

	RefineSingleConeSegmentation(left_labeling, right_labeling, laplacian_img, laplacian_img_dx, laplacian_img_dy, laplacian_img_dxy);
	RefineExtraSegmentations();

}
