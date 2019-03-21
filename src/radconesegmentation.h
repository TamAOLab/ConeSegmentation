/*
 *  radConeSegmentation.h
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef radConeSegmentation_H
#define radConeSegmentation_H

#include <QObject>

#include "radImgFunc.h"
#include <itkVoronoiDiagram2DGenerator.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <itkFastMarchingImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include "itkConfidenceConnectedImageFilter.h"
#include <itkJoinSeriesImageFilter.h>
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkListSample.h"
#include "itkHistogram.h"
#include "itkSampleToHistogramFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryContourImageFilter.h"

#include <itkLabelMap.h>
#include <itkShapeLabelObject.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkBinaryFillholeImageFilter.h>

// Math
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPolygon.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkConvexHull2D.h>

#include "radScaleSpace.h"

// using namespace cv;

typedef float HistogramMeasurementType;
const unsigned int NumberOfHistogramComponents = 1;
typedef itk::Statistics::Histogram< HistogramMeasurementType > HistogramType;

class radConeSegmentation : public QObject
{
	Q_OBJECT
private:

	std::pair<FloatImageType::Pointer, FloatImageType::Pointer> InputImages;

	BinaryImageType::Pointer ConvexHullImage, FinalSegmentedImage;
	BinaryImageType::Pointer ShapePriorImage;
	ShortImageType::Pointer ObjectLabelImage; //used for the detecting the issue of segmentation touches
	std::pair<BinaryImageType::Pointer, BinaryImageType::Pointer> HessianImages;
	std::pair<BinaryImageType::Pointer, BinaryImageType::Pointer> DarkSegmentationBlobs;
	DoublePointArray SeedPoints;
	vector<float> SeedScales;
	vector< DoublePointArray > SegmentationContours;

	double DiffusionTimeStep;
	unsigned int NumberofSteps;
	double DiffusionWeight;
	float m_Scale; // scale for detecting gradient information via DoG
	double HessianThreshold;
	double DarkConeThresholdRatio;
	int PairGACIterationNumber, SingleGACIterationNumber;

	void LevelsetSegmentation(FloatImageType::Pointer, BinaryImageType::Pointer, int, BinaryImageType::Pointer);

	FloatImageType::Pointer ComputeImageMagn(FloatImageType::Pointer, float);
	FloatImageType::Pointer ComputeLOGImage(FloatImageType::Pointer, double, bool normalized = true);

	//determine paired blobs
	int PairDistance;
	DoublePointArray LeftSegmentationCenters, RightSegmentationCenters;
	vector<float> LeftBlobRadius, RightBlobRadius;
	vector< pair<unsigned int, unsigned int> > SegmentationPairs;
	double AveragePairDistance;
	double DarkConeThreshold, DarkConeThreshold1;
	
	void ExpandSegmentationPairs(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer, 
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer);
	void ComputeSegmentationCenters(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer, 
		DoublePointArray &, vector<float> &);
	void DecidePairsAccordingToCenters();
	int SelectLeftBlob(int, int, vector<bool> &);
	itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer LabelObjectSegmentation(BinaryImageType::Pointer);
	void FillBlobsBasedOnDistance(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer, 
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer);

	float ComputeDistanceBetweenTwoBlobs(BinaryImageType::Pointer, FloatImageType::Pointer); 
	void FillGapBetweenBlobs(FloatImageType::Pointer, FloatImageType::Pointer, BinaryImageType::Pointer, float, 
		DoublePointType &, DoublePointType &, float, float);
	//connect paired blobs

	//Refine cone segmentation, only work on cones with dark and bright regions in this step
	double ConeRadiusThreshold, ConeRoundnessThreshold, BlobSizeThreshold;
	vector< vtkSmartPointer<vtkPoints> > ConvexHullPoints;
	void RefineConePairSegmentation(FloatImageType::Pointer, FloatImageType::Pointer, 
		FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
		FloatImageType::Pointer laplacian_dxy);
	bool FindNearestPointPair(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType*, 
		BinaryImageType::Pointer, int, BinaryImageType::IndexType &, BinaryImageType::IndexType &);
	void ComputeConeThresholds();
	void FillConvexHullImage(BinaryImageType::Pointer);
	void FillConvexHullImage(BinaryImageType::Pointer, vtkSmartPointer<vtkPoints>);
	void FillContourImage(DoublePointArray &, BinaryImageType::Pointer);
	void FillImageRegion(BinaryImageType::IndexType &, BinaryImageType::IndexType &, BinaryImageType::Pointer);
	bool ImproveConeSegmentation(FloatImageType::Pointer, BinaryImageType::Pointer, FloatImageType::Pointer, 
		DoublePointArray &, FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
		FloatImageType::Pointer laplacian_dxy, ShortImageType::Pointer label_img = NULL, 
		int label = -1, bool determination_flag = false);
	bool IsSegmentationValid(BinaryImageType::Pointer, BinaryImageType::Pointer); //check segmentation based on the overlap with the existing segmentation
	bool IsSegmentationValid(BinaryImageType::Pointer, double, double); //check segmentation based on the region size and roundness
	bool IsSegmentationValid(BinaryImageType::Pointer, 
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType*); //check segmentation based on the blob sizes that are used for convex hull
	void LabelBinaryImage(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType*, short,
		ShortImageType::Pointer);
	void GetRegionLabelsBelongingToCurrentObject(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType*, 
		vector<unsigned int> &);
	void AssignLabelImage(BinaryImageType::Pointer, int, ShortImageType::Pointer); //used for assign images
	void ReassignBinaryImageFromLabelImage(ShortImageType::Pointer, int, BinaryImageType::Pointer);

	float BilinearInterpolate(FloatImageType::Pointer img, double x, double y);
	float BicubicInterpolate(FloatImageType::Pointer img, FloatImageType::Pointer img_dx, 
		FloatImageType::Pointer img_dy, FloatImageType::Pointer img_dxy, double x, double y);
	bool SearchForSubPixelLocation(BinaryImageType::IndexType &, BinaryImageType::IndexType &, 
		FloatImageType::Pointer, FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
		FloatImageType::Pointer laplacian_dxy, int, int, double, DoublePointType &);
	void SubPixelImprovement(FloatImageType::Pointer, FloatImageType::Pointer, BinaryImageType::Pointer, DoublePointArray &, 
		FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
		FloatImageType::Pointer laplacian_dxy);
	bool IsCurrentIndexInRightPosition(BinaryImageType::IndexType, vector< BinaryImageType::IndexType > &);
	double ComputeDarkBlobHistogram();
	void RefineSingleConeSegmentation(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer, 
		itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer, FloatImageType::Pointer, 
		FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
		FloatImageType::Pointer laplacian_dxy);
	float SearchLocalDarkPoint(DoublePointType &, FloatImageType::Pointer, BinaryImageType::Pointer, 
		bool, ShortPointType &);
	void RegionGrowingSegmentation(ShortPointType &, FloatImageType::Pointer, BinaryImageType::Pointer);
	void ComputeExtraConvexHullPoints(itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::LabelObjectType*, 
		BinaryImageType::Pointer, vtkSmartPointer<vtkPoints>);
	bool ImproveSingleConeSegmentation(BinaryImageType::Pointer, FloatImageType::Pointer, 
		FloatImageType::Pointer laplacian_dx, FloatImageType::Pointer laplacian_dy,
		FloatImageType::Pointer laplacian_dxy);
	void ExtractBoundaryPixels(BinaryImageType::Pointer, DoublePointArray &, int);

	//snake improvement for subpixel refine
	vnl_matrix<double> ComputeP(double, double, double, double) throw (int);
	vnl_vector<double> SampleImage(vnl_vector<double>, vnl_vector<double>, VectorImageType::Pointer, int);
	bool SnakeRefinement(FloatImageType::Pointer, DoublePointArray &);

	void DetectShadingBlobs(FloatImageType::Pointer, BinaryImageType::Pointer);

	void ExpandConvexHull(BinaryImageType::Pointer, itk::LabelMap< itk::ShapeLabelObject<unsigned long, 2> >::Pointer);

	void RefineExtraSegmentations();

signals:
	void tick();
	
public:
	
	radConeSegmentation();
	~radConeSegmentation();

	void SetScale(float scale) {m_Scale = scale;}
	void SegmentRgbImage(RGBImageType::Pointer rgbImg);
	void SegmentImage(std::pair<FloatImageType::Pointer, FloatImageType::Pointer> &, DoublePointArray &);
	void SetParameters(double time_step = 1.0, unsigned int iteration_step = 10, double weight = 1.0);
	void SetParameters(ConeSegmentationParameters & params);
	inline BinaryImageType::Pointer GetConvexHullImage() {return ConvexHullImage;}
	inline BinaryImageType::Pointer GetShapePriorImage() {return ShapePriorImage;}
	inline BinaryImageType::Pointer GetFinalSegmentedImage() {return FinalSegmentedImage;}
	
	inline void SetSeedScales(vector<float> & s) {SeedScales = s;}
	inline void SetHessianThreshold(double val) {HessianThreshold = val;}
	inline void SetDarkConeThresholdRatio(double val) {DarkConeThresholdRatio = val;}
	inline void SetPairGACIterationNumber(int val) {PairGACIterationNumber = val;}
	inline void SetSingleGACIterationNumber(int val) {SingleGACIterationNumber = val;}
	inline void SetPairDistance(int val) {PairDistance = val;}

	DoublePointArray & GetLeftSegmentationCenters() {return LeftSegmentationCenters;} 
	DoublePointArray & GetRightSegmentationCenters() {return RightSegmentationCenters;}
	vector< pair<unsigned int, unsigned int> > & GetSegmentationPairs() {return SegmentationPairs;}
	vector< vtkSmartPointer<vtkPoints> > & GetConvexHullPoints() {return ConvexHullPoints;}
	vector< DoublePointArray > & GetSegmentationContours() {return SegmentationContours;}

};

#endif // segFileIO_H

