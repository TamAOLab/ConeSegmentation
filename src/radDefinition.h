/*
 *  radDefinition.h
 *  
 *
 *  Created by jianfei liu on 8/24/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef RADDEFINITION_H
#define RADDEFINITION_H

#include <QColor>

#include <itkRGBPixel.h>
#include "itkImage.h"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkNeighborhoodIterator.h>
#include <itkPoint.h>
#include <itkPointSet.h>
#include <vtkSmartPointer.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <vector>
#include <string>
#include <stdio.h>
#include <stack>

using namespace std;

typedef float                               FloatPixelType;
typedef unsigned char                       BinaryPixelType;
typedef short                               ShortPixelType;
typedef itk::RGBPixel<unsigned char>        RGBPixelType;

typedef itk::CovariantVector<float, 2>		VectorPixelType;

typedef itk::Image<FloatPixelType, 2>		FloatImageType;
typedef itk::Image<BinaryPixelType, 2>		BinaryImageType;
typedef itk::Image<ShortPixelType, 2>		ShortImageType;
typedef itk::Image<RGBPixelType, 2>         RGBImageType;

typedef itk::Image<FloatPixelType, 3>		FloatImageType3D;

typedef itk::Image< VectorPixelType, 2>	VectorImageType;

typedef itk::ImageRegionIteratorWithIndex<FloatImageType>       FloatIteratorType;
typedef itk::ImageRegionConstIteratorWithIndex<FloatImageType>  FloatConstIteratorType;
typedef itk::ImageRegionIteratorWithIndex<BinaryImageType>      BinaryIteratorType;
typedef itk::ImageRegionIteratorWithIndex<ShortImageType>       ShortIteratorType;
typedef itk::ImageRegionConstIteratorWithIndex<RGBImageType>  RGBConstIteratorType;

typedef itk::Point<short, 2>                ShortPointType;
typedef std::vector< ShortPointType >		ShortPointArray;

typedef itk::Point<double, 2>               DoublePointType;
typedef std::vector< DoublePointType >		DoublePointArray;

struct ConeSegmentationParameters
{
	double HessianThreshold;				// 100 (10..500, 10) Hessian Response
	int GACIterationNumber;					// 10 (1..150, 1) Pair/Single GAC Iteration #
	inline void loadDefaults()
	{
		HessianThreshold = 100.;
		GACIterationNumber = 10;
	}
};

namespace ConeSegmentation {
	const double diffusionTimeStep = 1.0;
	const unsigned int diffusionNumberOfSteps = 10;
	const double diffusionWeight = 1.0;
	const float scale = 2.0;
	const double DarkConeThresholdRatio = 1.5;
	const int PairDistance = 30;
}

enum MouseOperation { Mouse_Normal, Mouse_Add_Point, Mouse_Delete_Point, Mouse_Add_Contour, Mouse_Delete_Contour, 
	Mouse_Delete_Single_Contour, Mouse_Edit_Contour };
const double PI = 3.141592653589793238463;
enum ShapeOperator { Shape_Delete, Shape_Add };

class MarkerBase
{
public:
	MarkerBase()
	{
		Initialize();
	}

	void Initialize()
	{
		is_visible = true;
		marker_center[0] = marker_center[1] = 0;
	}

	bool IsPointInsideBox(DoublePointType & pt, DoublePointType & low_box_pt, DoublePointType & high_box_pt)
	{
		if (pt[0] >= low_box_pt[0] && pt[0] <= high_box_pt[0]
			&& pt[1] >= low_box_pt[1] && pt[1] <= high_box_pt[1])
			return true;
		else
			return false;
	}

	bool is_visible;
	DoublePointType marker_center;
};

class ContourMarker : public MarkerBase
{
public:
	ContourMarker()
	{
		Initialize();
	}

	bool is_edited;
	DoublePointType marker_bounding_pts[2]; //low left, and top right points
	DoublePointArray marker_contours;
	vector<BinaryImageType::IndexType> marker_regions;

	void Initialize()
	{
		is_edited = false;
		marker_bounding_pts[0][0] = marker_bounding_pts[0][1]
			= marker_bounding_pts[1][0] = marker_bounding_pts[1][1] = 0;
		marker_contours.clear();
		marker_regions.clear();
	}

	void ComputeShapeInfor() //include center, bounding box, etc
	{
		if (marker_contours.empty())
			return;

		marker_center[0] = marker_center[1] = 0;
		marker_bounding_pts[0][0] = marker_bounding_pts[0][1] = std::numeric_limits<double>::max();
		marker_bounding_pts[1][0] = marker_bounding_pts[1][1] = std::numeric_limits<double>::min();

		DoublePointArray::iterator it;
		for (it = marker_contours.begin(); it != marker_contours.end(); ++it)
		{
			for (unsigned int i = 0; i<2; i++)
			{
				marker_center[i] += (*it)[i];

				if (marker_bounding_pts[0][i] >(*it)[i])
					marker_bounding_pts[0][i] = (*it)[i];
				if (marker_bounding_pts[1][i] < (*it)[i])
					marker_bounding_pts[1][i] = (*it)[i];
			}
		}

		for (unsigned int i = 0; i<2; i++)
			marker_center[i] /= marker_contours.size();
	}

	template <class TImageType>
	void ComputeShapeRegion(typename TImageType::Pointer img)
	{
		if (marker_contours.empty() || !img)
			return;

		//step 1: add new contour to cone segmentation
		BinaryImageType::IndexType pixelIndex;
		vector<float> tmp_contour_pts;

		for (DoublePointArray::iterator it = marker_contours.begin(); it != marker_contours.end(); ++it)
		{
			img->TransformPhysicalPointToIndex(*it, pixelIndex);
			for (unsigned int k = 0; k<2; k++)
				tmp_contour_pts.push_back(pixelIndex[k]);
		}

		int  nodes, i1, j1, swap;
		int polyCorners = tmp_contour_pts.size() / 2;
		int nodeX[1500];
		int IMAGE_RIGHT = img->GetLargestPossibleRegion().GetSize()[0], IMAGE_LEFT = 0;
		//  Loop through the rows of the image.
		marker_regions.clear();
		for (unsigned int y = 0; y<img->GetLargestPossibleRegion().GetSize()[1]; y++)
		{
			//  Build a list of nodes.
			nodes = 0; j1 = polyCorners - 1;
			for (i1 = 0; i1<polyCorners; i1++)
			{
				if ((tmp_contour_pts[2 * i1 + 1] < (double)y && tmp_contour_pts[2 * j1 + 1] >= (double)y)
					|| (tmp_contour_pts[2 * j1 + 1] < (double)y && tmp_contour_pts[2 * i1 + 1] >= (double)y))
				{
					nodeX[nodes++] = (int)(tmp_contour_pts[2 * i1] + (y - tmp_contour_pts[2 * i1 + 1])
						/ (tmp_contour_pts[2 * j1 + 1] - tmp_contour_pts[2 * i1 + 1])
						*(tmp_contour_pts[2 * j1] - tmp_contour_pts[2 * i1]));
				}

				j1 = i1;
			}

			//  Sort the nodes, via a simple ?Bubble? sort.
			i1 = 0;
			while (i1<nodes - 1)
			{
				if (nodeX[i1]>nodeX[i1 + 1])
				{
					swap = nodeX[i1];
					nodeX[i1] = nodeX[i1 + 1];
					nodeX[i1 + 1] = swap;
					if (i1) i1--;
				}
				else
				{
					i1++;
				}
			}

			//  Fill the pixels between node pairs.
			for (i1 = 0; i1<nodes; i1 += 2)
			{
				if (nodeX[i1] >= IMAGE_RIGHT)
					break;

				if (nodeX[i1 + 1]> IMAGE_LEFT)
				{
					if (nodeX[i1] < IMAGE_LEFT)
						nodeX[i1] = IMAGE_LEFT;

					if (nodeX[i1 + 1]> IMAGE_RIGHT)
						nodeX[i1 + 1] = IMAGE_RIGHT;

					for (j1 = nodeX[i1]; j1<nodeX[i1 + 1]; j1++)
					{
						pixelIndex[0] = j1;
						pixelIndex[1] = y;

						//img->SetPixel(pixelIndex, 255);	
						marker_regions.push_back(pixelIndex);
					}
				}
			}
		}
	}

	void UpdateShapeRegion(BinaryImageType::Pointer img, int val)
	{
		BinaryImageType::IndexType pixelIndex;
		for (unsigned int i = 0; i<marker_regions.size(); i++)
		{
			pixelIndex = marker_regions[i];
			img->SetPixel(pixelIndex, val);
		}
	}

	bool withinImageRegion(itk::ImageRegion<2U> & rgn) {
		double right = (double)rgn.GetSize()[0] - 1.;
		double bottom = (double)rgn.GetSize()[1] - 1.;
		if (marker_bounding_pts[0][0] <= 0 || marker_bounding_pts[0][1] <= 0 ||
			marker_bounding_pts[1][0] > right || marker_bounding_pts[1][1] > bottom)
			return false;
		return true;
	}
};

class StackInformation
{
public:
	ShapeOperator shape_operator;
	vector<int> shape_list;

	StackInformation(ShapeOperator op, int id)
	{
		shape_operator = op;
		shape_list.push_back(id);
	}

	StackInformation(ShapeOperator op, vector<int> ids)
	{
		shape_operator = op;
		shape_list.clear();
		shape_list.assign(ids.begin(), ids.end());
	}
};


class MarkerInformation 
{
public:

	MarkerInformation()
	{
		Initialize();
		ClearStack(contour_operator_stack);
	}
	~MarkerInformation()
	{
	}

	RGBImageType::Pointer split_image;
	pair<string, string> split_file_names;
	vector< ContourMarker > cone_contour_markers;

	stack<StackInformation> contour_operator_stack;

	ConeSegmentationParameters segment_params;
	
	void ClearStack(stack<StackInformation> &marker_stack)
	{
		while ( !marker_stack.empty() )
		{
			marker_stack.pop();
		}
	}

	void Initialize() 
	{
		split_image = NULL;
		split_file_names.first.clear();
		split_file_names.second.clear();
		cone_contour_markers.clear();
		segment_params.loadDefaults();
	}
};

class radColor : public QColor
{
private:
	double arr[3];
public:
	double *toDoubleArr()
	{
		arr[0] = red() / 255.;
		arr[1] = green() / 255.;
		arr[2] = blue() / 255.;
		return arr;
	}
	std::string toStdString()
	{
		return name().toStdString();
	}
	radColor & operator = (const QColor &color)
	{
		if (color.isValid())
			setRgb(color.red(), color.green(), color.blue());
		return *this;
	}
	radColor & operator = (const QString & name)
	{
		QColor color(name);
		if (color.isValid())
			setRgb(color.red(), color.green(), color.blue());
		return *this;
	}
};

struct MarkerSystemSettings
{
	int visibility_contour, visibility_region, visibility_center;
	int size_contour;
	double size_center;
	radColor contour_color;
	radColor region_color;
	radColor center_color;
	double region_opacity;

	MarkerSystemSettings()
	{
		visibility_contour = true;
		visibility_center = true;
		visibility_region = false;

		size_contour = 2;
		size_center = 4.;
		contour_color.setNamedColor(QString("#00ff00"));
		center_color.setNamedColor(QString("#00ff00"));
		region_color.setNamedColor(QString("#756bb1"));
		region_opacity = 0.33;
	}
};

#endif // RADDEFINITION_H
