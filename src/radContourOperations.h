/*
 *  radContourOperations.h
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef radContourOperations_H
#define radContourOperations_H

#include "radImgFunc.h"
#include <itkSpatialObjectToImageFilter.h>

// Math
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

struct ContourSegmentInformation 
{
	int contour_id;
	vector< vector<unsigned int> > connected_contour_segments;
	vector<int> contour_segment_labels;

	void Initialize()
	{
		contour_id = -1;
		connected_contour_segments.clear();
		contour_segment_labels.clear();
	}
};

class radContourOperations
{
private:
	
	static void FindIntersectedContours(DoublePointArray &, MarkerInformation &, 
		vector<int> &);
    static bool IsPartlyIntersected(DoublePointArray &, DoublePointArray &);

public:
	
	radContourOperations()
	{
	}

	~radContourOperations()
	{
	}
   
//	static int RemoveConeMarker(DoublePointType &, double, MarkerInformation &);
//	static void AddConeMarker(DoublePointType &, int, MarkerInformation &);
	static void AddConeContour(DoublePointArray &, MarkerInformation &);
	static void RemoveConeContour(DoublePointArray &, MarkerInformation &, vector<int> &);
	static int GetEditedConeContour(DoublePointType &, MarkerInformation &);
	static void UpdateEditedConeContour(DoublePointArray &);
};

#endif // segFileIO_H

