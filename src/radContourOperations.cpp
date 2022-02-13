/*
 *  radContourOperations.cpp
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "itkLabelContourImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

#include "radContourOperations.h"
#include "QuickView.h"

 //-- "Winding Number" algorithm for the inclusion of a point in polygon, adapted from:
 // http://geomalgorithms.com/a03-inclusion.html

 // tests if a point is Left|On|Right of an infinite line.
 //    Input:  three points P0, P1, and P2
 //    Return: >0 for P2 left of the line through P0 and P1
 //            =0 for P2  on the line
 //            <0 for P2  right of the line
static inline int isLeft(DoublePointType P0, DoublePointType P1, DoublePointType P2)
{
	return ((P1[0] - P0[0]) * (P2[1] - P0[1])
		- (P2[0] - P0[0]) * (P1[1] - P0[1]));
}

static bool IsPointInsidePolygon(DoublePointType P, DoublePointArray& vertices)
{
	int n = (int)vertices.size();
	if (n == 0) return false;

	// V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
	std::vector<DoublePointType> V(vertices);
	V.push_back(vertices[0]);

	int wn = 0;    // the  winding number counter (=0 only when P is outside)

	// loop through all edges of the polygon
	for (int i = 0; i < n; i++) {		// edge from V[i] to  V[i+1]
		if (V[i][1] <= P[1]) {			// start y <= P.y
			if (V[i + 1][1] > P[1])		// an upward crossing
				if (isLeft(V[i], V[i + 1], P) > 0)	// P left of  edge
					++wn;				// have  a valid up intersect
		}
		else {							// start y > P.y (no test needed)
			if (V[i + 1][1] <= P[1])		// a downward crossing
				if (isLeft(V[i], V[i + 1], P) < 0)	// P right of  edge
					--wn;				// have  a valid down intersect
		}
	}
	return wn != 0;
}

//-- End of "Winding Number" algorithm

void radContourOperations::AddConeContour(DoublePointArray &contour_pts, MarkerInformation &marker_infor)
{
	ContourMarker marker_contour;
	marker_contour.Initialize();
	
	//assign cone contours
	marker_contour.marker_contours = contour_pts;
	marker_contour.ComputeShapeInfor();
	marker_contour.ComputeShapeRegion<RGBImageType>(marker_infor.split_image);
	marker_infor.cone_contour_markers.push_back(marker_contour);
}

int radContourOperations::GetEditedConeContour(DoublePointType &pt, MarkerInformation &marker_infor)
{
	vector<unsigned int> edited_contour_ids;
	ContourMarker *marker_contour;
	for (unsigned int i=0; i<marker_infor.cone_contour_markers.size(); i++)
	{
		marker_contour = &(marker_infor.cone_contour_markers[i]);
		if (!marker_contour->is_visible)
			continue;

		if (IsPointInsidePolygon(pt, marker_contour->marker_contours))
		{
			edited_contour_ids.push_back(i);
		}	
	}

	if (edited_contour_ids.empty())
		return -1;

	//search for the contour with the shortest distance from the center to the pt
	double min_dist = std::numeric_limits<double>::max();
	double cur_dist;
	unsigned int res_id;
	for (unsigned int i=0; i<edited_contour_ids.size(); i++)
	{
		marker_contour = &(marker_infor.cone_contour_markers[edited_contour_ids[i]]);
		marker_contour->ComputeShapeInfor();

		cur_dist = pt.EuclideanDistanceTo(marker_contour->marker_center);
		if (min_dist > cur_dist)
		{
			min_dist = cur_dist;
			res_id = edited_contour_ids[i];
		}
	}

	return res_id;
}

void radContourOperations::RemoveConeContour(DoublePointArray &contour_pts, MarkerInformation &marker_infor, 
											 vector<int> &deleted_contour_ids)
{
	deleted_contour_ids.clear();
	ContourMarker *marker_contour;
	for (unsigned int i=0; i<marker_infor.cone_contour_markers.size(); i++)
	{
		marker_contour = &(marker_infor.cone_contour_markers[i]);
		if (!marker_contour->is_visible)
			continue;

		/*for (unsigned int j=0; j<marker_contour->marker_contours.size(); j++)
		{
			if (IsPointInsidePolygon(marker_contour->marker_contours[j], contour_pts))
			{
				deleted_contour_ids.push_back(i);
				break;
			}
		}*/

		if (IsPointInsidePolygon(marker_contour->marker_center, contour_pts))
			deleted_contour_ids.push_back(i);
	}

	if (deleted_contour_ids.empty())
		return;

	//first update all deleted contours;
	for (unsigned int i=0; i<deleted_contour_ids.size(); i++)
	{
		marker_contour = &(marker_infor.cone_contour_markers[deleted_contour_ids[i]]);
		marker_infor.cone_contour_markers[deleted_contour_ids[i]].is_visible = false;
	}
}

bool radContourOperations::IsPartlyIntersected(DoublePointArray & contour1, DoublePointArray & contour2)
{
	bool intersection_flag1 = false, intersection_flag2 = false;
	for (unsigned int i=0; i<contour1.size(); i++)
	{
		if (IsPointInsidePolygon(contour1[i], contour2))
			intersection_flag1 = true;
	}

	for (unsigned int i=0; i<contour2.size(); i++)
	{
		if (IsPointInsidePolygon(contour2[i], contour1))
			intersection_flag2 = true;
	}

	if (intersection_flag1 && intersection_flag2)
		return true;
	else 
		return false;
}

void radContourOperations::FindIntersectedContours(DoublePointArray & contour_pts,
												   MarkerInformation &marker_infor,
												   vector<int> & res_contours)
{
	res_contours.clear();
	
	//if the current contour totally inside or outside the existing contours, remaining...
	ContourMarker *marker_contour;

	for (unsigned int i=0; i<marker_infor.cone_contour_markers.size(); i++)
	{
		marker_contour = &(marker_infor.cone_contour_markers[i]);
		if (!marker_contour->is_visible)
			continue;

		//check the intersection
		if (!IsPartlyIntersected(contour_pts, marker_contour->marker_contours))
			continue;

		//if yes, we will determine the intesected contours & regions
		res_contours.push_back(i);
	}
}

void radContourOperations::UpdateEditedConeContour(DoublePointArray & contour_pts)
{
	if (contour_pts.size() < 5)
		return;

	const int num_of_control_pts = 20;
	double perimeter = 0;
	for (unsigned int i=0; i<contour_pts.size(); i++)
	{
		if (i == contour_pts.size()-1)
			perimeter += contour_pts[i].EuclideanDistanceTo(contour_pts[0]);
		else
			perimeter += contour_pts[i].EuclideanDistanceTo(contour_pts[i+1]);
	}
	double average_spacing = perimeter/num_of_control_pts;

	DoublePointArray updated_contour_pts;
	updated_contour_pts.push_back(contour_pts[0]);
	double cur_dist = 0;
	DoublePointType cur_pt;
	cur_pt = contour_pts[0];

	for (unsigned int i=1; i<contour_pts.size(); i++)
	{
		cur_dist += contour_pts[i-1].EuclideanDistanceTo(contour_pts[i]);
		if (cur_dist > average_spacing)
		{
			updated_contour_pts.push_back(contour_pts[i]);
			cur_dist = 0;
		}
	}

	contour_pts.clear();
	contour_pts.assign(updated_contour_pts.begin(), updated_contour_pts.end());
}