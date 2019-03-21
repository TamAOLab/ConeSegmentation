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

double Angle2D(DoublePointType p1, DoublePointType p2)
{
	double dtheta, theta1, theta2;

	theta1 = atan2(p1[1], p1[0]);
	theta2 = atan2(p2[1], p2[0]);
	dtheta = theta2 - theta1;
   
	while (dtheta > PI)
		dtheta -= 2*PI;
	while (dtheta < -PI)
		dtheta += 2*PI;

   return dtheta;
}

bool IsPointInsidePolygon(DoublePointType pt, DoublePointArray &contour)
{
	double angle=0;
	DoublePointType p1, p2;
	unsigned int pt_num = contour.size();

	for (unsigned int i=0; i<pt_num; i++) 
	{
		for (unsigned int j=0; j<2; j++)
		{
			p1[j] = contour[i][j] - pt[j];
			p2[j] = contour[(i+1)%pt_num][j] - pt[j];
			angle += Angle2D(p1, p2);
		}
   }

   if (fabs(angle) < PI)
	   return false;
   else
	   return true;
}

/*
int radContourOperations::RemoveConeMarker(DoublePointType &deleted_pt, double toleracne_ratio, MarkerInformation & marker_infor)
{
	assert(toleracne_ratio > 0.0);
	double tolerance_error = toleracne_ratio*(marker_infor.split_image->GetSpacing()[0]+marker_infor.split_image->GetSpacing()[1]);
	double min_dist = std::numeric_limits<double>::max(), dist;
	int min_id = -1;
	for (unsigned int i=0; i<marker_infor.cone_point_markers.size(); i++)
	{
		if (marker_infor.cone_point_markers[i].is_visible)
			//&& fabs(marker_infor.cone_marker_shapes[i].shape_center[0]-deleted_pt[0]) < 0.01
			//&& fabs(marker_infor.cone_marker_shapes[i].shape_center[1]-deleted_pt[1]) < 0.01)
		{
			dist = marker_infor.cone_point_markers[i].marker_center.EuclideanDistanceTo(deleted_pt);
			if (dist < min_dist && dist < tolerance_error)
			{
				min_dist = dist;
				min_id = i;
			}
		}
	}

	if (min_id != -1)
		marker_infor.cone_point_markers[min_id].is_visible = false;
	return min_id;

}

void radContourOperations::AddConeMarker(DoublePointType &pt, int detection_category, MarkerInformation &marker_infor)
{
	PointMarker marker_pt;
	marker_pt.Initialize();

	marker_pt.marker_center[0] = pt[0];
	marker_pt.marker_center[1] = pt[1];
	marker_pt.marker_type = detection_category;

	marker_infor.cone_point_markers.push_back(marker_pt);
}
*/

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