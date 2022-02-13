/*
 *  radFileIO.cpp
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "radFileIO.h"
#include <itkPoint.h>
#include <itkPolygonSpatialObject.h>
#include <vtkPointData.h>
#include <QDebug>
#include "QuickView.h"
#include "itkRGBToLuminanceImageFilter.h"

static QJsonArray pt2json(DoublePointType & pt)
{
	QJsonArray res;
	res.append(QJsonValue(pt[0]));
	res.append(QJsonValue(pt[1]));
	return res;
}

void ContourMarkersToJSON(MarkerInformation & split_infor, QFile & fo)
{
	vector<ContourMarker> & markers = split_infor.cone_contour_markers;
	QJsonArray jMarkers;
	for (ContourMarker & marker : markers)
	{
		if (!marker.is_visible) continue;
		QJsonObject jMarker;
		jMarker["center"] = pt2json(marker.marker_center);
		QJsonArray jBoundingPts;
		jBoundingPts.append(pt2json(marker.marker_bounding_pts[0]));
		jBoundingPts.append(pt2json(marker.marker_bounding_pts[1]));
		jMarker["bounding_pts"] = jBoundingPts;
		QJsonArray jContours;
		for (auto & pt : marker.marker_contours)
			jContours.append(pt2json(pt));
		jMarker["contours"] = jContours;
		jMarkers.append(jMarker);
	}
	QJsonObject jobj;
	jobj["filename"] = split_infor.split_file_names.first.c_str();
	jobj["type"] = "MarkerContours";
	jobj["markers"] = jMarkers;
	QJsonDocument json = QJsonDocument(jobj);
	if (fo.open(QIODevice::WriteOnly)) {
		fo.write(json.toJson(QJsonDocument::Compact));
		fo.close();
	}
}

static DoublePointType json2pt(QJsonArray jpt)
{
	DoublePointType pt;
	pt[0] = jpt[0].toDouble();
	pt[1] = jpt[1].toDouble();
	return pt;
}

static bool jsonToContourMarker(QJsonObject & jobj, MarkerInformation & split_infor)
{
	if (!jobj["type"].isString()) return false;
	QString jType = jobj["type"].toString();
	if (jType != "MarkerContours") {
		std::cout << "Wrong marker type: " << jType.toStdString() << std::endl;
		return false;
	}
	if (!jobj["markers"].isArray()) return false;
	split_infor.ClearStack(split_infor.contour_operator_stack);
	QJsonArray jMarkers = jobj["markers"].toArray();
	split_infor.cone_contour_markers.resize((size_t)jMarkers.size());
	size_t iiMark = 0;
	itk::ImageRegion<2U> rgn = split_infor.split_image->GetLargestPossibleRegion();
	for (int iMark = 0; iMark < jMarkers.size(); iMark++) {
		QJsonObject jMarker = jMarkers[iMark].toObject();
		ContourMarker & marker = split_infor.cone_contour_markers[iiMark];
		marker.is_visible = true;
		marker.marker_center = json2pt(jMarker["center"].toArray());
		QJsonArray jBoundingPts = jMarker["bounding_pts"].toArray();
		marker.marker_bounding_pts[0] = json2pt(jBoundingPts[0].toArray());
		marker.marker_bounding_pts[1] = json2pt(jBoundingPts[1].toArray());
		QJsonArray jContours = jMarker["contours"].toArray();
		if (!marker.withinImageRegion(rgn))
			// skip markers that cross image boundaries
			continue;
		marker.marker_contours.resize((size_t)jContours.size());
		for (int ic = 0; ic < jContours.size(); ic++) {
			marker.marker_contours[(size_t)ic] = json2pt(jContours[ic].toArray());
		}
		iiMark++;
	}
	split_infor.cone_contour_markers.resize(iiMark);
	return true;
}

int ContourMarkersFromJSON(QFile & fi, std::vector<MarkerInformation> & SplitMarkerInfor, bool ignore_path)
{
	if (fi.open(QIODevice::ReadOnly)) {
		QJsonDocument json = QJsonDocument::fromJson(fi.readAll());
		fi.close();

		QJsonObject jobj = json.object();
		if (!jobj["filename"].isString()) return -1;
		std::string filename = jobj["filename"].toString().toStdString();
		if (ignore_path) {
			QFileInfo qfp(filename.c_str());
			filename = qfp.completeBaseName().toStdString();
		}
		for (size_t id = 0; id < SplitMarkerInfor.size(); id++) {
			if (ignore_path && SplitMarkerInfor[id].split_file_names.second != filename) continue;
			if (!ignore_path && SplitMarkerInfor[id].split_file_names.first != filename) continue;
			if (jsonToContourMarker(jobj, SplitMarkerInfor[id]))
				return (int)id;
			else
				return -1;
		}
	}
	return -1;
}

static double ComputeContourDiameter(ContourMarker & marker, double *p_area)
{
	typedef itk::PolygonSpatialObject<3> PolygonType;

	PolygonType::Pointer polygon = PolygonType::New();
	PolygonType::PointType point;

	for (int i = 0; i < marker.marker_contours.size(); i++)
	{
		point[0] = marker.marker_contours[i][0];
		point[1] = marker.marker_contours[i][1];
		point[2] = 0;
		polygon->AddPoint(point);
	}
	double area = polygon->MeasureArea();
	double diameter = sqrt(area / PI) * 2;
	if (p_area) *p_area = area;
	return diameter;
}

void SaveContourMarkers(MarkerInformation & split_infor, QString jsonfile)
{
	// Save JSON
    QFile qjf(jsonfile);
	ContourMarkersToJSON(split_infor, qjf);

	vector<ContourMarker> & markers = split_infor.cone_contour_markers;

	// Generate CSV file names
	QString baseName = jsonfile;
	if (baseName.endsWith(QString(".json"), Qt::CaseInsensitive))
		baseName.remove(baseName.length() - 5, 5);
	QString contoursName = baseName + "_contours.csv";
	QString centersName = baseName + "_detections.csv";
	QString measurementsName = baseName + "_measurements.csv";
	// Save contours
#ifdef _WIN32
	ofstream foutContours(contoursName.toStdWString(), std::ofstream::out);
#else
	ofstream foutContours(contoursName.toStdString(), std::ofstream::out);
#endif
	if (foutContours.is_open()) {
		foutContours << "# Contour point coordinate pairs: X1 Y1 X2 Y2 X3 Y3 ..." << std::endl;
		for (ContourMarker & marker : markers)
		{
			if (!marker.is_visible) continue;
			size_t ncont = marker.marker_contours.size();
			for (size_t j = 0; j < ncont; j++) {
				foutContours << marker.marker_contours[j][0] << "," << marker.marker_contours[j][1];
				if (j == ncont - 1)
					foutContours << std::endl;
				else
					foutContours << ",";
			}
		}
		foutContours.close();
	}
	// Save center points
#ifdef _WIN32
	ofstream foutCenters(centersName.toStdWString(), std::ofstream::out);
#else
	ofstream foutCenters(centersName.toStdString(), std::ofstream::out);
#endif
	if (foutCenters.is_open()) {
		foutCenters << "# Contour center point coordinates: XCenter YCenter" << std::endl;
		for (ContourMarker & marker : markers)
		{
			if (!marker.is_visible) continue;
			foutCenters << marker.marker_center[0] << "," << marker.marker_center[1] << std::endl;
		}
		foutCenters.close();
	}
	// Save areas and diameters
#ifdef _WIN32
	ofstream foutMeasurements(measurementsName.toStdWString(), std::ofstream::out);
#else
	ofstream foutMeasurements(measurementsName.toStdString(), std::ofstream::out);
#endif
	if (foutMeasurements.is_open()) {
		foutMeasurements << "# Contour measurements: Area Diameter" << std::endl;
		for (ContourMarker & marker : markers)
		{
			if (!marker.is_visible) continue;
			double area;
			double diameter = ComputeContourDiameter(marker, &area);
			foutMeasurements << area << "," << diameter << std::endl;
		}
		foutMeasurements.close();
	}
}


radFileIO::radFileIO()
{
}

radFileIO::~radFileIO()
{ 
}

bool radFileIO::IsFileExisted(string full_file_name, vector< MarkerInformation > & split_infor)
{
	for (unsigned int i=0; i<split_infor.size(); i++)
	{
		if (split_infor[i].split_file_names.first.compare(full_file_name) == 0)
			return true;
	}

	return false;
}

RGBImageType::Pointer radFileIO::ReadSplitImage(string fileName)
{
	typedef itk::ImageFileReader<RGBImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	RGBImageType::Pointer res = NULL;

	reader->SetFileName(fileName.c_str());
	try {
		reader->Update();
		res = reader->GetOutput();
	}
	catch (...) {}

	return res;
}

