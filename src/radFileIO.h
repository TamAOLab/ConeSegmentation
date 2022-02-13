/*
 *  radFileIO.h
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef radFileIO_H
#define radFileIO_H

#define _CRT_SECURE_NO_WARNINGS

#include "radImgFunc.h"
#include <itkSpatialObjectToImageFilter.h>
#include <itkPasteImageFilter.h>
#include <itkJoinSeriesImageFilter.h>
#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QProgressBar>
#include <QApplication>
#include <QImageReader>
#include <QImage>
#include <QJsonObject>
#include <QJsonArray>
#include <QJsonDocument>

// Math
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

void ContourMarkersToJSON(MarkerInformation & split_infor, QFile & fo);
int ContourMarkersFromJSON(QFile & fi, std::vector<MarkerInformation> & SplitMarkerInfor, bool ignore_path=false);

void SaveContourMarkers(MarkerInformation & split_infor, QString jsonfile);

class radFileIO
{
private:
	
public:
	
	radFileIO();
	~radFileIO();
    
	bool IsFileExisted(string, vector< MarkerInformation > &);
	RGBImageType::Pointer ReadSplitImage(string);
};

#endif // segFileIO_H

