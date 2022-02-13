/*
 *  radBackup.h
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef radBackup_H
#define radBackup_H

#include "radImgFunc.h"
#include <QJsonObject>
#include <QJsonDocument>
#include <QDir>
#include <QDirIterator>

// Math
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>

class radBackup
{
private:

	template <class T>
	void ExtractNumsFromString(string &, vector<T> &, int data_type);

	bool IsDiretoryExisted(pair<string, string> &);
	bool IsRecordExisted(string &, string &);
	int GetDirectoryIndex(pair<string, string> &);
	int GetMaximumDirectoryIndex(pair<string, string> &);
	string BackupDir;

	void WriteContourMarkers(string &, vector< ContourMarker > &);
	void ReadContourMarkers(string &, vector< ContourMarker > &, itk::ImageRegion<2U>);
	bool ReadTxtFile(const char *, vector<string> &);
	bool IsValidHistoryFile(vector<string> &);
	
public:
	
	radBackup();
	~radBackup();

	string WriteSegmentationParameters(string &file_dir, MarkerInformation  & split_infor);
	void ReadSegmentationParameters(string &file_dir, MarkerInformation  & split_infor);

	void WriteBackup(MarkerInformation &);
	bool ReadBackup(MarkerInformation &);

	bool RemoveDir(const QString &dirName);

	inline void SetBackupDir(string str) {BackupDir.assign(str);}
};

#endif // radBackup_H

