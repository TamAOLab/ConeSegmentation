/*
 *  radBackup.cpp
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "radBackup.h"
#include <QDebug>
#include "QuickView.h"

radBackup::radBackup()
{
	
}

radBackup::~radBackup()
{ 
}

template <class T>
void radBackup::ExtractNumsFromString(string &record, vector<T> &arr, int data_type)
{
	//record.erase(remove_if(record.begin(), record.end(), isspace), record.end());

	arr.clear();
	if (record.length() == 0)
	{
		return;
	}

	size_t found, found1;
    string value_str;
    found = record.find(',');
    found1 = 0;
	
    while (found != string::npos) 
	{
        if (found1 == 0) 
		{
			value_str.assign(record.begin()+found1, record.begin()+found);
        }
        else
            value_str.assign(record.begin()+found1+1, record.begin()+found);
        
		if (data_type == 0) //int/short
			arr.push_back(atoi(value_str.c_str()));
		else if (data_type == 1) //float/double
			arr.push_back(atof(value_str.c_str()));
		else if (data_type == 2) //unsigned int
		{
			arr.push_back(atol(value_str.c_str()));
		}

		found1 = found;
		found = record.find(',', found+1);
    }
    
    value_str.assign(record.begin()+found1+1, record.end());
    arr.push_back(atof(value_str.c_str()));
}

bool radBackup::ReadTxtFile(const char *filename, vector<string> & str_list)
{
	string line;
	ifstream myfile(filename);

	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			if (!line.empty())
				str_list.push_back(line);
		}
		myfile.close();
		return true;
	}
	else 
		return false;
}

bool radBackup::IsValidHistoryFile(vector<string> & str_list)
{
	if (str_list.size() != 3)
		return false;

	if (str_list[0].find(".tif") == string::npos)
		return false;

	if (str_list[2].find(".txt") == string::npos)
		return false;

	return true;
}

void radBackup::ReadContourMarkers(string & file_name, vector< ContourMarker > & cone_shapes, itk::ImageRegion<2U> rgn)
{
	cone_shapes.clear();

	string line;
	ifstream myfile(file_name.c_str());
	vector<double> tmp_array;
	bool shape_valid_flag = true;
	ContourMarker tmp_shape;
	int loop_id = 0;

	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			if (!line.empty())
			{
				if (line.find("contours") != string::npos)
					continue;

				if (loop_id % 3 == 0)
				{
					ExtractNumsFromString<double>(line, tmp_array, 1);
					if (tmp_array.size() == 2)
					{
						tmp_shape.marker_center[0] = tmp_array[0];
						tmp_shape.marker_center[1] = tmp_array[1];
					}
					else
						shape_valid_flag = false;
				}
				else if (loop_id % 3 == 1)
				{
					ExtractNumsFromString<double>(line, tmp_array, 1);
					if (tmp_array.size() == 4)
					{
						tmp_shape.marker_bounding_pts[0][0] = tmp_array[0];
						tmp_shape.marker_bounding_pts[0][1] = tmp_array[1];
						tmp_shape.marker_bounding_pts[1][0] = tmp_array[2];
						tmp_shape.marker_bounding_pts[1][1] = tmp_array[3];
						// skip shapes that cross image boundaries
						shape_valid_flag = tmp_shape.withinImageRegion(rgn);
					}
					else
						shape_valid_flag = false;
				}
				else if (loop_id % 3 == 2)
				{
					ExtractNumsFromString<double>(line, tmp_array, 1);
					if (tmp_array.size() < 6)
					{
						shape_valid_flag = false;
					}
					else
					{
						DoublePointType pt;
						for (unsigned int i = 0; i<tmp_array.size() / 2; i++)
						{
							pt[0] = tmp_array[2 * i];
							pt[1] = tmp_array[2 * i + 1];
							tmp_shape.marker_contours.push_back(pt);
						}
					}

					if (shape_valid_flag)
					{
						cone_shapes.push_back(tmp_shape);
						tmp_shape.Initialize();
					}

					shape_valid_flag = true;
				}
				loop_id++;
			}
		}
	}
}

void radBackup::ReadBackup(MarkerInformation & split_infor)
{
	//get all sub-directories
	QDirIterator it(BackupDir.c_str(), QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot, QDirIterator::Subdirectories);
	vector<string> str_list;
	string img_file_path;
	
	while (it.hasNext()) 
	{
		str_list.clear();
		string file_dir = it.next().toStdString();
		string filename = file_dir+"/"+split_infor.split_file_names.second+".txt";

		if (ReadTxtFile(filename.c_str(), str_list) && IsValidHistoryFile(str_list))
		{
			if (str_list[0].compare(split_infor.split_file_names.first) != 0)
				continue;

			ReadSegmentationParameters(file_dir, split_infor);

			//load existing marker shapes
			ReadContourMarkers(str_list[2], split_infor.cone_contour_markers, split_infor.split_image->GetLargestPossibleRegion());
			for (unsigned int i=0; i<split_infor.cone_contour_markers.size(); i++)
			{
				split_infor.cone_contour_markers[i].ComputeShapeRegion<RGBImageType>(split_infor.split_image);
			}
			
		}
	}
}

void radBackup::WriteContourMarkers(string & file_name, vector< ContourMarker > & cone_shapes)
{
	ofstream myfile;

	myfile.open(file_name.c_str());
	myfile << "contours" << std::endl;
	for (unsigned int i = 0; i<cone_shapes.size(); i++)
	{
		if (cone_shapes[i].is_visible)
		{
			myfile << cone_shapes[i].marker_center[0] << ", "
				<< cone_shapes[i].marker_center[1] << std::endl;

			myfile << cone_shapes[i].marker_bounding_pts[0][0] << ", "
				<< cone_shapes[i].marker_bounding_pts[0][1] << ", " 
				<< cone_shapes[i].marker_bounding_pts[1][0] << ", "
				<< cone_shapes[i].marker_bounding_pts[1][1] << std::endl;

			for (unsigned int j = 0; j<cone_shapes[i].marker_contours.size() - 1; j++)
				myfile << cone_shapes[i].marker_contours[j][0] << ", " << cone_shapes[i].marker_contours[j][1] << ", ";
			unsigned int lid = cone_shapes[i].marker_contours.size() - 1;
			myfile << cone_shapes[i].marker_contours[lid][0] << ", " << cone_shapes[i].marker_contours[lid][1] << std::endl;
		}
	}

	myfile.close();
}

bool radBackup::RemoveDir(const QString &dirName)
{
	bool result = true;
    QDir dir(dirName);
 
    if (dir.exists(dirName)) 
	{
        Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst)) 
		{
            if (info.isDir()) 
			{
                result = RemoveDir(info.absoluteFilePath());
            }
            else 
			{
                result = QFile::remove(info.absoluteFilePath());
            }
 
            if (!result) 
			{
                return result;
            }
        }
        result = dir.rmdir(dirName);
    }
 
    return result;
}

bool radBackup::IsRecordExisted(string & file_name, string & searched_str)
{
	string line;
	ifstream myfile(file_name.c_str());

	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			if (line.compare(searched_str) == 0)
			{
				myfile.close();
				return true;
			}
		}
		myfile.close();
		return false;
	}
	else 
		return false;
}

bool radBackup::IsDiretoryExisted(pair<string, string> & dir_names)
{
	QDir cur_dir(BackupDir.c_str());
	cur_dir.setFilter(QDir::AllDirs | QDir::NoDotAndDotDot | QDir::NoSymLinks);
    QStringList dirList = cur_dir.entryList();

	for (int i=0; i<dirList.size(); i++)
	{
		if (dirList[i].toStdString().find(dir_names.second) != string::npos)
		{
			QString newPath = QString("%1/%2").arg(cur_dir.absolutePath()).arg(dirList.at(i));
			string filename = newPath.toStdString()+"/"+dir_names.second+".txt";
			if (IsRecordExisted(filename, dir_names.first))
				return true;
		}
	}

	return false;
}

int radBackup::GetDirectoryIndex(pair<string, string> & dir_names)
{
	QDir cur_dir(BackupDir.c_str());
	cur_dir.setFilter(QDir::AllDirs | QDir::NoDotAndDotDot | QDir::NoSymLinks);
    QStringList dirList = cur_dir.entryList();

	for (int i=0; i<dirList.size(); i++)
	{
		if (dirList[i].toStdString().find(dir_names.second) != string::npos)
		{
			QString newPath = QString("%1/%2").arg(cur_dir.absolutePath()).arg(dirList.at(i));
			string filename = newPath.toStdString()+"/"+dir_names.second+".txt";
			if (IsRecordExisted(filename, dir_names.first))
			{
				if (dirList[i].toStdString().length() == dir_names.second.length())
					return 0;
				else
				{
					string tmp_str = dirList[i].toStdString();
					tmp_str.assign(tmp_str.begin()+dir_names.second.length(), tmp_str.end());
					return atoi(tmp_str.c_str());
				}
			}
		}
	}

	return -1;
}

int radBackup::GetMaximumDirectoryIndex(pair<string, string> & dir_names)
{
	QDir cur_dir(BackupDir.c_str());
	cur_dir.setFilter(QDir::AllDirs | QDir::NoDotAndDotDot | QDir::NoSymLinks);
    QStringList dirList = cur_dir.entryList();
	int result_id = -1;

	for (int i=0; i<dirList.size(); i++)
	{
		if (dirList[i].toStdString().find(dir_names.second) != string::npos)
		{
			if (dirList[i].toStdString().length() != dir_names.second.length())
			{
				string tmp_str = dirList[i].toStdString();
				tmp_str.assign(tmp_str.begin()+dir_names.second.length(), tmp_str.end());
				if (result_id < atoi(tmp_str.c_str()))
					result_id = atoi(tmp_str.c_str());
			}
			else
			{
				if (result_id < 0)
					result_id = 0;
			}
		}
	}

	return result_id;
}

string radBackup::WriteSegmentationParameters(string & file_dir, MarkerInformation & split_infor)
{
	ConeSegmentationParameters & params = split_infor.segment_params;
	QJsonObject jobj;
	jobj["HessianThreshold"] = QJsonValue(params.HessianThreshold);
	jobj["GACIterationNumber"] = QJsonValue(params.GACIterationNumber);

	string filename = file_dir + "/" + split_infor.split_file_names.second + "_param.json";
	QFile fo(filename.c_str());
	if (fo.open(QIODevice::WriteOnly)) {
		QJsonDocument json = QJsonDocument(jobj);
		fo.write(json.toJson(QJsonDocument::Indented));
		fo.close();
	}
	return filename;
}

void radBackup::ReadSegmentationParameters(string &file_dir, MarkerInformation  & split_infor)
{
	string filename = file_dir + "/" + split_infor.split_file_names.second + "_param.json";
	QFile fi(filename.c_str());
	if (!fi.open(QIODevice::ReadOnly))
		return;

	QJsonDocument json = QJsonDocument::fromJson(fi.readAll());
	fi.close();

	QJsonObject jobj = json.object();
	ConeSegmentationParameters & params = split_infor.segment_params;

	params.loadDefaults();
	if (jobj["HessianThreshold"].isDouble())
		params.HessianThreshold = jobj["HessianThreshold"].toDouble();
	if (jobj["GACIterationNumber"].isDouble())
		params.GACIterationNumber = jobj["GACIterationNumber"].toInt();
}


void radBackup::WriteBackup(MarkerInformation & split_infor)
{
	string filename, filename1;
	string file_dir;

	//new version - to handle the existence of files with duplicated names
	if (IsDiretoryExisted(split_infor.split_file_names))
	{
		int dir_id = GetDirectoryIndex(split_infor.split_file_names);
		if (dir_id <= 0)
		{
			file_dir = BackupDir + "/" + split_infor.split_file_names.second;
		}
		else
			file_dir = BackupDir + "/" + split_infor.split_file_names.second+QString::number(dir_id).toStdString();
	}
	else
	{
		int dir_id = GetMaximumDirectoryIndex(split_infor.split_file_names);
		if (dir_id < 0)
		{
			file_dir = BackupDir + "/" + split_infor.split_file_names.second;
		}
		else
			file_dir = BackupDir + "/" + split_infor.split_file_names.second+QString::number(dir_id+1).toStdString();

		QDir infor_dir(file_dir.c_str());
		infor_dir.mkpath(".");
	}	

	ofstream myfile;

	filename.assign(file_dir+"/"+split_infor.split_file_names.second+".txt");
	myfile.open(filename.c_str());
	myfile << split_infor.split_file_names.first << std::endl;

	filename1 = WriteSegmentationParameters(file_dir, split_infor);
	myfile << WriteSegmentationParameters(file_dir, split_infor) << std::endl;
	WriteSegmentationParameters(file_dir, split_infor);

	filename1.assign(file_dir + "/" + split_infor.split_file_names.second + "_contours.txt");
	myfile << filename1 << std::endl;
	WriteContourMarkers(filename1, split_infor.cone_contour_markers);

	myfile.close();
}