/*
 *  radImageView.h
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef radImageView_H
#define radImageView_H

// Math
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "radImgFunc.h"
#include "radVoronoi.h"

#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkShortArray.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkImageActor.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkObjectFactory.h>
#include <QVTKInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkGlyph3D.h>
#include <vtkCellArray.h>
#include <vtkActor.h>
#include <vtkGlyphSource2D.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkImageMapToColors.h>
#include <vtkRegularPolygonSource.h>
#include <vtkPointPicker.h>
#include <vtkObjectFactory.h>
#include <vtkRendererCollection.h>
#include <vtkSphereSource.h>
#include <vtkTubeFilter.h>
#include <vtkContourFilter.h>
#include <vtkTransform.h>
#include <vtkTriangleFilter.h>
#include <vtkContourWidget.h>
#include <vtkOrientedGlyphContourRepresentation.h>
#include <vtkLinearContourLineInterpolator.h>

//#include "radWatershed.h"

// Define interaction style
class radMouseInteractorStylePP : public vtkInteractorStyleImage
{
private:
	DoublePointArray contour_pts;
	bool LeftMousedPressed;
	bool ControlDown, ShiftDown, MouseScroll;
	bool MouseIn;

	int img_dims[3];
	int last_pick_value;
	ColorInfo ci;

	void add_contour_pt(double picked[3]);
	void closest_border(double picked[3], DoublePointType pt, int* pidx = NULL);
	void closest_border_2(double picked[3], DoublePointType pt);

public:
    static radMouseInteractorStylePP* New();
    vtkTypeMacro(radMouseInteractorStylePP, vtkInteractorStyleImage);

	radMouseInteractorStylePP()
	{
		LeftMousedPressed = false;
		ControlDown = false;
		ShiftDown = false;
		MouseScroll = false;
		MouseIn = false;
		last_pick_value = 0;
		img_dims[0] = img_dims[1] = img_dims[2] = 0;
	}
 
	virtual void OnChar() override;
	virtual void OnKeyDown() override;
	virtual void OnKeyUp() override;
	virtual void OnEnter() override;
	virtual void OnLeave() override;
	virtual void OnLeftButtonDown() override;
	virtual void OnMouseMove() override;
	virtual void OnLeftButtonUp() override;
	virtual void OnMiddleButtonDown() override;
	virtual void OnMiddleButtonUp() override;
	virtual void OnRightButtonDown() override;
};

class callbackContourWidget : public vtkCommand
{
public:

    static callbackContourWidget *New()
    {
      return new callbackContourWidget;
    }
    callbackContourWidget(){}

    virtual void Execute(vtkObject *caller, unsigned long, void*);
};


//ground truth color
const double GroundTruthColor[3] = {221.0/255.0, 28.0/255.0, 119.0/255.0};
const double ResultColor[3] = {0.0/255.0, 255.0/255.0, 0.0/255.0};
const double RegionColor[3] = {117.0/255.0, 107.0/255.0, 177.0/255.0};
const double SmallDisplacement = -0.01;

class radImageView
{
private:
	bool interpolationFlag = false;
	bool voronoiFlag = false;
	
	vtkSmartPointer<vtkImageData> ImageData;
	vtkSmartPointer<vtkImageActor> ImageActor;
	void DrawInputImage();

	vtkSmartPointer<vtkPoints> InteractiveContourPoints;
	vtkSmartPointer<vtkCellArray> InteractiveContourCells;
	vtkSmartPointer<vtkPolyData> InteractiveContourPolydata;
	vtkSmartPointer<vtkPolyDataMapper> InteractiveContourMapper;
	vtkSmartPointer<vtkActor> InteractiveContourActor;
	void DrawInteractiveContours();

	vtkSmartPointer<vtkPoints> CenterPoints;
	vtkSmartPointer<vtkPolyData> CenterPolydata;
	vtkSmartPointer<vtkGlyphSource2D> CenterGlyphSource;
	vtkSmartPointer<vtkGlyph3D> CenterGlyph;
	vtkSmartPointer<vtkDataSetMapper> CenterMapper;
	vtkSmartPointer<vtkActor> CenterActor;

	vtkSmartPointer<vtkPoints> ContourPoints, RegionPoints;
	vtkSmartPointer<vtkCellArray> ContourCells, RegionCells;
	vtkSmartPointer<vtkPolyData> ContourPolydata, RegionPolyData;
	vtkSmartPointer<vtkTriangleFilter> RegionTriangleFilter;
	vtkSmartPointer<vtkPolyDataMapper> ContourMapper, RegionMapper;
	vtkSmartPointer<vtkActor> ContourActor, RegionActor;
	void DrawContours();

	vtkSmartPointer<vtkPoints> VoronoiContourPoints;
	vtkSmartPointer<vtkCellArray> VoronoiContourCells;
	vtkSmartPointer<vtkPolyData> VoronoiContourPolydata;
	vtkSmartPointer<vtkPolyDataMapper> VoronoiContourMapper;
	vtkSmartPointer<vtkActor> VoronoiContourActor;
	void DrawVoronoiContours();

	vtkSmartPointer<vtkPoints> EditedContourPoints;
	vtkSmartPointer<vtkCellArray> EditedContourCells;
	vtkSmartPointer<vtkPolyData> EditedContourPoly;
	vtkSmartPointer<vtkContourWidget> EditedContourWidget;
	vtkSmartPointer<vtkOrientedGlyphContourRepresentation> EditedContourRepresentation;
	void DrawEditedContours();

	vtkSmartPointer<vtkRenderer> ImageRender;
	vtkSmartPointer<radMouseInteractorStylePP> ImageStyle;
    vtkSmartPointer<vtkRenderWindowInteractor> WinInteractor;
    vtkSmartPointer<vtkRenderWindow> RenderWin;

	void ChangeCameraOrientation(vtkSmartPointer<vtkRenderer>);

	MarkerSystemSettings *SystemSettings;
	void ApplySystemSettings();

public:
	
	radImageView(MarkerSystemSettings *);
	~radImageView();
    
	void GetImageDimensions(int dims[3]) {
		ImageData->GetDimensions(dims);
	}

	void ResetView(bool camera_flag = true); //used for initialization

    vtkSmartPointer<vtkRenderWindow> GetRenderWin() {return RenderWin;}

	void SetSplitImage(RGBImageType::Pointer);
	void SetColorInfo(ColorInfo ci);
	void SetColorInfo(double color_level, double color_window);
	ColorInfo GetColorInfo();

	void InitializeView();
	void SetConeContourVisibility(bool flag);
	void SetConeCenterVisibility(bool flag);
	void SetConeRegionVisibility(bool flag);
	void SetConeRegionOpacity(double op);
	void SetCenterGlyphScale(double scale);

	bool GetInterpolation() { return interpolationFlag; }
	void SetInterpolation(bool flag);
	
	void SetConeContourColor(double rgb[3]);
	void SetConeCenterColor(double rgb[3]);
	void SetConeRegionColor(double rgb[3]);

	void SetContourWidth(int);
	void SetInteractiveContours(DoublePointArray &, bool ending_flag = false);
	void SetContourMarkers(vector< ContourMarker > &);
	void GetContourPoints(DoublePointArray& markers);

	bool getVoronoi() { return voronoiFlag; }
	void setVoronoi(bool flag);

	void updateVoronoiDiagram();

	void EnableEditedContour(DoublePointArray &);
	void DisableEditedContour();

	void SetVoronoiContours(std::vector< DoublePointArray>& contours);
};

#endif // radImageView_H

