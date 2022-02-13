/*
 *  radImageView.cpp
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <QApplication>
#include <QCursor>
#include "radImageView.h"
#include "radmainwindow.h"
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkCellData.h>
#include <vtkTextProperty.h>
#include <vtkImageMapper3D.h>
#include <vtkImageProperty.h>
#include <vtkNamedColors.h>

vtkStandardNewMacro(radMouseInteractorStylePP);

void radMouseInteractorStylePP::add_contour_pt(double picked[3])
{
	DoublePointType pt;
	pt[0] = picked[0];
	pt[1] = picked[1];
	contour_pts.push_back(pt);
	radMainWindow::GetPointer()->GetImageView()->SetInteractiveContours(contour_pts);
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radMouseInteractorStylePP::closest_border(double picked[3], DoublePointType pt, int* pidx)
{
	double dx0 = fabs(pt[0]);
	double dx1 = fabs(pt[0] - img_dims[0]);
	double dy0 = fabs(pt[1]);
	double dy1 = fabs(pt[1] - img_dims[1]);
	int idx = 0;
	if (dx0 <= dx1 && dx0 <= dy0 && dx0 <= dy1) {
		picked[0] = 0.;
		picked[1] = pt[1];
	}
	else if (dx1 <= dy0 && dx1 <= dy1) {
		picked[0] = double(img_dims[0]);
		picked[1] = pt[1];
	}
	else if (dy0 <= dy1) {
		picked[0] = pt[0];
		picked[1] = 0.;
		idx = 1;
	}
	else {
		picked[0] = pt[0];
		picked[1] = double(img_dims[1]);
		idx = 1;
	}
	if (pidx) *pidx = idx;
	// std::cout << "closest_border(): " << picked[0] << " " << picked[1] << std::endl;
}

void radMouseInteractorStylePP::closest_border_2(double picked[3], DoublePointType pt)
{
	int idx1;
	closest_border(picked, pt, &idx1);
	if (!contour_pts.empty()) {
		int idx0;
		double corner[3];
		closest_border(corner, contour_pts[contour_pts.size() - 1], &idx0);
		if (idx0 != idx1) {
			pt[idx0] = corner[idx0];
			pt[idx1] = picked[idx1];
			// std::cout << "corner pt: " << pt[0] << " " << pt[1] << std::endl;
			contour_pts.push_back(pt);
		}
	}
}

void radMouseInteractorStylePP::OnChar()
{
	vtkRenderWindowInteractor* rwi = this->Interactor;
	char c = rwi->GetKeyCode();
	if (!strchr("wWrRxXyYzZfF", c)) {
		vtkInteractorStyleImage::OnChar();
	}
}

void radMouseInteractorStylePP::OnKeyDown()
{
	while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
	if (!MouseIn) {
		vtkInteractorStyleImage::OnKeyDown();
		return;
	}
	vtkRenderWindowInteractor *rwi = this->Interactor;
	std::string key = rwi->GetKeySym();
	if (key == "Up") {
		radMainWindow::GetPointer()->PreviousImage();
		return;
	}
	else if (key == "Down") {
		radMainWindow::GetPointer()->NextImage();
		return;
	}
	if (key == "Shift_L") {
		ShiftDown = true;
		if (!MouseScroll) {
			QApplication::setOverrideCursor(Qt::OpenHandCursor);
		}
	}
	else {
		MouseOperation op = radMainWindow::GetPointer()->MouseOperationType;
		if (op == Mouse_Add_Contour || op == Mouse_Edit_Contour)
		{
			if (key == "Control_L") {
				ControlDown = true;
				while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
				QApplication::setOverrideCursor(Qt::CrossCursor);
				return;
			}
		}
	}
	vtkInteractorStyleImage::OnKeyDown();
}

void radMouseInteractorStylePP::OnKeyUp()
{
	while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
	vtkRenderWindowInteractor *rwi = this->Interactor;
	std::string key = rwi->GetKeySym();
	if (key == "Up" || key == "Down") {
		return;
	}
	if (ControlDown) {
		if (key == "Control_L") {
			ControlDown = false;
			return;
		}
	}
	if (key == "Shift_L") {
		ShiftDown = false;
	}
	vtkInteractorStyleImage::OnKeyUp();
}

void radMouseInteractorStylePP::OnEnter()
{
	while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
	MouseOperation op = radMainWindow::GetPointer()->MouseOperationType;
	if (!MouseIn) {
		MouseIn = true;
		if (ShiftDown) QApplication::setOverrideCursor(Qt::OpenHandCursor);
		else if (ControlDown && (op == Mouse_Add_Contour || op == Mouse_Edit_Contour)) {
			QApplication::setOverrideCursor(Qt::CrossCursor);
			return;
		}
	}
	else {
		// If a second OnEnter() received without a matching OnLeave(), the QVTKWidget
		// does not have keyboard focus and it won't receive OnKeyUp() for Ctrl/Shift either.
		// Like the user tried to drag the mouse from another widget while holding Shift down.
		ControlDown = ShiftDown = false;
	}
	vtkInteractorStyleImage::OnEnter();
}

void radMouseInteractorStylePP::OnLeave()
{
	MouseIn = false;
	while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
	vtkInteractorStyleImage::OnLeave();
}

void radMouseInteractorStylePP::OnLeftButtonDown()
{
	while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
	if (!MouseIn) {
		ControlDown = ShiftDown = false;
		vtkInteractorStyleImage::OnLeftButtonDown();
		return;
	}
	if (ShiftDown) {
		MouseScroll = true;
		QApplication::setOverrideCursor(Qt::SizeAllCursor);
		vtkInteractorStyleImage::OnLeftButtonDown();
		return;
	}
	LeftMousedPressed = true;
	MouseOperation op = radMainWindow::GetPointer()->MouseOperationType;

	if (this->Interactor->GetControlKey() && (op == Mouse_Add_Contour || op == Mouse_Edit_Contour)) {
		// Override Add/Edit if Ctrl is down
		op = Mouse_Delete_Single_Contour;
		ControlDown = true;
		QApplication::setOverrideCursor(Qt::CrossCursor);
	}
	if (op == Mouse_Add_Contour || op == Mouse_Delete_Contour)
	{
		int pick_value =  this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0], 
			this->Interactor->GetEventPosition()[1], 0,  // always zero.
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);
		
		//std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << ", " << pick_value << std::endl;
		radMainWindow::GetPointer()->GetImageView()->GetImageDimensions(img_dims);
		if (pick_value != 0)
		{
			add_contour_pt(picked);
		}
		last_pick_value = pick_value;
	}
	else if (op == Mouse_Edit_Contour)
	{
		int pick_value =  this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0], 
			this->Interactor->GetEventPosition()[1], 0,  // always zero.
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);
		
		//std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << ", " << pick_value << std::endl;
		if (pick_value != 0)
		{
			radMainWindow::GetPointer()->EditConeContours(picked[0], picked[1], picked[2]);
		}
	}
	else if (op == Mouse_Delete_Single_Contour)
	{
		int pick_value = this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
			this->Interactor->GetEventPosition()[1], 0,  // always zero.
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);

		//std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << ", " << pick_value << std::endl;
		if (pick_value != 0)
		{
			radMainWindow::GetPointer()->RemoveSingleConeContour(picked[0], picked[1], picked[2]);
			radMainWindow::GetPointer()->GetImageView()->ResetView(false);
		}
	} else
		vtkInteractorStyleImage::OnLeftButtonDown();
}

void radMouseInteractorStylePP::OnMouseMove()
{
	if (ShiftDown) {
		vtkInteractorStyleImage::OnMouseMove();
		return;
	}

	if ((radMainWindow::GetPointer()->MouseOperationType == Mouse_Add_Contour 
		|| radMainWindow::GetPointer()->MouseOperationType == Mouse_Delete_Contour) 
		&& LeftMousedPressed)
	{
		int pick_value =  this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0], 
			this->Interactor->GetEventPosition()[1], 0,  // always zero.
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);
		
		//std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << ", " << pick_value << std::endl;
		if (pick_value == 0) {
			if (last_pick_value != 0 && !contour_pts.empty()) {
				closest_border(picked, contour_pts[contour_pts.size() - 1]);
				add_contour_pt(picked);
			}
		}
		else {
			if (last_pick_value == 0) {
				DoublePointType pt;
				pt[0] = picked[0];
				pt[1] = picked[1];
				double bord[3];
				closest_border_2(bord, pt);
				add_contour_pt(bord);
			}
			add_contour_pt(picked);
		}
		last_pick_value = pick_value;
		return;
	}
	vtkInteractorStyleImage::OnMouseMove();
}

void radMouseInteractorStylePP::OnLeftButtonUp()
{
	while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
	if (MouseScroll) {
		MouseScroll = false;
		if (ShiftDown) QApplication::setOverrideCursor(Qt::OpenHandCursor);
		vtkInteractorStyleImage::OnLeftButtonUp();
		LeftMousedPressed = false;
		return;
	}
	MouseOperation op = radMainWindow::GetPointer()->MouseOperationType;
	if (this->Interactor->GetControlKey() && (op == Mouse_Add_Contour || op == Mouse_Edit_Contour)) {
		// Override Add/Edit if Ctrl is down
		op = Mouse_Delete_Single_Contour;
		ControlDown = true;
		QApplication::setOverrideCursor(Qt::CrossCursor);
		LeftMousedPressed = false;
		return;
	}
	if ((op == Mouse_Add_Contour || op == Mouse_Delete_Contour)
		&& LeftMousedPressed && !ControlDown)
	{
		int pick_value =  this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0], 
			this->Interactor->GetEventPosition()[1], 0,  // always zero.
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);
		
		//std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << ", " << pick_value << std::endl;
		if (pick_value == 0 && !contour_pts.empty()) {
			closest_border_2(picked, contour_pts[0]);
			pick_value = 1;
		}
		if (pick_value != 0)
		{
			DoublePointType pt;
			pt[0] = picked[0];
			pt[1] = picked[1];
			contour_pts.push_back(pt);
			radMainWindow::GetPointer()->GetImageView()->SetInteractiveContours(contour_pts, true);

			if (radMainWindow::GetPointer()->MouseOperationType == Mouse_Add_Contour)
				radMainWindow::GetPointer()->AddConeContours(contour_pts);
			else if (radMainWindow::GetPointer()->MouseOperationType == Mouse_Delete_Contour)
			{
				radMainWindow::GetPointer()->RemoveConeContours(contour_pts);
			}

			radMainWindow::GetPointer()->GetImageView()->ResetView(false);
		}

		contour_pts.clear();
	} else
		vtkInteractorStyleImage::OnLeftButtonUp();
	LeftMousedPressed = false;
}

void radMouseInteractorStylePP::OnMiddleButtonDown()
{
	while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
	if (MouseIn) {
		MouseScroll = true;
		QApplication::setOverrideCursor(Qt::SizeAllCursor);
	}
	vtkInteractorStyleImage::OnMiddleButtonDown();
}
void radMouseInteractorStylePP::OnMiddleButtonUp()
{
	while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
	MouseScroll = false;
	if (ShiftDown) QApplication::setOverrideCursor(Qt::OpenHandCursor);
	vtkInteractorStyleImage::OnMiddleButtonUp();
}

void radMouseInteractorStylePP::OnRightButtonDown()
{
	if (this->Interactor->GetAltKey())
	{
		int pick_value = this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
			this->Interactor->GetEventPosition()[1], 0,  // always zero.
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
		double picked[3];
		this->Interactor->GetPicker()->GetPickPosition(picked);

		//std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << ", " << pick_value << std::endl;
		if (pick_value != 0)
		{
			radMainWindow::GetPointer()->RemoveSingleConeContour(picked[0], picked[1], picked[2]);
			radMainWindow::GetPointer()->GetImageView()->ResetView(false);
		}
		return;
	}
	vtkInteractorStyleImage::OnRightButtonDown();
}

void callbackContourWidget::Execute(vtkObject *caller, unsigned long, void*)
{
	vtkContourWidget *contourWidget = reinterpret_cast<vtkContourWidget*>(caller);
	vtkContourRepresentation* rep = static_cast<vtkContourRepresentation*>(contourWidget->GetRepresentation());
		
	vtkPolyData *polyData = rep->GetContourRepresentationAsPolyData();
	radMainWindow::GetPointer()->UpdateConeContours(polyData);
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radImageView::DrawInputImage()
{
	interpolationFlag = true;
	ImageData = vtkSmartPointer<vtkImageData>::New();
	ImageData->SetDimensions(1, 1, 1);
	ImageData->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
	ImageActor = vtkSmartPointer<vtkImageActor>::New();
	ImageActor->SetInputData(ImageData);
	if (interpolationFlag)
		ImageActor->InterpolateOn();
	else
		ImageActor->InterpolateOff();
}

void radImageView::DrawInteractiveContours()
{
	//==========================used for contour selection=========================================
	InteractiveContourPoints = vtkSmartPointer<vtkPoints>::New();
	InteractiveContourCells = vtkSmartPointer<vtkCellArray>::New();
	InteractiveContourPolydata = vtkSmartPointer<vtkPolyData>::New();
	InteractiveContourPolydata->SetPoints(InteractiveContourPoints);
	InteractiveContourPolydata->SetLines(InteractiveContourCells);
	InteractiveContourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	InteractiveContourMapper->SetInputData(InteractiveContourPolydata);
	InteractiveContourMapper->ScalarVisibilityOff();
	InteractiveContourActor = vtkSmartPointer<vtkActor>::New();
	InteractiveContourActor->SetMapper(InteractiveContourMapper);
	InteractiveContourActor->GetProperty()->SetColor(217/255.0, 95.0/255.0, 14.0/255.0);
	InteractiveContourActor->GetProperty()->SetLineWidth(SystemSettings->size_contour);
	//=============================================================================================
}

void radImageView::DrawContours()
{
	//draw cone contours
	ContourPoints = vtkSmartPointer<vtkPoints>::New();
	ContourCells = vtkSmartPointer<vtkCellArray>::New();
	ContourPolydata = vtkSmartPointer<vtkPolyData>::New();
	ContourPolydata->SetPoints(ContourPoints);
	ContourPolydata->SetLines(ContourCells);
	ContourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	ContourMapper->SetInputData(ContourPolydata);
	ContourMapper->ScalarVisibilityOff();
	ContourActor = vtkSmartPointer<vtkActor>::New();
	ContourActor->SetMapper(ContourMapper);
	ContourActor->GetProperty()->SetColor(ResultColor[0], ResultColor[1], ResultColor[2]);
	ContourActor->GetProperty()->SetLineWidth(SystemSettings->size_contour);

	//draw cone centerpoints
	CenterPoints = vtkSmartPointer<vtkPoints>::New();
	CenterPolydata = vtkSmartPointer<vtkPolyData>::New();
	CenterPolydata->SetPoints(CenterPoints);

	CenterGlyphSource = vtkSmartPointer<vtkGlyphSource2D>::New();
	CenterGlyphSource->SetGlyphTypeToCross();
	CenterGlyphSource->SetScale(4);

	CenterGlyph = vtkSmartPointer<vtkGlyph3D>::New();
	CenterGlyph->SetSourceConnection(CenterGlyphSource->GetOutputPort());
	CenterGlyph->SetInputData(CenterPolydata);
	CenterGlyph->SetRange(0, 1);
	CenterGlyph->SetColorModeToColorByScalar();
	CenterGlyph->SetScaleModeToDataScalingOff();

	CenterMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	CenterMapper->SetInputConnection(CenterGlyph->GetOutputPort());
	CenterMapper->ScalarVisibilityOff();
	CenterMapper->SetScalarRange(0, 1);
	CenterMapper->SetColorModeToMapScalars();
	CenterMapper->SetScalarModeToUsePointData();

	CenterActor = vtkSmartPointer<vtkActor>::New();
	CenterActor->SetMapper(CenterMapper);
	CenterActor->GetProperty()->SetColor(ResultColor[0], ResultColor[1], ResultColor[2]);

	//draw cone regions
	RegionPoints = vtkSmartPointer<vtkPoints>::New();
	RegionCells = vtkSmartPointer<vtkCellArray>::New();
	RegionPolyData = vtkSmartPointer<vtkPolyData>::New();
	RegionPolyData->SetPoints(RegionPoints);
	RegionPolyData->SetPolys(RegionCells);
	RegionTriangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	RegionTriangleFilter->SetInputData(RegionPolyData);
	RegionMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	RegionMapper->SetInputConnection(RegionTriangleFilter->GetOutputPort());
	RegionMapper->ScalarVisibilityOff();
	RegionActor = vtkSmartPointer<vtkActor>::New();
	RegionActor->SetMapper(RegionMapper);
	RegionActor->GetProperty()->SetColor(RegionColor[0], RegionColor[1], RegionColor[2]);
	RegionActor->SetVisibility(false);

}

void radImageView::DrawEditedContours()
{
	//add segmentation contour widget
	EditedContourPoints = vtkSmartPointer<vtkPoints>::New();
	EditedContourCells = vtkSmartPointer<vtkCellArray>::New();
	EditedContourPoly = vtkSmartPointer<vtkPolyData>::New();
	EditedContourPoly->SetPoints(EditedContourPoints);
	EditedContourPoly->SetLines(EditedContourCells);

	EditedContourWidget = vtkSmartPointer<vtkContourWidget>::New();
	EditedContourWidget->ContinuousDrawOff();
	EditedContourRepresentation = vtkSmartPointer<vtkOrientedGlyphContourRepresentation>::New();
    EditedContourRepresentation->GetLinesProperty()->SetColor(1., 1., 0);	// Edited contour lines
	EditedContourRepresentation->GetActiveProperty()->SetColor(1, 0.5, 0.5);	// Outlines around active contour edges
	EditedContourRepresentation->GetProperty()->SetColor(1., 0.5, 0.5);		// Edited contour edges
	EditedContourRepresentation->SetLineInterpolator(vtkLinearContourLineInterpolator::New());

    EditedContourWidget->SetRepresentation(EditedContourRepresentation);
}

radImageView::radImageView(MarkerSystemSettings * settings)
{
	SystemSettings = settings;
	
	//segmentation visualization
    DrawInputImage();

	//manual input
	DrawInteractiveContours();
	
	//draw contours
	DrawContours();

	//draw edited contour
	DrawEditedContours();

	//=============================================================================================
	ImageRender = vtkSmartPointer<vtkRenderer>::New();
	//add dose contour actor
	ImageRender->SetBackground( 0.0f, 0.0f, 0.0f );
	ImageRender->AddActor(ImageActor);
	ImageRender->AddActor(CenterActor);
	ImageRender->AddActor(InteractiveContourActor);
	ImageRender->AddActor(ContourActor);
	ImageRender->AddActor(RegionActor);

	ImageStyle = vtkSmartPointer<radMouseInteractorStylePP>::New();
	RenderWin = vtkSmartPointer<vtkRenderWindow>::New();
    RenderWin->AddRenderer(ImageRender);
	WinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    WinInteractor->SetInteractorStyle(ImageStyle);
    WinInteractor->SetRenderWindow(RenderWin);

	EditedContourWidget->SetInteractor(WinInteractor);
	EditedContourWidget->On();
	EditedContourWidget->Initialize();
	EditedContourWidget->Off();
	vtkSmartPointer<callbackContourWidget> callback = vtkSmartPointer<callbackContourWidget>::New();
	EditedContourWidget->AddObserver(vtkCommand::InteractionEvent,callback);
	
	WinInteractor->Initialize();

	ApplySystemSettings();
}

radImageView::~radImageView()
{
}

void radImageView::ApplySystemSettings()
{
	SetConeContourColor(SystemSettings->contour_color.toDoubleArr());
	SetConeCenterColor(SystemSettings->center_color.toDoubleArr());
	SetConeRegionColor(SystemSettings->region_color.toDoubleArr());
	SetConeContourVisibility(SystemSettings->visibility_contour);
	SetConeCenterVisibility(SystemSettings->visibility_center);
	SetConeRegionVisibility(SystemSettings->visibility_region);
	SetConeRegionOpacity(SystemSettings->region_opacity);
	SetCenterGlyphScale(SystemSettings->size_center);
}

void radImageView::SetSplitImage(RGBImageType::Pointer img)
{
	if (!img)
		return;

	ImageData->Initialize();
	ConvertITKImageToVTKImage<RGBImageType, unsigned char>(img, 3, ImageData, VTK_UNSIGNED_CHAR);
	ImageData->Modified();
}

void radImageView::SetColorInfo(ColorInfo ci)
{
	ImageActor->GetProperty()->SetColorLevel(ci.color_level);
	ImageActor->GetProperty()->SetColorWindow(ci.color_window);
}
ColorInfo radImageView::GetColorInfo()
{
	ColorInfo ci;
	ci.color_level = ImageActor->GetProperty()->GetColorLevel();
	ci.color_window = ImageActor->GetProperty()->GetColorWindow();
	return ci;
}

void radImageView::InitializeView()
{
	InteractiveContourPoints->Initialize();
	InteractiveContourCells->Initialize();
	InteractiveContourPoints->Modified();
	InteractiveContourCells->Modified();
	InteractiveContourPolydata->Modified();

	ContourPoints->Initialize();
	ContourCells->Initialize();
	ContourPoints->Modified();
	ContourCells->Modified();
	ContourPolydata->Modified();

	RegionPoints->Initialize();
	RegionCells->Initialize();
	RegionPoints->Modified();
	RegionCells->Modified();
	RegionPolyData->Modified();

	EditedContourPoints->Initialize();
	EditedContourCells->Initialize();
	EditedContourPoints->Modified();
	EditedContourCells->Modified();
	EditedContourPoly->Modified();

	CenterPoints->Initialize();
	CenterPoints->SetNumberOfPoints(0);
	CenterPolydata->Modified();
	// CenterGlyph->Update();
}

void radImageView::SetInterpolation(bool flag)
{
	interpolationFlag = flag;
	if (interpolationFlag)
		ImageActor->InterpolateOn();
	else
		ImageActor->InterpolateOff();
}

void radImageView::ChangeCameraOrientation(vtkSmartPointer<vtkRenderer> img_ren)
{
	img_ren->ResetCamera();
	double* fp = img_ren->GetActiveCamera()->GetFocalPoint();
	double* p = img_ren->GetActiveCamera()->GetPosition();
	double dist = sqrt( (p[0]-fp[0])*(p[0]-fp[0]) + (p[1]-fp[1])*(p[1]-fp[1]) 
		+ (p[2]-fp[2])*(p[2]-fp[2]) );
	img_ren->GetActiveCamera()->SetPosition(fp[0], fp[1], fp[2]-dist);
	img_ren->GetActiveCamera()->SetViewUp(0.0, -1.0, 0.0);
}

void radImageView::ResetView(bool camera_flag)
{
	if (camera_flag)
		ChangeCameraOrientation(ImageRender);
    RenderWin->Render();
}

void radImageView::SetInteractiveContours(DoublePointArray & pts, bool ending_flag)
{
	InteractiveContourPoints->Initialize();
	InteractiveContourPoints->SetNumberOfPoints(pts.size());
	InteractiveContourCells->Initialize();
	InteractiveContourCells->InsertNextCell(pts.size());

	for (unsigned int i=0; i<pts.size(); i++)
	{
		InteractiveContourPoints->SetPoint(i, pts[i][0], pts[i][1], SmallDisplacement);
		InteractiveContourCells->InsertCellPoint(i);
	}
	InteractiveContourPoints->Modified();
	InteractiveContourCells->Modified();
	InteractiveContourPolydata->Modified();

	if (ending_flag)
	{
		InteractiveContourPoints->Initialize();
		InteractiveContourCells->Initialize();
		InteractiveContourPoints->Modified();
		InteractiveContourCells->Modified();
		InteractiveContourPolydata->Modified();
	}
}

void radImageView::SetContourMarkers(vector< ContourMarker > &cone_marker_shapes)
{
	RegionPoints->Initialize();
	RegionCells->Initialize();
	ContourPoints->Initialize();
	ContourCells->Initialize();
	CenterPoints->Initialize();
	CenterPoints->SetNumberOfPoints(cone_marker_shapes.size());

	unsigned jj = 0;
	for (unsigned int i = 0; i<cone_marker_shapes.size(); i++)
	{
		if (cone_marker_shapes[i].is_visible)
		{
			// Add contour centers
			DoublePointType &cpt = cone_marker_shapes[i].marker_center;
			CenterPoints->SetPoint(jj, cpt[0], cpt[1], -0.05);
			jj++;

			//add contour points
			ContourCells->InsertNextCell(cone_marker_shapes[i].marker_contours.size() + 1);
			RegionCells->InsertNextCell(cone_marker_shapes[i].marker_contours.size());
			unsigned int current_pt_num = ContourPoints->GetNumberOfPoints();
			for (unsigned int j = 0; j<cone_marker_shapes[i].marker_contours.size(); j++)
			{
				ContourPoints->InsertNextPoint(cone_marker_shapes[i].marker_contours[j][0],
					cone_marker_shapes[i].marker_contours[j][1], SmallDisplacement);
				RegionPoints->InsertNextPoint(cone_marker_shapes[i].marker_contours[j][0],
					cone_marker_shapes[i].marker_contours[j][1], SmallDisplacement);
				ContourCells->InsertCellPoint(j + current_pt_num);
				RegionCells->InsertCellPoint(j + current_pt_num);
			}
			ContourCells->InsertCellPoint(current_pt_num);
		}
	}

	CenterPoints->SetNumberOfPoints(jj);

	RegionPoints->Modified();
	RegionCells->Modified();
	RegionPolyData->Modified();
	ContourPoints->Modified();
	ContourCells->Modified();
	CenterPoints->Modified();
	CenterPolydata->Modified();
}

void radImageView::SetConeContourColor(double rgb[3])
{
	ContourActor->GetProperty()->SetColor(rgb);
}

void radImageView::SetConeCenterColor(double rgb[3])
{
	CenterActor->GetProperty()->SetColor(rgb);
}

void radImageView::SetConeContourVisibility(bool flag)
{
	ContourActor->SetVisibility(flag);
}

void radImageView::SetConeCenterVisibility(bool flag)
{
	CenterActor->SetVisibility(flag);
}

void radImageView::SetCenterGlyphScale(double scale)
{
	CenterGlyphSource->SetScale(scale);
}

void radImageView::SetConeRegionVisibility(bool flag)
{
	RegionActor->SetVisibility(flag);
}

void radImageView::SetConeRegionOpacity(double op)
{
	RegionActor->GetProperty()->SetOpacity(op);
}

void radImageView::SetConeRegionColor(double rgb[3])
{
	RegionActor->GetProperty()->SetColor(rgb);
}

void radImageView::SetContourWidth(int width)
{
	ContourActor->GetProperty()->SetLineWidth(width);
	InteractiveContourActor->GetProperty()->SetLineWidth(width);
}

void radImageView::EnableEditedContour(DoublePointArray &contour_pts)
{
	EditedContourPoints->Initialize();
	EditedContourCells->Initialize();
	
	EditedContourPoints->SetNumberOfPoints(contour_pts.size());
	EditedContourCells->InsertNextCell(contour_pts.size()+1);
	for (unsigned int i=0; i<contour_pts.size(); i++)
	{
		EditedContourPoints->SetPoint(i, contour_pts[i][0], contour_pts[i][1], SmallDisplacement);
		EditedContourCells->InsertCellPoint(i);
	}
	EditedContourCells->InsertCellPoint(0);
	EditedContourPoly->Modified();

	EditedContourWidget->On();
	EditedContourWidget->Initialize(EditedContourPoly);
}

void radImageView::DisableEditedContour()
{
	EditedContourPoints->Initialize();
	EditedContourCells->Initialize();
	EditedContourPoly->Modified();

	EditedContourWidget->Off();
}
