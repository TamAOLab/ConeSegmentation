/*
 *  radSettingsDialog.cpp
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "radSettingsDialog.h"
#include "radmainwindow.h"


radSettingsDialog::radSettingsDialog(MarkerSystemSettings *settings, QWidget *parent)
    : QDialog(parent)
{
	SystemSettings = settings;

	//create input group
	CreateWidgetGroup();
	//this->setWindowFlags( ( (this->windowFlags() | Qt::CustomizeWindowHint) & ~Qt::WindowCloseButtonHint) );
	setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
	setWindowTitle(tr("Display Settings"));
}

radSettingsDialog::~radSettingsDialog()
{ 
}

void radSettingsDialog::updateSettings()
{
	contour_color->setColor(SystemSettings->contour_color);
	contour_visibility_box->setChecked(SystemSettings->visibility_contour);
	center_color->setColor(SystemSettings->center_color);
	center_visibility_box->setChecked(SystemSettings->visibility_center);
	contour_width_input->setValue(SystemSettings->size_contour);
	center_size_input->setValue(SystemSettings->size_center);
	region_color->setColor(SystemSettings->region_color);
	region_visibility_box->setChecked(SystemSettings->visibility_region);
	contour_opacity->setValue((int)(SystemSettings->region_opacity * 100.));
	UpdateOpacityText();
}

void radSettingsDialog::UpdateOpacityText()
{
	int v = contour_opacity->value();
	opacity_text->setText(std::to_string(v).c_str());
}

void radSettingsDialog::CreateWidgetGroup()
{
	QGroupBox *contourBox = new QGroupBox(tr("Contours"));
	QGridLayout *contourLayout = new QGridLayout();
	contourBox->setLayout(contourLayout);

	contour_visibility_box = new QCheckBox("Visible", this);
	contour_visibility_box->setCheckable(true);
	contour_visibility_box->setChecked(SystemSettings->visibility_contour);
	connect(contour_visibility_box, SIGNAL(toggled(bool)), this, SLOT(SetContourVisibility(bool)));

	contour_color = new radColorButton(" Color", this, SystemSettings->contour_color);
	contour_color->setAutoDefault(false);
	connect(contour_color, SIGNAL(colorChanged(QColor)), this, SLOT(SetContourColor(QColor)));

	QLabel *contour_width_label = new QLabel(tr("Width:"));
	contour_width_input = new QSpinBox;
	contour_width_input->setRange(1, 50);
	contour_width_input->setSingleStep(1);
	contour_width_input->setValue(SystemSettings->size_contour);
	contour_width_input->setAlignment(Qt::AlignCenter);
	contour_width_input->setMinimumWidth(120);
	connect(contour_width_input, SIGNAL(valueChanged(int)), this, SLOT(SetContourWidth(int)));

	contourLayout->addWidget(contour_visibility_box, 0, 0, 1, 2);
	contourLayout->addWidget(contour_color, 1, 0, 1, 2, Qt::AlignLeft);
	contourLayout->addWidget(contour_width_label, 2, 0);
	contourLayout->addWidget(contour_width_input, 2, 1);

	QGroupBox *centerBox = new QGroupBox(tr("Center Points"));
	QGridLayout *centerLayout = new QGridLayout();
	centerBox->setLayout(centerLayout);

	center_visibility_box = new QCheckBox("Visible", this);
	center_visibility_box->setCheckable(true);
	center_visibility_box->setChecked(SystemSettings->visibility_center);
	connect(center_visibility_box, SIGNAL(toggled(bool)), this, SLOT(SetCenterVisibility(bool)));

	center_color = new radColorButton(" Color", this, SystemSettings->center_color);
	center_color->setAutoDefault(false);
	connect(center_color, SIGNAL(colorChanged(QColor)), this, SLOT(SetCenterColor(QColor)));

	QLabel *center_size_label = new QLabel(tr("Size:"));
	center_size_input = new QDoubleSpinBox;
	center_size_input->setRange(0.5, 10.);
	center_size_input->setSingleStep(0.5);
	center_size_input->setValue(SystemSettings->size_center);
	center_size_input->setAlignment(Qt::AlignCenter);
	center_size_input->setMinimumWidth(120);
	connect(center_size_input, SIGNAL(valueChanged(double)), this, SLOT(SetCenterSize(double)));

	centerLayout->addWidget(center_visibility_box, 0, 0, 1, 2);
	centerLayout->addWidget(center_color, 1, 0, 1, 2, Qt::AlignLeft);
	centerLayout->addWidget(center_size_label, 2, 0);
	centerLayout->addWidget(center_size_input, 2, 1);

	QGroupBox *regionBox = new QGroupBox(tr("Regions"));
	QGridLayout *regionLayout = new QGridLayout();
	regionBox->setLayout(regionLayout);

	region_visibility_box = new QCheckBox("Visible", this);
	region_visibility_box->setCheckable(true);
	region_visibility_box->setChecked(SystemSettings->visibility_region);
	connect(region_visibility_box, SIGNAL(toggled(bool)), this, SLOT(SetRegionVisibility(bool)));

	region_color = new radColorButton(" Color", this, SystemSettings->region_color);
	region_color->setAutoDefault(false);
	connect(region_color, SIGNAL(colorChanged(QColor)), this, SLOT(SetRegionColor(QColor)));

	QLabel * opacity_label = new QLabel(tr("Opacity:"));
	contour_opacity = new QSlider(Qt::Horizontal, this);
	contour_opacity->setRange(0, 100);
	contour_opacity->setValue(100);
	contour_opacity->setTickInterval(10);
	contour_opacity->setMinimumWidth(150);
	connect(contour_opacity, SIGNAL(valueChanged(int)), this, SLOT(SetRegionOpacity(int)));
	opacity_text = new QLabel(tr(" 100 "));
	opacity_text->setMinimumWidth(60);

	regionLayout->addWidget(region_visibility_box, 0, 0, 1, 3);
	regionLayout->addWidget(region_color, 1, 0, 1, 3, Qt::AlignLeft);
	regionLayout->addWidget(opacity_label, 2, 0);
	regionLayout->addWidget(contour_opacity, 2, 1);
	regionLayout->addWidget(opacity_text, 2, 2);

	QWidget *ButtonGroup = new QWidget();
	QGridLayout *ButtonLayout = new QGridLayout();
	ButtonGroup->setLayout(ButtonLayout);

	QPushButton *DefaultButton = new QPushButton(tr("Restore Defaults"));
	DefaultButton->setAutoDefault(false);
	ButtonLayout->addWidget(DefaultButton, 0, 0);
	connect(DefaultButton, SIGNAL(clicked()), this, SLOT(RestoreDefaults()));
	QPushButton *CloseButton = new QPushButton(tr("Close"));
	CloseButton->setDefault(true);
	CloseButton->setAutoDefault(true);
	ButtonLayout->addWidget(CloseButton, 0, 1);
	connect(CloseButton, SIGNAL(clicked()), this, SLOT(close()));

	QGridLayout *view_layout = new QGridLayout;
	view_layout->addWidget(contourBox, 0, 0);
	view_layout->addWidget(centerBox, 0, 1);
	view_layout->addWidget(regionBox, 0, 2);
	view_layout->addWidget(ButtonGroup, 1, 0, 1, 3);
	setLayout(view_layout);
}

void radSettingsDialog::SetRegionVisibility(bool flag)
{
	SystemSettings->visibility_region = flag;
	radMainWindow::GetPointer()->GetImageView()->SetConeRegionVisibility(flag);
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::SetContourVisibility(bool flag)
{
	SystemSettings->visibility_contour = flag;
	radMainWindow::GetPointer()->GetImageView()->SetConeContourVisibility(flag);
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::SetCenterVisibility(bool flag)
{
	SystemSettings->visibility_center = flag;
	radMainWindow::GetPointer()->GetImageView()->SetConeCenterVisibility(flag);
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::SetContourColor(QColor color)
{
	SystemSettings->contour_color = color;
	radMainWindow::GetPointer()->GetImageView()->SetConeContourColor(SystemSettings->contour_color.toDoubleArr());
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::SetCenterColor(QColor color)
{
	SystemSettings->center_color = color;
	radMainWindow::GetPointer()->GetImageView()->SetConeCenterColor(SystemSettings->center_color.toDoubleArr());
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::SetRegionColor(QColor color)
{
	SystemSettings->region_color = color;
	radMainWindow::GetPointer()->GetImageView()->SetConeRegionColor(SystemSettings->region_color.toDoubleArr());
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::SetRegionOpacity(int op)
{
	UpdateOpacityText();
	SystemSettings->region_opacity = op / 100.;
	radMainWindow::GetPointer()->GetImageView()->SetConeRegionOpacity(SystemSettings->region_opacity);
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::SetContourWidth(int contour_width)
{
	SystemSettings->size_contour = contour_width;
	radMainWindow::GetPointer()->GetImageView()->SetContourWidth(contour_width);
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::SetCenterSize(double value)
{
	SystemSettings->size_center = value;
	radMainWindow::GetPointer()->GetImageView()->SetCenterGlyphScale(value);
	radMainWindow::GetPointer()->GetImageView()->ResetView(false);
}

void radSettingsDialog::RestoreDefaults()
{
	MarkerSystemSettings settings;
	SetContourWidth(settings.size_contour);
	SetContourColor(settings.contour_color);
	SetContourVisibility(settings.visibility_contour);
	SetRegionColor(settings.region_color);
	SetRegionVisibility(settings.visibility_region);
	SetRegionOpacity((int)(settings.region_opacity * 100.));
	SetCenterVisibility(settings.visibility_center);
	SetCenterColor(settings.center_color);
	SetCenterSize(settings.size_center);
	updateSettings();
}
