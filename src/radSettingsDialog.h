/*
 *  radSettingsDialog.h
 *  
 *
 *  Created by Jianfei Liu on 9/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef radSettingsDialog_H
#define radSettingsDialog_H

#include "radDefinition.h"

#include <QtGui>
#include <QLabel>
#include <QGroupBox>
#include <QCheckBox>
#include <QPushButton>
#include <QAction>
#include <QMenu>
#include <QSpinBox>
#include <QColorDialog>
#include <QSlider>

#include "radColorButton.h"

class radSettingsDialog : public QDialog
{
  Q_OBJECT
public:
	
	radSettingsDialog(MarkerSystemSettings *, QWidget *parent = 0);
	~radSettingsDialog();
	void updateSettings();

private slots:
	void SetContourVisibility(bool);
	void SetCenterVisibility(bool);
	void SetRegionVisibility(bool);

	void SetContourWidth(int);
	void SetCenterSize(double);
	void SetContourColor(QColor color);
	void SetCenterColor(QColor color);
	void SetRegionColor(QColor color);
	void SetRegionOpacity(int op);
	void RestoreDefaults();

private:

	MarkerSystemSettings *SystemSettings;

	QCheckBox *contour_visibility_box;
	QCheckBox *center_visibility_box;
	QCheckBox *region_visibility_box;
	radColorButton *contour_color;
	radColorButton *center_color;
	radColorButton *region_color;
	QSpinBox *contour_width_input;
	QDoubleSpinBox *center_size_input;
	QSlider *contour_opacity;
	QLabel *opacity_text;

	void CreateWidgetGroup();
	void UpdateOpacityText();
};

#endif // radProgressDialog

