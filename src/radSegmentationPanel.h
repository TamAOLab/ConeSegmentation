#ifndef radSegmentationPanel_h
#define radSegmentationPanel_h

#include <QDialog>
#include <QGroupBox>
#include <QLabel>
#include <QGridLayout>
#include <QGridLayout>
#include <QPushButton>
#include <QDoubleSpinBox>

#include "radDefinition.h"

class radSegmentationPanel : public QDialog
{
	Q_OBJECT
public:

	radSegmentationPanel(QWidget *parent = NULL);
	~radSegmentationPanel();

	void SetParameters(ConeSegmentationParameters &params);
	ConeSegmentationParameters & GetParameters() { return segmentationParameters; }

private slots:
	void ChangeHessianResponse(double value);
	void ChangeGac(int value);

	void ClickedSegmentCurrent();
	void ClickedSegmentAll();
	void ClickedRestoreDefaults();

signals:
	void launchSegmentCurrent();
	void launchSegmentAll();

private:
	ConeSegmentationParameters segmentationParameters;

	QGroupBox *SegmentationSetupGroup;
	QGridLayout * SegmentationSetupLayout;

	QLabel * HessianResponseLabel;
	QDoubleSpinBox * HessianResponseInput;
	QLabel * GacLabel;
	QSpinBox * GacInput;

	QWidget * SegmentationLaunchGroup;
	QGridLayout * SegmentationLaunchLayout;

	QPushButton * RestoreDefaultsButton;

	QPushButton * LaunchCurrentButton;
	QPushButton * LaunchAllButton;

	QVBoxLayout * ViewLayout;

	void UpdateParameters();
	void CreateInputGroup();
};

#endif
