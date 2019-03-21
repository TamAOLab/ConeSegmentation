#include "radSegmentationPanel.h"

radSegmentationPanel::radSegmentationPanel(QWidget *parent)
	: QDialog(parent)
{
	setWindowTitle(tr("Split Image Cone Segmentation"));
	setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
	segmentationParameters.loadDefaults();
	CreateInputGroup();
}

radSegmentationPanel::~radSegmentationPanel()
{

}

void radSegmentationPanel::SetParameters(ConeSegmentationParameters &params)
{
	segmentationParameters = params;
	UpdateParameters();
}

void radSegmentationPanel::UpdateParameters()
{
	HessianResponseInput->setValue(segmentationParameters.HessianThreshold);
	GacInput->setValue(segmentationParameters.GACIterationNumber);
}

void radSegmentationPanel::CreateInputGroup()
{
	SegmentationSetupGroup = new QGroupBox(tr("Cone Segmentation Setup"));

	HessianResponseLabel = new QLabel(tr("Hessian Response:"));
	HessianResponseInput = new QDoubleSpinBox;
	HessianResponseInput->setRange(10.0, 500.);
	HessianResponseInput->setSingleStep(10.0);
	HessianResponseInput->setAlignment(Qt::AlignCenter);
	connect(HessianResponseInput, SIGNAL(valueChanged(double)), this, SLOT(ChangeHessianResponse(double)));

	GacLabel = new QLabel(tr("GAC Iteration #:"));
	GacInput = new QSpinBox;
	GacInput->setRange(1, 150);
	GacInput->setSingleStep(1);
	GacInput->setAlignment(Qt::AlignCenter);
	connect(GacInput, SIGNAL(valueChanged(int)), this, SLOT(ChangeGac(int)));

	SegmentationSetupLayout = new QGridLayout;
	SegmentationSetupLayout->setHorizontalSpacing(30);
	SegmentationSetupLayout->setColumnMinimumWidth(1, 150);
	SegmentationSetupLayout->setColumnMinimumWidth(3, 150);

	SegmentationSetupLayout->addWidget(HessianResponseLabel, 0, 0);
	SegmentationSetupLayout->addWidget(HessianResponseInput, 0, 1);
	SegmentationSetupLayout->addWidget(GacLabel, 0, 2);
	SegmentationSetupLayout->addWidget(GacInput, 0, 3);

	QFrame *line = new QFrame(this);
	line->setFrameShape(QFrame::HLine);
	line->setFrameShadow(QFrame::Sunken);
	SegmentationSetupLayout->addWidget(line, 1, 0, 1, 4);

	RestoreDefaultsButton = new QPushButton("  Restore Defaults  ");
	SegmentationSetupLayout->addWidget(RestoreDefaultsButton, 2, 2, 1, 2, Qt::AlignRight | Qt::AlignBottom);
	connect(RestoreDefaultsButton, SIGNAL(clicked()), this, SLOT(ClickedRestoreDefaults()));

	SegmentationSetupGroup->setLayout(SegmentationSetupLayout);

	SegmentationLaunchGroup = new QWidget();
	SegmentationLaunchLayout = new QGridLayout;
	SegmentationLaunchGroup->setLayout(SegmentationLaunchLayout);

	LaunchCurrentButton = new QPushButton(tr("Segment Current"));
	SegmentationLaunchLayout->addWidget(LaunchCurrentButton, 0, 0);
	connect(LaunchCurrentButton, SIGNAL(clicked()), this, SLOT(ClickedSegmentCurrent()));
	LaunchAllButton = new QPushButton(tr("Segment All"));
	SegmentationLaunchLayout->addWidget(LaunchAllButton, 0, 1);
	connect(LaunchAllButton, SIGNAL(clicked()), this, SLOT(ClickedSegmentAll()));

	ViewLayout = new QVBoxLayout;
	ViewLayout->addWidget(SegmentationSetupGroup);
	ViewLayout->addWidget(SegmentationLaunchGroup);
	setLayout(ViewLayout);

	RestoreDefaultsButton->setAutoDefault(false);
	LaunchCurrentButton->setAutoDefault(false);
	LaunchAllButton->setAutoDefault(true);

	UpdateParameters();
}

void radSegmentationPanel::ChangeHessianResponse(double value)
{
	segmentationParameters.HessianThreshold = value;
}

void radSegmentationPanel::ChangeGac(int value)
{
	segmentationParameters.GACIterationNumber = value;
}

void radSegmentationPanel::ClickedSegmentCurrent()
{
	close();
	emit launchSegmentCurrent();
}

void radSegmentationPanel::ClickedSegmentAll()
{
	close();
	emit launchSegmentAll();
}

void radSegmentationPanel::ClickedRestoreDefaults()
{
	segmentationParameters.loadDefaults();
	UpdateParameters();
}

