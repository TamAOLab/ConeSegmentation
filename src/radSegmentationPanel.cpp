#include "radSegmentationPanel.h"

radSegmentationPanel::radSegmentationPanel(QWidget *parent)
	: QDialog(parent)
{
	setWindowTitle(tr("Split Image Cone Segmentation"));
	setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
	segmentationParameters.loadDefaults();
	CreateInputGroup();

	normal = QFont(this->font());
	bold = QFont(normal);
	bold.setBold(true);
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
	SegmentationLaunchLayout->setColumnStretch(0, 10);
	SegmentationLaunchLayout->setColumnStretch(1, 10);
	SegmentationLaunchLayout->setColumnStretch(2, 10);

	SegmentationLaunchGroup->setLayout(SegmentationLaunchLayout);

	LaunchCheckedButton = new QPushButton(tr("Segment Checked"));
	SegmentationLaunchLayout->addWidget(LaunchCheckedButton, 0, 1);
	connect(LaunchCheckedButton, SIGNAL(clicked()), this, SLOT(ClickedSegmentChecked()));
	CancelButton = new QPushButton(tr("Cancel"));
	SegmentationLaunchLayout->addWidget(CancelButton, 0, 2);
	connect(CancelButton, SIGNAL(clicked()), this, SLOT(close()));

	imageTable = new QTableWidget(0, 2);

	QStringList headers;
	headers << "\xE2\x88\x9A" << "Split File Name";

	imageTable->setHorizontalHeaderLabels(headers);
	imageTable->setColumnWidth(0, 12);
	imageTable->verticalHeader()->setVisible(false);
	imageTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
	imageTable->setSelectionBehavior(QAbstractItemView::SelectRows);
	imageTable->setSelectionMode(QAbstractItemView::SingleSelection);
	imageTable->setShowGrid(false);

	QHeaderView *hdr = imageTable->horizontalHeader();
	hdr->setSectionResizeMode(0, QHeaderView::ResizeToContents);
	hdr->setSectionResizeMode(1, QHeaderView::Stretch);
	hdr->setSectionsClickable(true);
	connect(hdr, SIGNAL(sectionClicked(int)), SLOT(onHeaderClicked(int)));

	ViewLayout = new QGridLayout();

	ViewLayout->setColumnStretch(0, 0);
	ViewLayout->setColumnStretch(1, 10);
	ViewLayout->setRowStretch(0, 10);
	ViewLayout->setRowStretch(1, 0);
	ViewLayout->setRowStretch(2, 0);
	ViewLayout->addWidget(imageTable, 0, 0, 1, 2);
	ViewLayout->addWidget(SegmentationSetupGroup, 1, 0);
	ViewLayout->addWidget(SegmentationLaunchGroup, 2, 0);
	setLayout(ViewLayout);

	RestoreDefaultsButton->setAutoDefault(false);
	CancelButton->setAutoDefault(false);
	LaunchCheckedButton->setAutoDefault(true);

	UpdateParameters();
}

void radSegmentationPanel::SetItemList(QStringList &items)
{
	imageTable->setRowCount(items.size());
	for (int row = 0; row < items.size(); row++) {
		QCheckBox *cb = new QCheckBox();
		cb->setContentsMargins(8, 2, 2, 0);
		imageTable->setCellWidget(row, 0, cb);
		imageTable->setItem(row, 1, new QTableWidgetItem(items[row]));
	}
	imageTable->resizeColumnsToContents();
	imageTable->resizeRowsToContents();
}
void radSegmentationPanel::SetCheckedRows(QList<int> &rows)
{
	checkedRows = rows;
	for (int row = 0; row < imageTable->rowCount(); row++) {
		((QCheckBox *)(imageTable->cellWidget(row, 0)))->setChecked(rows.indexOf(row) >= 0);
	}
}
void radSegmentationPanel::SetHighlightedRow(int row)
{
	for (int i = 0; i < imageTable->rowCount(); i++) {
		imageTable->item(i, 1)->setFont(normal);
	}
	if (row < 0 || row >= imageTable->rowCount()) return;
	imageTable->selectRow(row);
	imageTable->item(row, 1)->setFont(bold);
	imageTable->scrollToItem(imageTable->item(row, 1));
}
void radSegmentationPanel::onHeaderClicked(int hdr)
{
	if (hdr != 0 || imageTable->rowCount() == 0) return;
	int nchecked = 0;
	for (int row = 0; row < imageTable->rowCount(); row++) {
		if (((QCheckBox *)(imageTable->cellWidget(row, 0)))->isChecked()) {
			++nchecked;
			((QCheckBox *)(imageTable->cellWidget(row, 0)))->setChecked(false);
		}
	}
	if (nchecked > 0) return;
	for (int row = 0; row < imageTable->rowCount(); row++) {
		((QCheckBox *)(imageTable->cellWidget(row, 0)))->setChecked(true);
	}
}

void radSegmentationPanel::closeEvent(QCloseEvent *event)
{
	dlgeom = geometry();
	dlgeomset = true;
	QDialog::closeEvent(event);
}
void radSegmentationPanel::showEvent(QShowEvent *event)
{
	if (dlgeomset) setGeometry(dlgeom);
	QDialog::showEvent(event);
}

void radSegmentationPanel::ChangeHessianResponse(double value)
{
	segmentationParameters.HessianThreshold = value;
}

void radSegmentationPanel::ChangeGac(int value)
{
	segmentationParameters.GACIterationNumber = value;
}

void radSegmentationPanel::ClickedSegmentChecked()
{
	QList<int> checked;
	for (int row = 0; row < imageTable->rowCount(); row++) {
		if (((QCheckBox *)(imageTable->cellWidget(row, 0)))->isChecked()) {
			checked << row;
		}
	}
	if (checked.size() == 0) {
		return;
	}
	close();
	emit launchSegmentChecked(checked);
}

void radSegmentationPanel::ClickedRestoreDefaults()
{
	segmentationParameters.loadDefaults();
	UpdateParameters();
	SetCheckedRows(checkedRows);
}

