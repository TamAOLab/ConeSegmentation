#ifndef radSegmentationPanel_h
#define radSegmentationPanel_h

#include <QDialog>
#include <QGroupBox>
#include <QLabel>
#include <QGridLayout>
#include <QPushButton>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QTableWidget>
#include <QHeaderView>

#include "radDefinition.h"

class radSegmentationPanel : public QDialog
{
	Q_OBJECT
public:

	radSegmentationPanel(QWidget *parent = NULL);
	~radSegmentationPanel();

	void SetParameters(ConeSegmentationParameters &params);
	ConeSegmentationParameters & GetParameters() { return segmentationParameters; }

	void SetItemList(QStringList &items);
	void SetCheckedRows(QList<int> &rows);
	void SetHighlightedRow(int row);

protected:
	virtual void closeEvent(QCloseEvent *event) override;
	virtual void showEvent(QShowEvent *event) override;

private slots:
	void ChangeHessianResponse(double value);
	void ChangeGac(int value);

	void ClickedSegmentChecked();
	void ClickedRestoreDefaults();
	void onHeaderClicked(int);

signals:
	void launchSegmentChecked(QList<int> checked);

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

	QPushButton * LaunchCheckedButton;
	QPushButton * CancelButton;

	QGridLayout *ViewLayout;
	QTableWidget *imageTable;

	QRect dlgeom;
	bool dlgeomset = false;
	QFont normal, bold;
	QList<int> checkedRows;

	void UpdateParameters();
	void CreateInputGroup();
};

#endif
