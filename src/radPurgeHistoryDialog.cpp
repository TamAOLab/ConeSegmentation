#include <iostream>
#include "radPurgeHistoryDialog.h"

radPurgeHistoryDialog::radPurgeHistoryDialog(QWidget *parent)
{
	setModal(true);
	setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
	setWindowIcon(QIcon(":SegmentIcon.png"));
	setWindowTitle(tr("Purge Cone Segmentation History"));

	layout = new QVBoxLayout();

	table = new QTableWidget();

	QStringList headers;
	headers << "Last Modified" << "Name";

	table->setColumnCount(2);
	table->setRowCount(0);
	table->setHorizontalHeaderLabels(headers);
	table->verticalHeader()->setVisible(false);
	table->setEditTriggers(QAbstractItemView::NoEditTriggers);
	table->setSelectionBehavior(QAbstractItemView::SelectRows);
	table->setSelectionMode(QAbstractItemView::MultiSelection);
	table->setShowGrid(false);

	QHeaderView *verticalHeader = table->verticalHeader();
	verticalHeader->setSectionResizeMode(QHeaderView::Fixed);

	int vsize = verticalHeader->defaultSectionSize() * 6 / 10;
	verticalHeader->setDefaultSectionSize(vsize);

	table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
	table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);

	layout->addWidget(table);

	lblInfo = new QLabel(tr(""), this);
	layout->addWidget(lblInfo);

	buttonGroup = new QWidget();
	buttonLayout = new QGridLayout();
	buttonGroup->setLayout(buttonLayout);

	btnSelectAll = new QPushButton(tr("Select All"), this);
	buttonLayout->addWidget(btnSelectAll, 0, 0);
	connect(btnSelectAll, SIGNAL(clicked()), SLOT(handleSelectAll()));

	btnDeselectAll = new QPushButton(tr("Deselect All"), this);
	buttonLayout->addWidget(btnDeselectAll, 0, 1);
	connect(btnDeselectAll, SIGNAL(clicked()), SLOT(handleDeselectAll()));

	btnDelete = new QPushButton(tr("Remove Selected"), this);
	buttonLayout->addWidget(btnDelete, 0, 2);
	connect(btnDelete, SIGNAL(clicked()), SLOT(handleDelete()));

	btnClose = new QPushButton(tr("Close"), this);
	buttonLayout->addWidget(btnClose, 0, 3);
	connect(btnClose, SIGNAL(clicked()), SLOT(close()));

	connect(table->selectionModel(), SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)), SLOT(handleSelectionChange(const QItemSelection &, const QItemSelection &)));

	layout->addWidget(buttonGroup);

	setLayout(layout);
}

radPurgeHistoryDialog::~radPurgeHistoryDialog()
{
}

void radPurgeHistoryDialog::handleSelectAll()
{
	table->selectAll();
	table->setFocus();
}

void radPurgeHistoryDialog::handleDeselectAll()
{
	table->clearSelection();
	table->setFocus();
}

void radPurgeHistoryDialog::handleDelete()
{
	QItemSelectionModel *select = table->selectionModel();
	if (!select->hasSelection()) return;
	QModelIndexList selection = select->selectedRows();
	QMessageBox msgBox;
	string msg = "You are about to delete " + std::to_string(selection.size()) + " entries from the segmentation history.\n"+
		"This will free up some disk space, but any unsaved data in these entries will be lost.";
	msgBox.setText(msg.c_str());
	msgBox.setInformativeText("Do you want to continue?");
	msgBox.setIcon(QMessageBox::Warning);
	msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
	msgBox.setDefaultButton(QMessageBox::Yes);
	if (msgBox.exec() == QMessageBox::Yes) {
		radBackup back_up;
		back_up.SetBackupDir(historyDir);
		QDir hdir(historyDir.c_str());
		for (auto idx : selection) {
			QString fn = table->item(idx.row(), 1)->text();
			QFileInfo fpath(hdir, fn);
			back_up.RemoveDir(fpath.absoluteFilePath());
		}
		loadHistory();
	}
}

void radPurgeHistoryDialog::handleSelectionChange(const QItemSelection &selected, const QItemSelection &deselected)
{
	QItemSelectionModel *select = table->selectionModel();
	if (select->hasSelection()) {
		QModelIndexList selection = select->selectedRows();
		std::string txt = std::to_string(selection.size()) + " of " + std::to_string(table->rowCount()) + " entries selected.";
		lblInfo->setText(txt.c_str());
		btnDelete->setEnabled(true);
	}
	else {
		lblInfo->setText("");
		btnDelete->setEnabled(false);
	}
	lblInfo->update();
}

void radPurgeHistoryDialog::loadHistory()
{
	namelist.clear();
	QDirIterator it(historyDir.c_str(), QDir::AllDirs | QDir::System | QDir::Hidden | QDir::NoDotAndDotDot);
	while (it.hasNext())
	{
		namelist.append(it.next());
	}
	table->setSortingEnabled(false);
	table->setRowCount(namelist.size());
	for (int iRow = 0; iRow < namelist.size(); iRow++) {
		QFileInfo info(namelist[iRow]);

		QDateTime dt = info.lastModified();

		table->setItem(iRow, 0, new QTableWidgetItem(dt.toString(tr("yyyy-MM-dd hh:mm:ss  "))));
		table->setItem(iRow, 1, new QTableWidgetItem(info.baseName()));
	}
	table->sortByColumn(0, Qt::DescendingOrder);
	table->setSortingEnabled(true);

	btnDelete->setEnabled(false);
}

void radPurgeHistoryDialog::showHistory(std::string historyDir)
{
	this->historyDir = historyDir;
	loadHistory();
	exec();
}