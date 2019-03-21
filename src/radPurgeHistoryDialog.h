#pragma once
#include <QtGlobal>
#include <QtGui>
#include <QtWidgets>
#include <QVBoxLayout>
#include <qdialog.h>

#include "radBackup.h"

class radPurgeHistoryDialog :
	public QDialog
{
Q_OBJECT
public:
	radPurgeHistoryDialog(QWidget *parent);
	virtual ~radPurgeHistoryDialog();

	void showHistory(std::string historyDir);
private:
	QVBoxLayout *layout;
	QTableWidget *table;

	QLabel *lblInfo;

	QWidget *buttonGroup;
	QGridLayout *buttonLayout;
	QPushButton *btnSelectAll;
	QPushButton *btnDeselectAll;
	QPushButton *btnDelete;
	QPushButton *btnClose;

	QList<QString> namelist;
	std::string historyDir;

	void loadHistory();

private slots:
	void handleDelete();
	void handleSelectAll();
	void handleDeselectAll();
	void handleSelectionChange(const QItemSelection &selected, const QItemSelection &deselected);
};

