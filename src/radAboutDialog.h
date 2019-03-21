#pragma once
#include <QtGlobal>
#include <QtGui>
#include <QtWidgets>
#include <QVBoxLayout>
#include <qdialog.h>
class radAboutDialog :
	public QDialog
{
public:
	radAboutDialog(QWidget *parent);
	virtual ~radAboutDialog();
private:
	QVBoxLayout *layout;
	QLabel *lblAbout;
	QPushButton *btnOk;
};

