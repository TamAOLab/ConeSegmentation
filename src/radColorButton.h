#ifndef radColorButton_h
#define radColorButton_h

#include <QPushButton>
#include <QImage>
#include <QPixmap>
#include <QColorDialog>

class radColorButton : public QPushButton
{
Q_OBJECT
public:
	radColorButton(const char *text, QWidget *parent, QColor color);
	~radColorButton();

	void setColor(QColor color);
	QColor getColor() { return color; }
	void setIconSize(int iconWidth, int iconHeight);
signals:
	void colorChanged(QColor color);
private:
	QWidget *wparent;
	int iconWidth, iconHeight;
	QColor color;

	void updateIcon();
private slots:
	void showColorPicker();
};

#endif
