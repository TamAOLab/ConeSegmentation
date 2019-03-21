
#include "radColorButton.h"

radColorButton::radColorButton(const char *text, QWidget *parent, QColor color) :
	QPushButton(tr(text), parent)
{
	this->wparent = parent;
	this->iconWidth = this->iconHeight = 32;
	this->color = color;
	updateIcon();
	connect(this, SIGNAL(clicked()), this, SLOT(showColorPicker()));
}

radColorButton::~radColorButton()
{

}

void radColorButton::setColor(QColor color)
{
	this->color = color;
	updateIcon();
}

void radColorButton::setIconSize(int iconWidth, int iconHeight)
{
	this->iconWidth = iconWidth;
	this->iconHeight = iconHeight;
	updateIcon();
}

void radColorButton::updateIcon()
{
	unsigned char r = (unsigned char)color.red();
	unsigned char g = (unsigned char)color.green();
	unsigned char b = (unsigned char)color.blue();
	QImage img(iconWidth, iconHeight, QImage::Format_RGB888);
	img.fill(Qt::black);
	for (int row = 2; row < img.height() - 2; row++) {
		unsigned char *p = img.scanLine(row);
		p += 6;
		for (int col = 2; col < img.width() - 2; col++) {
			*p++ = r;
			*p++ = g;
			*p++ = b;
		}
	}
	QPixmap pixm(img.width(), img.height());
	pixm.convertFromImage(img);
	setIcon(QIcon(pixm));
}

void radColorButton::showColorPicker()
{
	QColor color = QColorDialog::getColor(getColor(), wparent);
	if (color.isValid()) {
		setColor(color);
		
		emit colorChanged(color);
	}
}