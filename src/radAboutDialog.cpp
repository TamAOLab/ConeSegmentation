#include "radAboutDialog.h"

extern std::string Segm_VERSION;

static const char *atext1 = "<table><tr>\
<td><img src=\":SegmentIcon256x256.png\">&nbsp;&nbsp;</td>\
<td><b>Cone Segmentation ";
static const char *atext2 = "</b><div>\
Tam lab<br>\
National Eye Institute<br>\
National Institutes of Health</div><div><br>\
If any portion of this software is used, please<br>\
cite the following paper in your publication:\
</div>\
</td></tr><tr><td colspan=2>\
<b>Jianfei Liu, HaeWon Jung, Alfredo Dubra, and Johnny Tam</b>, <br>\
\"Cone Photoreceptor Cell Segmentation and Diameter Measurement on <br>\
Adaptive Optics Images Using Circularly Constrained Active Contour Model,\"<br> \
<i>Investigative Ophthalmology & Visual Science</i>, \
59(11): 4639 - 4652, 2018</td></tr></table>";

radAboutDialog::radAboutDialog(QWidget *parent)
	: QDialog(parent)
{
	setModal(true);
	setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
	setWindowIcon(QIcon(":about.png"));
	setWindowTitle(tr("About Cone Segmentation"));

	layout = new QVBoxLayout();
	lblAbout = new QLabel(tr(atext1)+tr(Segm_VERSION.c_str())+tr(atext2));
	lblAbout->setTextFormat(Qt::RichText);
	btnOk = new QPushButton("OK");
	layout->addWidget(lblAbout);
	layout->addWidget(btnOk);

	connect(btnOk, SIGNAL(clicked()), this, SLOT(close()));
	
	setLayout(layout);
}


radAboutDialog::~radAboutDialog()
{
	delete layout;
	delete lblAbout;
	delete btnOk;
}
