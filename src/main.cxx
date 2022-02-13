
#include <QApplication>
#include <itkTextOutput.h>
#include "radmainwindow.h"

int main(int argc, char** argv)
{
	itk::OutputWindow::SetInstance(itk::TextOutput::New());
	QApplication app(argc, argv);

	// Opening files from the command line, will work for drag and drop Split/Segmentation files on the App icon
	QStringList inputNames;
	for (int i = 1; i < argc; i++) {
		inputNames.append(QString(argv[i]));
	}
	QStringList splitFileNames, segmentationFileNames;
	decodeFileList(inputNames, splitFileNames, segmentationFileNames);

	radMainWindow window;
	window.setWindowIcon(QIcon(":SegmentIcon.png"));
	window.show();

	if (!splitFileNames.isEmpty()) {
		window.openSplitImages(splitFileNames, true);
		window.loadSegmentations(segmentationFileNames);
	}
	else {
		window.checkWhatsNew();
	}

	app.exec();
  
	return EXIT_SUCCESS;
}
