
#pragma warning(disable : 4996)
#define _CRT_SECURE_NO_WARNINGS

#include <QtGlobal>
#if QT_VERSION >= 0x050000
#include "QuickView.h"
#include <QToolButton>
#include <QToolBar>
#include <QMenuBar>
#include <QFileDialog>
#include <QColorDialog>
#include <QButtonGroup>
#include <QComboBox>
#endif

#include "radmainwindow.h"
#include "radContourOperations.h"

#include "version.h"
std::string Segm_VERSION = __version__;

static void absFileList(QStringList &inputNames, QFileInfoList &outputFiles)
{
	QStringList filters;
	filters << "*.tif" << "*.json";
	for (int i = 0; i < inputNames.size(); i++) {
		QFileInfo file = QFileInfo(inputNames[i]);
		if (file.isDir()) {
			QDir dir = QDir(file.absoluteFilePath());
			QStringList entries = dir.entryList(filters, QDir::Files);
			for (int j = 0; j < entries.size(); j++) {
				outputFiles.append(QFileInfo(dir, entries[j]));
			}
			continue;
		}
		if (file.exists())
			outputFiles.append(file);
	}
}

// inputNames may contain *.tif/*.json files or directories
void decodeFileList(QStringList &inputNames, QStringList &splitFileNames, QStringList &detectionFileNames)
{
	QFileInfoList inputFiles;
	absFileList(inputNames, inputFiles);
	for (int i = 0; i < inputFiles.size(); i++) {
		QString fn = inputFiles[i].absoluteFilePath();
		if (fn.endsWith(".tif", Qt::CaseInsensitive))
			splitFileNames.append(fn);
		else if (fn.endsWith(".json", Qt::CaseInsensitive))
			detectionFileNames.append(fn);
	}
}

static QString GetListName(std::string &imgpath)
{
	QFileInfo qfi(imgpath.c_str());
	QString basename = qfi.completeBaseName();
	QFileInfo qcsv(qfi.dir(), basename + QString(".json"));
	if (qcsv.exists()) {
		return QString("\xE2\x88\x9A") + basename;
	}
	return QString(" ") + basename;
}

radMainWindow* radMainWindow::TheApp = NULL;

// Constructor
radMainWindow::radMainWindow() 
{
	SystemSettings = new MarkerSystemSettings();

	// Handle directories
	SegmHomeDir = QDir(QDir::home().filePath(".ConeSegmentation"));
	SegmStateFile = QFileInfo(SegmHomeDir, QString("state.json"));
	QDir historyDir = QDir(SegmHomeDir.filePath("History"));
	// create history directory
	if (!historyDir.exists())
		historyDir.mkpath(".");
	BackupDir = historyDir.path().toStdString();

	QDir systemDir = QDir(SegmHomeDir.filePath("System"));
	if (!systemDir.exists())
		systemDir.mkpath(".");
	SegmSystemFile = QFileInfo(systemDir, QString("settings.json"));

	if (SegmSystemFile.exists())
		ReadSystemSettings();
	else
		WriteSystemSettings();

	EditedContourId = -1;
	
	//create history directory
	QDir infor_dir(BackupDir.c_str());
	if (!infor_dir.exists())
		infor_dir.mkpath(".");

	MouseOperationType = Mouse_Normal;
	
	ImageView = new radImageView(SystemSettings);
	
	QScreen *screen = QGuiApplication::primaryScreen();
	QRect  screenGeometry = screen->geometry();
	screen_height = screenGeometry.height();
	screen_width = screenGeometry.width();

	createActions();
	createMenu();
	createView();
	createToolBar();

    setWindowTitle(tr("Cone Segmentation ver.") + QString(Segm_VERSION.c_str()));
	setMinimumSize(screen_width / 2, screen_height * 2 / 3);
	move(screen_width / 6, screen_height / 8);

	FileIO = new radFileIO;
	TheApp = this;

	SettingsDlg = new radSettingsDialog(SystemSettings, this);

	SegmentationPanel = new radSegmentationPanel(this);
	SegmentationPanel->setMinimumSize(screen_width / 5, screen_height / 2);
	connect(SegmentationPanel, SIGNAL(launchSegmentChecked(QList<int>)), this, SLOT(SegmentConesChecked(QList<int>)));

	createProgressDialog();
	createHelpWindow();

	purgeHistoryDialog = new radPurgeHistoryDialog(this);
	purgeHistoryDialog->setMinimumSize(screen_width / 2, screen_height / 3);

	setAcceptDrops(true);
	loadState();
}

radMainWindow::~radMainWindow()
{
	delete ImageView;
	ImageView = NULL;
	delete FileIO;
	FileIO = NULL;
	delete SettingsDlg;
	SettingsDlg = NULL;
	delete SegmentationPanel;
	SegmentationPanel = NULL;
	delete progressDialog;
	progressDialog = NULL;
	
	SplitMarkerInfor.clear();
	delete SystemSettings;
	SystemSettings = NULL;
}

void radMainWindow::closeEvent(QCloseEvent *event)
{
	helpWindow->close();
	BackupResults(CurrentImageIndex);
	WriteSystemSettings();
	event->accept();
}

void radMainWindow::dragEnterEvent(QDragEnterEvent * event)
{
	event->acceptProposedAction();
}

void radMainWindow::dropEvent(QDropEvent *e)
{
	QStringList inputNames;
	QStringList splitFileNames;
	QStringList segmentationFileNames;

	for (int i = 0; i < e->mimeData()->urls().size(); i++) {
		inputNames.append(e->mimeData()->urls()[i].toLocalFile());
	}
	decodeFileList(inputNames, splitFileNames, segmentationFileNames);

	openSplitImages(splitFileNames, true);
	loadSegmentations(segmentationFileNames);
}

void radMainWindow::createActions()
{
	openSplitImageAct = new QAction(tr("&Open"), this);
	openSplitImageAct->setIcon(QIcon(":open.png"));
	openSplitImageAct->setShortcuts(QKeySequence::Open);
	openSplitImageAct->setToolTip(tr("Open Split Image(s) (Ctrl+O)"));
	connect(openSplitImageAct, SIGNAL(triggered()), this, SLOT(openSplitImage()));

	loadSegmentationAct = new QAction(tr("&Open Segmentation Results"), this);
	connect(loadSegmentationAct, SIGNAL(triggered()), this, SLOT(loadSegmentation()));

	saveSegmentationAct = new QAction(tr("&Save"), this);
	saveSegmentationAct->setToolTip(tr("Save Current Segmentation Results (Ctrl+S)"));
	saveSegmentationAct->setShortcuts(QKeySequence::Save);
	saveSegmentationAct->setIcon(QIcon(":saveas.png"));
	connect(saveSegmentationAct, SIGNAL(triggered()), this, SLOT(saveSegmentation()));

	saveAllSegmentationsAct = new QAction(tr("Save All"), this);
	saveAllSegmentationsAct->setToolTip(tr("Save All Segmentation Results (Ctrl+A)"));
	saveAllSegmentationsAct->setIcon(QIcon(":saveall.png"));
	saveAllSegmentationsAct->setShortcut(QKeySequence(tr("Ctrl+A")));
	connect(saveAllSegmentationsAct, SIGNAL(triggered()), this, SLOT(saveAllSegmentations()));

	quitAct = new QAction(tr("&Close"), this);
	quitAct->setShortcuts(QKeySequence::Quit);
	quitAct->setToolTip(tr("Close the program (Alt+F4)"));
	connect(quitAct, SIGNAL(triggered()), this, SLOT(quit()));

	segmentConesAct = new QAction(tr("Se&gment"), this);
	segmentConesAct->setShortcut(QKeySequence(tr("Ctrl+G")));
	segmentConesAct->setIcon(QIcon(":segment.png"));
	segmentConesAct->setToolTip(tr("Segment Cones (Ctrl+G)"));
	connect(segmentConesAct, SIGNAL(triggered()), this, SLOT(showSegmentationPanel()));

	purgeHistoryAct = new QAction(tr("Purge History"), this);
	connect(purgeHistoryAct, SIGNAL(triggered()), this, SLOT(purgeHistoryFiles()));

	drawActionGroup = new QActionGroup(this);
	mouseAct = new QAction(tr("Default"), drawActionGroup);
	mouseAct->setIcon(QIcon(":mouse.png"));
	mouseAct->setShortcut(QKeySequence(tr("Ctrl+M")));
	mouseAct->setToolTip(tr("Default Mouse Mode (Ctrl+M)"));
	mouseAct->setCheckable(true);
	mouseAct->setChecked(true);
	connect(mouseAct, SIGNAL(triggered()), this, SLOT(SetMouseFlag()));

	drawContourAct = new QAction(tr("Draw"), drawActionGroup);
	drawContourAct->setIcon(QIcon(":draw_contour.png"));
	drawContourAct->setShortcut(QKeySequence(tr("Ctrl+C")));
	drawContourAct->setToolTip(tr("Draw Cone Contours (Ctrl+C)"));
	drawContourAct->setCheckable(true);
	drawContourAct->setChecked(false);
	connect(drawContourAct, SIGNAL(triggered()), this, SLOT(SetContourDrawingFlag()));

	editContourAct = new QAction(tr("Edit"), drawActionGroup);
	editContourAct->setIcon(QIcon(":edit.png"));
	editContourAct->setShortcut(QKeySequence(tr("Ctrl+E")));
	editContourAct->setToolTip(tr("Edit Cone Contours (Ctrl+E)"));
	editContourAct->setCheckable(true);
	editContourAct->setChecked(false);
	connect(editContourAct, SIGNAL(triggered()), this, SLOT(SetContourEditFlag()));

	eraseContourAct = new QAction(tr("Erase M"), drawActionGroup);
	eraseContourAct->setIcon(QIcon(":erase.png"));
	eraseContourAct->setShortcut(QKeySequence(tr("Ctrl+D")));
	eraseContourAct->setToolTip(tr("Erase Cone Contours (Ctrl+D)"));
	eraseContourAct->setCheckable(true);
	eraseContourAct->setChecked(false);
	connect(eraseContourAct, SIGNAL(triggered()), this, SLOT(SetContourEraseFlag()));

	eraseSingleContourAct = new QAction(tr("Erase S"), drawActionGroup);
	eraseSingleContourAct->setIcon(QIcon(":erase_contour.png"));
	eraseSingleContourAct->setShortcut(QKeySequence(tr("Ctrl+W")));
	eraseSingleContourAct->setToolTip(tr("Erase Single Cone Contour (Ctrl+W)"));
	eraseSingleContourAct->setCheckable(true);
	eraseSingleContourAct->setChecked(false);
	connect(eraseSingleContourAct, SIGNAL(triggered()), this, SLOT(SetSingleContourEraseFlag()));

	redoAct = new QAction(tr("Redo"), this);
	redoAct->setIcon(QIcon(":redo.png"));
	redoAct->setShortcut(QKeySequence(tr("Ctrl+Z")));
	redoAct->setToolTip(tr("Redo Cone Operations (Ctrl+Z)"));
	redoAct->setEnabled(false);
	connect(redoAct, SIGNAL(triggered()), this, SLOT(RedoConeOperations()));

	emptyAct = new QAction(tr("Clear"), this);
	emptyAct->setIcon(QIcon(":trash.png"));
	emptyAct->setToolTip(tr("Clear All Mouse Operations"));
	connect(emptyAct, SIGNAL(triggered()), this, SLOT(RemoveInvisibleShapes()));

	aboutAct = new QAction(tr("About"), this);
	aboutAct->setIcon(QIcon(":about.png"));
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(ShowAboutDialog()));

	settingsAct = new QAction(tr("Settings"), this);
	settingsAct->setIcon(QIcon(":settings.png"));
	settingsAct->setToolTip("Image Display Settings");
	connect(settingsAct, SIGNAL(triggered()), SLOT(ShowSettingDialog()));

	toggleVisibilityAct = new QAction(tr("Toggle Visibility"), this);
	toggleVisibilityAct->setCheckable(true);
	toggleVisibilityAct->setChecked(true);
	toggleVisibilityAct->setShortcut(QKeySequence(tr("F2")));
	toggleVisibilityAct->setToolTip(tr("Toggle Segmented Cone Visibility (F2)"));
	connect(toggleVisibilityAct, SIGNAL(triggered()), this, SLOT(ToggleVisibility()));

	toggleInterpolationAct = new QAction(tr("Toggle Interpolation"), this);
	toggleInterpolationAct->setCheckable(true);
	toggleInterpolationAct->setChecked(true);
	toggleInterpolationAct->setShortcut(QKeySequence(tr("Ctrl+I")));
	toggleInterpolationAct->setToolTip(tr("Toggle Image Scale Pixel Interpolation (Ctrl+I)"));
	connect(toggleInterpolationAct, SIGNAL(triggered()), this, SLOT(ToggleInterpolation()));

	helpAct = new QAction(tr("Help"), this);
	helpAct->setShortcut(QKeySequence(tr("F1")));
	helpAct->setToolTip(tr("Display help screen (F1)"));
	helpAct->setIcon(QIcon(":help.png"));
	connect(helpAct, SIGNAL(triggered()), this, SLOT(ShowHelpWindow()));

	whatsNewAct = new QAction(tr("What's new?"), this);
	whatsNewAct->setIcon(QIcon(":help.png"));
	connect(whatsNewAct, SIGNAL(triggered()), this, SLOT(ShowWhatsNewWindow()));

	nextImageAct = new QAction(tr("Next Image"), this);
	nextImageAct->setShortcut(QKeySequence(tr("Down")));
	connect(nextImageAct, SIGNAL(triggered()), this, SLOT(OnNextImage()));
	prevImageAct = new QAction(tr("Previous Image"), this);
	prevImageAct->setShortcut(QKeySequence(tr("Up")));
	connect(prevImageAct, SIGNAL(triggered()), this, SLOT(OnPreviousImage()));
}

void radMainWindow::createMenu()
{
	fileMenu = menuBar()->addMenu(tr("&File"));
	SplitImageInputMenu = fileMenu->addMenu(tr("Split &Images"));
	SplitImageInputMenu->addAction(openSplitImageAct);

	saveMenu = fileMenu->addMenu(tr("Segmentation Results"));
	saveMenu->addAction(loadSegmentationAct);
	saveMenu->addAction(saveSegmentationAct);
	saveMenu->addAction(saveAllSegmentationsAct);

	fileMenu->addAction(settingsAct);
	fileMenu->addSeparator();
	fileMenu->addAction(quitAct);

	SplitMenu = menuBar()->addMenu(tr("&Split"));
	SplitMenu->addAction(segmentConesAct);
	SplitMenu->addAction(toggleVisibilityAct);
	SplitMenu->addAction(toggleInterpolationAct);
	SplitMenu->addSeparator();
	SplitMenu->addAction(purgeHistoryAct);

	helpMenu = menuBar()->addMenu(tr("&Help"));
	helpMenu->addAction(helpAct);
	helpMenu->addAction(whatsNewAct);
	helpMenu->addSeparator();
	helpMenu->addAction(aboutAct);
}

void radMainWindow::createView()
{
	QWidget *centralwidget = new QWidget;
	centralwidget->setAttribute(Qt::WA_DeleteOnClose, true);
    centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
    setCentralWidget(centralwidget);
    
    ImageWidget = new QVTKWidget(centralwidget);
	ImageWidget->SetRenderWindow(ImageView->GetRenderWin());

	SplitFileListWidget = new QListWidget();
	SplitFileListWidget->setToolTip("Image Files");
	SplitFileListWidget->setToolTipDuration(5000);
	SplitFileListWidget->setSelectionMode(QAbstractItemView::SingleSelection);
	connect(SplitFileListWidget, SIGNAL(currentItemChanged(QListWidgetItem*, QListWidgetItem*)), this, SLOT(SwitchSplitFile(QListWidgetItem*, QListWidgetItem*)));

	QVBoxLayout *FileListLayout = new QVBoxLayout;
	FileListLayout->addWidget(SplitFileListWidget, 4);
	
	QGridLayout *viewLayout = new QGridLayout(centralwidget);
	viewLayout->setGeometry(QRect(10, 40, 40, 300));
    viewLayout->setObjectName(QString::fromUtf8("viewLayout"));
	viewLayout->addWidget(ImageWidget, 0, 0);
	viewLayout->addLayout(FileListLayout, 0, 1);
	viewLayout->setColumnStretch(0, 5);
	viewLayout->setColumnStretch(1, 1);
}

void radMainWindow::createToolBar()
{
	QToolBar *drawToolBar = addToolBar(tr("&Draw"));
	drawToolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
	drawToolBar->addAction(openSplitImageAct);
	drawToolBar->addAction(saveSegmentationAct);
	drawToolBar->addAction(saveAllSegmentationsAct);
	drawToolBar->addAction(segmentConesAct);

	drawToolBar->addSeparator();
	drawToolBar->addAction(mouseAct);
	drawToolBar->addAction(drawContourAct);
	drawToolBar->addAction(editContourAct);
	drawToolBar->addAction(eraseContourAct);
	drawToolBar->addAction(eraseSingleContourAct);

	drawToolBar->addSeparator();
	drawToolBar->addAction(redoAct);
	drawToolBar->addAction(emptyAct);
	drawToolBar->addAction(settingsAct);

	drawToolBar->addSeparator();
	drawToolBar->addAction(helpAct);
}

void radMainWindow::createProgressDialog()
{
	progressDialog = new QProgressDialog("", "", 0, 100, this);
    Qt::WindowFlags flags = progressDialog->windowFlags();
    progressDialog->setWindowFlags(flags & ~(Qt::WindowCloseButtonHint |
		Qt::WindowMaximizeButtonHint | Qt::WindowMinimizeButtonHint |
		Qt::WindowSystemMenuHint | Qt::WindowContextHelpButtonHint));
	progressDialog->setCancelButton(0);
	progressDialog->reset();
	progressDialog->setWindowIcon(QIcon(":segment.png"));
	progressDialog->setWindowTitle("Cone Segmentation Progress");
	progressDialog->setMinimumSize(screen_width / 3, screen_height / 12);
	progressDialog->setMinimumDuration(10);
	connect(this, SIGNAL(sendFinishSegmentation()), this, SLOT(receiveFinishSegmentation()));
	connect(this, &radMainWindow::updateProgressText, progressDialog, &QProgressDialog::setLabelText);
}

void radMainWindow::createHelpWindow()
{
	helpWindow = new QWidget();
	helpWindow->setWindowIcon(QIcon(":help.png"));
	helpLayout = new QVBoxLayout();
	helpBrowser = new QTextBrowser();
	helpBrowser->setOpenExternalLinks(true);
	helpLayout->addWidget(helpBrowser);
	helpWindow->setLayout(helpLayout);

	QDir appDir(QCoreApplication::applicationDirPath());
	helpDir.setPath(appDir.absoluteFilePath(tr("Help")));
	helpFile = QFileInfo(helpDir.path(), tr("segment.html"));
	if (!helpFile.exists()) {
		appDir.cdUp();
		helpDir.setPath(appDir.absoluteFilePath(tr("Help")));
		helpFile = QFileInfo(helpDir.path(), tr("segment.html"));
	}

	helpWindow->setMinimumSize(screen_width * 55 / 100, screen_height * 50 / 100);
	helpWindow->move(screen_width * 20 / 100, screen_height * 25 / 100);
}

void radMainWindow::ShowAboutDialog()
{
	radAboutDialog dlg(this);
	dlg.exec();
}

void radMainWindow::ShowHelpWindow()
{
	helpWindow->setWindowTitle(tr("Help on Cone Segmentation"));
	if (helpFile.exists()) {
		helpBrowser->setSource(QUrl::fromLocalFile(helpFile.filePath()));
	}
	else {
		helpBrowser->setText("Sorry, no help available at this time.");
	}
	helpWindow->showNormal();
	helpWindow->activateWindow();
}

void radMainWindow::ShowWhatsNewWindow()
{
	QFileInfo whatsNewFile(helpDir.path(), tr("whatsnew.html"));
	if (!whatsNewFile.exists()) {
		return;
	}
	helpWindow->setWindowTitle(tr("What's new in Cone Segmentation"));
	helpBrowser->setSource(QUrl::fromLocalFile(whatsNewFile.filePath()));
	helpWindow->showNormal();
	helpWindow->activateWindow();
}

void radMainWindow::OnNextImage()
{
	if (SplitFileListWidget->count() == 0) return;
	int curRow = SplitFileListWidget->currentRow();
	if (curRow < 0) curRow = 0;
	else if (curRow + 1 < SplitFileListWidget->count()) ++curRow;
	if (curRow != SplitFileListWidget->currentRow()) {
		SplitFileListWidget->setCurrentRow(curRow);
		SplitFileListWidget->item(curRow)->setSelected(true);
	}
}
void radMainWindow::OnPreviousImage()
{
	if (SplitFileListWidget->count() == 0) return;
	int curRow = SplitFileListWidget->currentRow();
	if (curRow < 0) curRow = 0;
	else if (curRow > 0) --curRow;
	if (curRow != SplitFileListWidget->currentRow()) {
		SplitFileListWidget->setCurrentRow(curRow);
		SplitFileListWidget->item(curRow)->setSelected(true);
	}
}

void radMainWindow::openSplitImages(QStringList & fileNames, bool save_state)
{
	if (!GetVisibility()) SetVisibility(true);
	if (fileNames.isEmpty()) return;

	ClearSplitFileList();
	SplitMarkerInfor.resize(fileNames.size());
	qApp->processEvents(QEventLoop::ExcludeUserInputEvents);

	int j = 0;
	for (int i = 0; i < fileNames.size(); i++)
	{
		SplitMarkerInfor[i].Initialize();
		RGBImageType::Pointer res = FileIO->ReadSplitImage(fileNames[i].toStdString());
		if (res.GetPointer())
		{
			res = CreateTiffImage<RGBImageType>(res);
			SplitMarkerInfor[j].split_image = res;
			SplitMarkerInfor[j].split_file_names.first = fileNames[i].toStdString();
			SplitMarkerInfor[j].color_info.reset();
			// LoadBackupResults(j);
			++j;
		}
	}

	SplitMarkerInfor.resize(j);
	if (j == 0) return;

	UpdateSplitFileList();

	if (SplitMarkerInfor.size() > 0) {
		radBackup back_up;
		back_up.SetBackupDir(BackupDir);

		for (int id = 0; size_t(id) < SplitMarkerInfor.size(); id++) {
			if (!back_up.ReadBackup(SplitMarkerInfor[id])) {
				// First time -- try to open accompanying .json
				QFileInfo qimgfi(SplitMarkerInfor[id].split_file_names.first.c_str());
				QDir qimgdir = qimgfi.dir();
				QString jsonfn = qimgfi.completeBaseName() + ".json";
				QFileInfo qjsonfi(qimgdir, jsonfn);
				if (qjsonfi.exists()) {
					QFile qfi(qjsonfi.canonicalFilePath());
					ContourMarkersFromJSON(qfi, SplitMarkerInfor, true);
					back_up.WriteBackup(SplitMarkerInfor[id]);
				}
			}
		}
	}

	if (save_state && SplitMarkerInfor.size() > 0) {
		QFileInfo qimgfi(SplitMarkerInfor[0].split_file_names.first.c_str());
		saveDir = qimgfi.dir();
		saveState();
	}

	LoadSplitFile(0);
	SplitFileListWidget->setCurrentRow(CurrentImageIndex);
	SplitFileListWidget->item(CurrentImageIndex)->setSelected(true);
	BackupResults(CurrentImageIndex);
}

void radMainWindow::openSplitImage()
{
	QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::ExistingFiles);
	QStringList filters;
	filters << "TIFF (*.tif)";
	// filters << "PNG (*.png)";
	// filters << "JPEG (*.jpg)";
	// filters << "BMP (*.bmp)";
	dialog.setNameFilters(filters);
	dialog.setWindowTitle("Open Split Image");
	dialog.restoreState(fileDialogState);
	dialog.setDirectory(loadDir);

	QStringList fileNames; 
	if (dialog.exec())
	{
		fileNames = dialog.selectedFiles();
	}

    if ( !fileNames.isEmpty() )
    {

		saveDir = loadDir = dialog.directory();
		fileDialogState = dialog.saveState();
		saveState();

		openSplitImages(fileNames);
	}
}

void radMainWindow::loadSegmentations(QStringList & fileNames)
{
	if (SplitMarkerInfor.size() == 0) return;
	if (fileNames.isEmpty()) return;

	int firstRow = -1;
	// CurrentImageIndex
	for (int i = 0; i < fileNames.size(); i++) {
		QFileInfo fpath = QFileInfo(fileNames[i]);
        QFile qf(fpath.absoluteFilePath());
        int curRow = ContourMarkersFromJSON(qf, SplitMarkerInfor);
		if (curRow < 0) continue;
		if (firstRow < 0 || curRow == CurrentImageIndex)
			firstRow = curRow;
		BackupResults(curRow);
	}
	if (firstRow < 0) return;
	if (firstRow == CurrentImageIndex) {
		qApp->processEvents(QEventLoop::ExcludeUserInputEvents);
		ImageView->SetContourMarkers(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers);
		ImageView->ResetView(false);
		redoAct->setEnabled(false);
	}
	else {
		SplitFileListWidget->setCurrentRow(firstRow);
		SplitFileListWidget->item(firstRow)->setSelected(true);
	}
}

void radMainWindow::loadSegmentation()
{
	if (SplitMarkerInfor.size() == 0) return;

	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::ExistingFiles);
	dialog.setNameFilter(tr("Cone Segmentations (*.json)"));
	dialog.setWindowTitle("Open Cone Segmentations");
	dialog.restoreState(fileDialogState);
	dialog.setDirectory(saveDir);

	QStringList fileNames;
	if (dialog.exec()) {
		fileNames = dialog.selectedFiles();
		if (fileNames.isEmpty()) return;

		saveDir = dialog.directory();
		fileDialogState = dialog.saveState();
		saveState();

		loadSegmentations(fileNames);
	}
}

void radMainWindow::saveSegmentation()
{
	if (SplitFileListWidget->count() == 0 || SplitFileListWidget->currentRow() < 0
		|| SplitFileListWidget->currentRow() >= SplitFileListWidget->count())
		return;

	QString preferredName = QString::fromStdString(SplitMarkerInfor[SplitFileListWidget->currentRow()].split_file_names.second) + ".json";

	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::AnyFile);
	dialog.setAcceptMode(QFileDialog::AcceptSave);
	dialog.setWindowTitle("Save Cone Segmentations");
	dialog.selectFile(preferredName);
	dialog.selectNameFilter("Cone Segmentations (*.json)");
	dialog.restoreState(fileDialogState);
	dialog.setDirectory(saveDir);

	if (!dialog.exec()) return;

	QString filename = dialog.selectedFiles()[0];
	if (!filename.isNull()) {
		QFileInfo saveFile = QFileInfo(filename);
		saveDir.setPath(saveFile.dir().path());
		fileDialogState = dialog.saveState();
		saveState();
		SaveContourMarkers(SplitMarkerInfor[CurrentImageIndex], filename);
		UpdateSplitFileList(false);
	}
}

void radMainWindow::saveAllSegmentations()
{
	if (SplitFileListWidget->count() == 0)
		return;

	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::Directory);
	dialog.setWindowTitle("Select Output Directory");
	dialog.setOptions(QFileDialog::DontUseNativeDialog | QFileDialog::DontResolveSymlinks);
	// ! Non-native dialog: shows files, but save/restore state does not work
	// dialog.restoreState(fileDialogState);
	dialog.setMinimumSize(screen_width * 3 / 7, screen_height * 2 / 5);
	dialog.setDirectory(saveDir);

	if (!dialog.exec()) return;

	QString dir = dialog.selectedFiles()[0];
	if (dir.isNull()) return;

	saveDir.setPath(dir);
	// fileDialogState = dialog.saveState();
	saveState();

	vector<pair<QFileInfo, size_t>> todo;
	string existing = "";
	int cnt = 0;

	for (int row = 0; row < SplitFileListWidget->count(); row++) {
		QString qfn = QString::fromStdString(SplitMarkerInfor[row].split_file_names.second) + ".json";
		QFileInfo saveFile = QFileInfo(saveDir, qfn);
		if (saveFile.exists()) {
			++cnt;
			if (cnt <= 10) {
				if (existing.size() > 0) existing = existing + "\n";
				existing = existing + qfn.toStdString();
			}
		}
		else if (!SplitMarkerInfor[row].cone_contour_markers.empty()) {
			todo.push_back(std::make_pair(saveFile, (size_t)row));
		}
	}

	if (existing.size() > 0) {
		QMessageBox msgBox;
		string msg = "The following file(s) already exist(s):\n";
		msg = msg + existing;
		msgBox.setText(msg.c_str());
		msgBox.setInformativeText("Do you want to write over existing files?");
		msgBox.setIcon(QMessageBox::Warning);
		msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
		msgBox.setDefaultButton(QMessageBox::Yes);
		if (msgBox.exec() == QMessageBox::Yes)
			existing.clear();
	}
	if (existing.size() == 0) {
		for (std::pair<QFileInfo, size_t> & item : todo) {
			QFileInfo fn = item.first;
			SaveContourMarkers(SplitMarkerInfor[(size_t)item.second], fn.absoluteFilePath());
		}
		UpdateSplitFileList(false);
	}
}

void radMainWindow::purgeHistoryFiles()
{
	purgeHistoryDialog->showHistory(BackupDir);
}

void radMainWindow::LoadSplitFile(int id)
{
	CurrentImageIndex = id;
	ImageView->InitializeView();
	ImageView->SetSplitImage(SplitMarkerInfor[id].split_image);
	ImageView->SetColorInfo(SplitMarkerInfor[id].color_info);
	ImageView->SetContourMarkers(SplitMarkerInfor[id].cone_contour_markers);
	ImageView->ResetView();
}

void radMainWindow::SetDetectionCategory(int id)
{
	DetectionCategory = id;
}

void radMainWindow::quit()
{
    close();
}

void radMainWindow::UpdateSplitFileList(bool newlist)
{
	if (size_t(SplitFileListWidget->count()) != SplitMarkerInfor.size())
		newlist = true;
	if (newlist) {
		SplitFileListWidget->clear();
		for (int i = 0; i < SplitMarkerInfor.size(); i++)
		{
			QFileInfo qfi(SplitMarkerInfor[i].split_file_names.first.c_str());
			QString basename = qfi.completeBaseName();
			SplitMarkerInfor[i].split_file_names.second = basename.toStdString();
			SplitFileListWidget->addItem(new QListWidgetItem(GetListName(SplitMarkerInfor[i].split_file_names.first)));
		}
	}
	else {
		for (int i = 0; i < SplitMarkerInfor.size(); i++) {
			SplitFileListWidget->item(i)->setText(GetListName(SplitMarkerInfor[i].split_file_names.first));
		}
	}
}

void radMainWindow::ClearSplitFileList()
{
	SplitMarkerInfor.clear();
	SplitFileListWidget->clear();
}

void radMainWindow::SwitchSplitFile(QListWidgetItem *item, QListWidgetItem *previous)
{
	int i;

	if (previous) {
		std::string prev = previous->text().mid(1).toStdString();
		for (i = 0; i < SplitMarkerInfor.size(); i++) {
			if (SplitMarkerInfor.at(i).split_file_names.second.compare(prev) == 0) {
				SplitMarkerInfor[i].color_info = ImageView->GetColorInfo();
				BackupResults(i);
			}
		}
	}
	if (item) {
		std::string curr = item->text().mid(1).toStdString();
		for (i = 0; i < SplitMarkerInfor.size(); i++)
		{
			if (SplitMarkerInfor.at(i).split_file_names.second.compare(curr) == 0)
			{
				LoadSplitFile(i);
			}
		}

		qDebug() << item->text();
	}

}

void radMainWindow::BackupResults(int id)
{
	if (id < 0 || id >= (int)SplitMarkerInfor.size()) return;
	radBackup back_up;
	back_up.SetBackupDir(BackupDir);
	back_up.WriteBackup(SplitMarkerInfor[id]);
}

void radMainWindow::LoadBackupResults(int id)
{
	if (id < 0 || id >= (int)SplitMarkerInfor.size()) return;
	radBackup back_up;
	back_up.SetBackupDir(BackupDir);
	back_up.ReadBackup(SplitMarkerInfor[id]);
}

void radMainWindow::SetContourEraseFlag()
{
	MouseOperationType = Mouse_Delete_Contour;

	EditedContourId = -1;
	ImageView->DisableEditedContour();
	ImageView->ResetView(false);
}

void radMainWindow::SetSingleContourEraseFlag()
{
	MouseOperationType = Mouse_Delete_Single_Contour;
	EditedContourId = -1;
	ImageView->DisableEditedContour();
	ImageView->ResetView(false);
}

void radMainWindow::SetContourDrawingFlag()
{
	MouseOperationType = Mouse_Add_Contour;

	EditedContourId = -1;
	ImageView->DisableEditedContour();
	ImageView->ResetView(false);
}

void radMainWindow::SetContourEditFlag()
{
	MouseOperationType = Mouse_Edit_Contour;
	EditedContourId = -1;
}

void radMainWindow::SetMouseFlag()
{
	MouseOperationType = Mouse_Normal;

	EditedContourId = -1;
	ImageView->DisableEditedContour();
	ImageView->ResetView(false);
}

void radMainWindow::AddConeContours(DoublePointArray & contour_pts)
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size()
		|| contour_pts.size() < 3)
		return;

	radContourOperations::AddConeContour(contour_pts, SplitMarkerInfor[CurrentImageIndex]);
	ImageView->SetContourMarkers(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers);
	if (SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.empty())
		return;

	//update stack information
	SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.push(StackInformation(Shape_Add, 
		SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.size()-1));
	redoAct->setEnabled(!SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.empty());

	BackupResults(CurrentImageIndex);
}

void radMainWindow::RemoveConeContours(DoublePointArray &contour_pts)
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size()
		|| contour_pts.size() < 3)
		return;

	vector<int> deleted_contour_list;
	radContourOperations::RemoveConeContour(contour_pts, SplitMarkerInfor[CurrentImageIndex], 
		deleted_contour_list);

	//update image view
	if (!deleted_contour_list.empty())
	{
		ImageView->SetContourMarkers(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers);
		//update stack information
		SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.push(StackInformation(Shape_Delete,
			deleted_contour_list));
		redoAct->setEnabled(!SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.empty());
		BackupResults(CurrentImageIndex);
	}
}

void radMainWindow::EditConeContours(double xpos, double ypos, double zpos)
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size())
		return;

	if (SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.empty())
		return;

	DoublePointType pt;
	pt[0] = xpos; pt[1] = ypos;
	EditedContourId = radContourOperations::GetEditedConeContour(pt, SplitMarkerInfor[CurrentImageIndex]);

	if (EditedContourId != -1)
	{
		if (!SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[EditedContourId].is_edited)
		{
			radContourOperations::UpdateEditedConeContour(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[EditedContourId].marker_contours);
			SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[EditedContourId].is_edited = true;
		}
		ImageView->SetContourMarkers(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers);
		ImageView->EnableEditedContour(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[EditedContourId].marker_contours);
	}
	else {
		ImageView->DisableEditedContour();
	}
	ImageView->ResetView(false);
}

void radMainWindow::RemoveSingleConeContour(double xpos, double ypos, double zpos)
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size())
		return;

	DoublePointType deleted_pt;
	deleted_pt[0] = xpos;
	deleted_pt[1] = ypos;

	int res_id = radContourOperations::GetEditedConeContour(deleted_pt, SplitMarkerInfor[CurrentImageIndex]);

	//update stack information
	if (res_id != -1)
	{
		if (res_id == EditedContourId) {
			ImageView->DisableEditedContour();
			ImageView->ResetView(false);
		}
		vector<int> deleted_contour_list;
		deleted_contour_list.push_back(res_id);
		SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[res_id].is_visible = false;

		ImageView->SetContourMarkers(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers);
		//update stack information
		SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.push(StackInformation(Shape_Delete,
			deleted_contour_list));
		redoAct->setEnabled(!SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.empty());
		BackupResults(CurrentImageIndex);
	}
}

void radMainWindow::UpdateConeContours(vtkSmartPointer<vtkPolyData> poly)
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size())
		return;

	if (SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.empty() || EditedContourId == -1)
		return;

	SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[EditedContourId].marker_contours.clear();
	for (unsigned int i=0; i<poly->GetNumberOfPoints(); i++)
	{
		DoublePointType pt;
		pt[0] = poly->GetPoint(i)[0];
		pt[1] = poly->GetPoint(i)[1];
		SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[EditedContourId].marker_contours.push_back(pt);
	}
	SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[EditedContourId].ComputeShapeInfor();
	SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[EditedContourId].ComputeShapeRegion<RGBImageType>(
		SplitMarkerInfor[CurrentImageIndex].split_image);
	ImageView->SetContourMarkers(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers);

	BackupResults(CurrentImageIndex);
}

void radMainWindow::PushColorUndo(ColorInfo ci)
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size())
		return;

	SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.push(StackInformation(ci));
	redoAct->setEnabled(!SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.empty());
}

void radMainWindow::RemoveInvisibleShapes()
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size())
		return;

	vector<ContourMarker>::iterator it = SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.end();
	QApplication::setOverrideCursor(Qt::WaitCursor);

	//update here...
	vector<unsigned int> deleted_contours;
	int index = SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.size();
	while (it > SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.begin())
	{
		it--;
		index = std::distance(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.begin(), it);

		if (!it->is_visible)
		{
			it = SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.erase(it); //}
		}

		qApp->processEvents(QEventLoop::ExcludeUserInputEvents);
	}

	QApplication::restoreOverrideCursor();

	SplitMarkerInfor[CurrentImageIndex].ClearStack(SplitMarkerInfor[CurrentImageIndex].contour_operator_stack);
	redoAct->setEnabled(!SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.empty());

	BackupResults(CurrentImageIndex);
}

void radMainWindow::RedoConeOperations()
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size())
		return;

	if (SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.empty()
		|| SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.empty())
		return;

	StackInformation cur_stack_infor = SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.top();
	SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.pop();

	if (cur_stack_infor.shape_operator == Shape_Delete)
	{
		for (unsigned int i = 0; i<cur_stack_infor.shape_list.size(); i++)
		{
			if (cur_stack_infor.shape_list[i] >= 0 && cur_stack_infor.shape_list[i]
				< SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.size())
				SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[cur_stack_infor.shape_list[i]].is_visible = true;
		}
	}
	else if (cur_stack_infor.shape_operator == Shape_Add)
	{
		for (unsigned int i = 0; i<cur_stack_infor.shape_list.size(); i++)
		{
			if (cur_stack_infor.shape_list[i] >= 0 && cur_stack_infor.shape_list[i]
				< SplitMarkerInfor[CurrentImageIndex].cone_contour_markers.size())
				SplitMarkerInfor[CurrentImageIndex].cone_contour_markers[cur_stack_infor.shape_list[i]].is_visible = false;
		}
	}
	else if (cur_stack_infor.shape_operator == Image_Op) {
		ImageView->SetColorInfo(cur_stack_infor.ci);
	}

	redoAct->setEnabled(!SplitMarkerInfor[CurrentImageIndex].contour_operator_stack.empty());
	ImageView->SetContourMarkers(SplitMarkerInfor[CurrentImageIndex].cone_contour_markers);
	
	ImageView->ResetView(false);
	BackupResults(CurrentImageIndex);
}

void radMainWindow::ShowSettingDialog()
{
	if (!GetVisibility()) SetVisibility(true);
	SettingsDlg->updateSettings();
	SettingsDlg->show();
	WriteSystemSettings();
}

void radMainWindow::ToggleVisibility()
{
	SetVisibility(GetVisibility());
}

void radMainWindow::ToggleInterpolation()
{
	ImageView->SetInterpolation(toggleInterpolationAct->isChecked());
	ImageView->ResetView(false);
	saveState();
}

void radMainWindow::SetVisibility(bool flag)
{
	toggleVisibilityAct->setChecked(flag);
	if (flag) {
		ImageView->SetConeContourVisibility(SystemSettings->visibility_contour);
		ImageView->SetConeCenterVisibility(SystemSettings->visibility_center);
		ImageView->SetConeRegionVisibility(SystemSettings->visibility_region);
	}
	else {
		ImageView->SetConeContourVisibility(false);
		ImageView->SetConeCenterVisibility(false);
		ImageView->SetConeRegionVisibility(false);
	}
	ImageView->ResetView(false);
}

void radMainWindow::SegmentConesChecked(QList<int> checked)
{
	checkedItems = checked;
	if (checked.size() == 0) return;
	if (!GetVisibility()) SetVisibility(true);
	progressDialog->reset();
	progressDialog->setMaximum(8 * checkedItems.size());
	QThread *thread = QThread::create([this] {
		for (int id : checkedItems) {
			SegmentCones(id);
		}
		emit sendFinishSegmentation();
	});
	thread->start();
	progressDialog->exec();
}

void radMainWindow::tac()
{
	progressDialog->setValue(progressDialog->value() + 1);
}

void radMainWindow::receiveFinishSegmentation()
{
	progressDialog->reset();

	if (checkedItems.size() > 0 && checkedItems.indexOf(CurrentImageIndex) < 0) {
		CurrentImageIndex = checkedItems[0];
		SplitFileListWidget->setCurrentRow(CurrentImageIndex);
	}

	BackupResults(CurrentImageIndex);

	qApp->processEvents(QEventLoop::ExcludeUserInputEvents);
	MarkerInformation & split_infor = SplitMarkerInfor[CurrentImageIndex];
	ImageView->SetContourMarkers(split_infor.cone_contour_markers);
	ImageView->ResetView(false);
	redoAct->setEnabled(false);
}

void radMainWindow::showSegmentationPanel()
{
	if (CurrentImageIndex < 0 || CurrentImageIndex >= SplitMarkerInfor.size())
		return;
	EditedContourId = -1;
	ImageView->DisableEditedContour();

	QStringList items;
	QList<int> checked;
	for (int row = 0; size_t(row) < SplitMarkerInfor.size(); row++) {
		LoadBackupResults(row);
		MarkerInformation & split_infor = SplitMarkerInfor[row];
		items << split_infor.split_file_names.second.c_str();
		if (split_infor.cone_contour_markers.size() == 0)
			checked << row;
	}
	SegmentationPanel->SetItemList(items);
	SegmentationPanel->SetCheckedRows(checked);
	SegmentationPanel->SetHighlightedRow(CurrentImageIndex);

	SegmentationPanel->SetParameters(SplitMarkerInfor[CurrentImageIndex].segment_params);
	SegmentationPanel->show();
}

void radMainWindow::SegmentCones(int id)
{
	if (id < 0 || id >= SplitMarkerInfor.size()) return;

	MarkerInformation & split_infor = SplitMarkerInfor[id];
	split_infor.segment_params = SegmentationPanel->GetParameters();

	emit updateProgressText(QString::fromStdString(split_infor.split_file_names.second));

	radConeSegmentation cone_segmentation;
	connect(&cone_segmentation, SIGNAL(tick()), this, SLOT(tac()));

	cone_segmentation.SetParameters(split_infor.segment_params);
	
	cone_segmentation.SegmentRgbImage(split_infor.split_image);

	std::vector<DoublePointArray> & contours = cone_segmentation.GetSegmentationContours();

	split_infor.ClearStack(split_infor.contour_operator_stack);
	split_infor.cone_contour_markers.resize(contours.size());
	unsigned ii = 0;
	itk::ImageRegion<2U> rgn = split_infor.split_image->GetLargestPossibleRegion();
	for (unsigned int i = 0; i < split_infor.cone_contour_markers.size(); i++)
	{
		ContourMarker & marker = split_infor.cone_contour_markers[ii];
		marker.is_visible = true;
		DoublePointArray & src = contours[i];
		DoublePointArray & tgt = marker.marker_contours;
		tgt.resize(src.size());
		for (unsigned int j = 0; j < src.size(); j++) {
			tgt[j] = src[j];
		}
		marker.ComputeShapeInfor();
		if (!marker.withinImageRegion(rgn))
			// skip markers that cross image boundaries
			continue;
		marker.ComputeShapeRegion<RGBImageType>(split_infor.split_image);
		ii++;
	}
	split_infor.cone_contour_markers.resize(ii);

	disconnect(&cone_segmentation, SIGNAL(tick()), 0, 0);

}

void radMainWindow::ReadSystemSettings()
{
	QFile fi(SegmSystemFile.filePath());
	if (fi.open(QIODevice::ReadOnly)) {
		QJsonDocument json = QJsonDocument::fromJson(fi.readAll());
		fi.close();

		QJsonObject jobj = json.object();
		// cout << json.toJson(QJsonDocument::Indented).toStdString().c_str() << std::endl;

		QJsonObject jobj_coloring = jobj["Coloring"].toObject();
		if (jobj_coloring["Contour"].isString())
			SystemSettings->contour_color = jobj_coloring["Contour"].toString();
		if (jobj_coloring["Region"].isString())
			SystemSettings->region_color = jobj_coloring["Region"].toString();
		if (jobj_coloring["Center"].isString())
			SystemSettings->center_color = jobj_coloring["Center"].toString();

		QJsonObject jobj_visiblity = jobj["Visibility"].toObject();
		SystemSettings->visibility_contour = (bool) jobj_visiblity["Contour"].toInt();
		SystemSettings->visibility_region = (bool) jobj_visiblity["Region"].toInt();
		SystemSettings->visibility_center = (bool) jobj_visiblity["Center"].toInt();
		if (jobj_visiblity["Opacity"].isDouble())
			SystemSettings->region_opacity = jobj_visiblity["Opacity"].toDouble();

		QJsonObject jobj_size = jobj["Size"].toObject();
		SystemSettings->size_contour = jobj_size["Contour"].toInt();
		if (jobj_size["Center"].isDouble())
			SystemSettings->size_center = jobj_size["Center"].toDouble();
	}
}

void radMainWindow::WriteSystemSettings()
{
	QJsonObject jobj;

	QJsonObject jobj_visiblity;
	jobj_visiblity["Contour"] = QJsonValue((int) SystemSettings->visibility_contour);
	jobj_visiblity["Center"] = QJsonValue((int) SystemSettings->visibility_center);
	jobj_visiblity["Region"] = QJsonValue((int) SystemSettings->visibility_region);
	jobj_visiblity["Opacity"] = QJsonValue(SystemSettings->region_opacity);
	jobj["Visibility"] = jobj_visiblity;

	QJsonObject jobj_size;
	jobj_size["Contour"] = SystemSettings->size_contour;
	jobj_size["Center"] = SystemSettings->size_center;
	jobj["Size"] = jobj_size;

	QJsonObject jobj_coloring;
	jobj_coloring["Contour"] = SystemSettings->contour_color.name();
	jobj_coloring["Center"] = SystemSettings->center_color.name();
	jobj_coloring["Region"] = SystemSettings->region_color.name();
	jobj["Coloring"] = jobj_coloring;

	QJsonDocument json = QJsonDocument(jobj);

	// cout << json.toJson(QJsonDocument::Indented).toStdString().c_str() << std::endl;
	QFile fo(SegmSystemFile.filePath());
	if (fo.open(QIODevice::WriteOnly)) {
		fo.write(json.toJson(QJsonDocument::Indented));
		fo.close();
	}
}

void radMainWindow::saveState() {
	QJsonObject jobj;
	jobj["loadDir"] = loadDir.path();
	jobj["saveDir"] = saveDir.path();
	jobj["fileDialogState"] = QString(fileDialogState.toBase64());
	jobj["interpolation"] = ImageView->GetInterpolation();
	jobj["version"] = QString(Segm_VERSION.c_str());
	QJsonDocument json = QJsonDocument(jobj);

	// cout << json.toJson(QJsonDocument::Indented).toStdString().c_str() << std::endl;
	QFile fo(SegmStateFile.filePath());
	if (fo.open(QIODevice::WriteOnly)) {
		fo.write(json.toJson(QJsonDocument::Indented));
		fo.close();
	}
}

void radMainWindow::loadState() {
	QFile fi(SegmStateFile.filePath());
	if (fi.open(QIODevice::ReadOnly)) {
		QJsonDocument json = QJsonDocument::fromJson(fi.readAll());
		fi.close();

		QJsonObject jobj = json.object();
		// cout << json.toJson(QJsonDocument::Indented).toStdString().c_str() << std::endl;
		if (jobj["loadDir"].isString())
			loadDir.setPath(jobj["loadDir"].toString());
		if (jobj["saveDir"].isString())
			saveDir.setPath(jobj["saveDir"].toString());
		if (jobj["fileDialogState"].isString())
			fileDialogState = QByteArray::fromBase64(QByteArray(jobj["fileDialogState"].toString().toStdString().c_str()));
		if (jobj["interpolation"].isBool()) {
			ImageView->SetInterpolation(jobj["interpolation"].toBool());
			toggleInterpolationAct->setChecked(ImageView->GetInterpolation());
		}
		if (jobj["version"].isString())
			lastVersion = jobj["version"].toString();
	}
}

static long long decode_version(const char* ver)
{
	int mj, mn, mc;
	if (sscanf(ver, "%d.%d.%d ", &mj, &mn, &mc) != 3)
		return 0L;
	return (long long)(mj) * 1000000000L + (long long)(mn) * 1000000L + (long long)(mc);
}

void radMainWindow::checkWhatsNew()
{
	long long cur_ver = decode_version(Segm_VERSION.c_str());
	long long old_ver = decode_version(lastVersion.toStdString().c_str());
	// std::cout << "cur_ver = " << cur_ver << " ; old_ver = " << old_ver << std::endl;
	if (cur_ver != old_ver)
		saveState();
	if (cur_ver > old_ver)
		ShowWhatsNewWindow();
}
