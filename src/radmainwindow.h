#ifndef RADMAINWINDOW_H
#define RADMAINWINDOW_H

#include <QMainWindow>
#include <QProgressDialog>
#include <QtGui>
#include "QVTKWidget.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QCheckBox>
#include <QLineEdit>
#include <QDoubleValidator>
#include <QTreeWidget>
#include "radImageView.h"
#include "radFileIO.h"
#include "QActionGroup.h"
#include "radBackup.h"
#include <QAction>
#include <QLabel>
#include <QMenu>
#include <QGridLayout>
#include <QWidgetAction>
#include <QListWidgetItem>
#include <QListWidget>
#include <QSpinBox>
#include <QMessageBox>

#include "radAboutDialog.h"
#include "radSettingsDialog.h"
#include "radSegmentationPanel.h"
#include "radconesegmentation.h"
#include "radPurgeHistoryDialog.h"

void decodeFileList(QStringList &inputNames, QStringList &splitFileNames, QStringList &detectionFileNames);

class QVTKWidget;
class radMainWindow : public QMainWindow
{
  Q_OBJECT
public:

    // Constructor/Destructor
    radMainWindow(); 
    ~radMainWindow();
    
    static radMainWindow* GetPointer() { return TheApp; }
	radImageView * GetImageView() {return ImageView;}
	radFileIO * GetFileIO() {return FileIO;}

	void saveState();
	void loadState();

	void SetVisibility(bool flag = true);
	bool GetVisibility() { return toggleVisibilityAct->isChecked(); }

	void openSplitImages(QStringList & fileNames, bool save_state=false);
	void loadSegmentations(QStringList & fileNames);
	void NextImage() { OnNextImage(); }
	void PreviousImage() { OnPreviousImage(); }

	void AddConeContours(DoublePointArray &);
	void RemoveConeContours(DoublePointArray &);
	void RemoveSingleConeContour(double, double, double);

	void EditConeContours(double, double, double);
	void UpdateConeContours(vtkSmartPointer<vtkPolyData>);
	void PushColorUndo(ColorInfo ci);

	MouseOperation MouseOperationType;

	void checkWhatsNew();

signals:
	void sendFinishSegmentation();
	void updateProgressText(QString text);

private slots:

	void ShowAboutDialog();
	void ShowHelpWindow();
	void ShowWhatsNewWindow();
	void openSplitImage();
	void loadSegmentation();
	void saveSegmentation();
	void saveAllSegmentations();
	void quit();

	void SwitchSplitFile(QListWidgetItem*, QListWidgetItem*);
	void SetDetectionCategory(int);

	void SetMouseFlag();
	void SetContourDrawingFlag();
	void SetContourEraseFlag();
	void SetSingleContourEraseFlag();
	void SetContourEditFlag();
	
	void RedoConeOperations();
	void RemoveInvisibleShapes();
	void ShowSettingDialog();
	void ToggleVisibility();
	void ToggleInterpolation();
	void ToggleVoronoi();

	void showSegmentationPanel();
	void purgeHistoryFiles();
	void SegmentConesChecked(QList<int> checked);
	void tac();
	void receiveFinishSegmentation();

	void OnNextImage();
	void OnPreviousImage();

protected:
	virtual void closeEvent(QCloseEvent *event);
	void dragEnterEvent(QDragEnterEvent * event);
	void dropEvent(QDropEvent *e);

private:
	int screen_width, screen_height;
    
    static radMainWindow* TheApp;

	QDir SegmHomeDir;
	QFileInfo SegmStateFile;
	QFileInfo SegmSystemFile;
	QByteArray fileDialogState;
	QDir saveDir;
	QDir loadDir;
	string BackupDir;
	QList<int> checkedItems;

	QDir helpDir;
	QFileInfo helpFile;
	QString lastVersion;

	QAction * openSplitImageAct;
	QAction * loadSegmentationAct;
	QAction * saveSegmentationAct;
	QAction * saveAllSegmentationsAct;
	QAction * quitAct;

	QAction *nextImageAct;
	QAction *prevImageAct;

	QAction * segmentConesAct;
	QAction * purgeHistoryAct;
	QActionGroup * drawActionGroup;
	QAction * mouseAct;
	QAction * drawContourAct;
	QAction * editContourAct;
	QAction * eraseContourAct;
	QAction * eraseSingleContourAct;
	QAction* redoAct;
	QAction* voronoiAct;
	QAction * emptyAct;
	QAction * settingsAct;
	QAction * toggleVisibilityAct;
	QAction * toggleInterpolationAct;

	QAction* aboutAct;
	QAction* helpAct;
	QAction* whatsNewAct;

	QMenu * fileMenu;
	QMenu * SplitImageInputMenu;
	QMenu * saveMenu;
	QMenu * SplitMenu;
	QMenu * helpMenu;

	void createActions();
    void createMenu();
    void createView();
	void createToolBar();
	void createProgressDialog();
	void createHelpWindow();

    QVTKWidget *ImageWidget;
	QListWidget *SplitFileListWidget;
	radImageView *ImageView;

	radSegmentationPanel *SegmentationPanel;
	radPurgeHistoryDialog *purgeHistoryDialog;

	radFileIO *FileIO;
	QProgressDialog * progressDialog;
	radSettingsDialog *SettingsDlg;
	int CurrentImageIndex;
	int EditedContourId;
	int DetectionCategory;
	vector< MarkerInformation > SplitMarkerInfor;
	MarkerSystemSettings *SystemSettings;
	
	QWidget *helpWindow;
	QVBoxLayout *helpLayout;
	QTextBrowser *helpBrowser;

	void UpdateSplitFileList(bool newlist = true);
	void ClearSplitFileList();
	void LoadSplitFile(int);

	void SegmentCones(int id);
	
	void BackupResults(int);
	void LoadBackupResults(int);
	void ReadSystemSettings();
	void WriteSystemSettings();
};

#endif
