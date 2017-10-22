#ifndef MAINMENU_H
#define MAINMENU_H

#include <QMainWindow>
#include <QFileDialog>
#include <fstream>
#include <QFuture>
#include <QUrl>
#include <QFile>
#include <QDesktopServices>
#include <QtConcurrentRun>
#include <QtConcurrent>
#include <QFutureWatcher>
#include <QApplication>
#include "helpviewer.h"
#include "Algorithm.h"
#include "maftranslator.h"
#include <QtWebEngine>
#ifndef DS
#include "DataStructures.h"
#define DS
#endif



namespace Ui {
class MainMenu;
}

class MainMenu : public QMainWindow
{
    Q_OBJECT


public slots:
    void openHelp_About_BML();
    void openHelp_BML_publication();
    void openHelp_Graphviz();
    void openHelp_Guide();
    void openHelp_Changelog();
    void openHelp_cBioPortal();
    void browseMatrixPressed();
    void cutoffAutoSet(bool);
    void bootstrapEnabledSet(bool);
    void runBML();
    void runMAF();
    void MAFtranslationComplete();
    void executeFinished();
    void wlog(std::string);

public:
    void log(std::string);
    explicit MainMenu(QWidget *parent = 0);
    bool file_exist(std::string fileName);
    ~MainMenu();

private:
    Ui::MainMenu *ui;
    void openHelp(std::string);
    bool inputValuesPreliminaryCheck();


};

#endif // MAINMENU_H
