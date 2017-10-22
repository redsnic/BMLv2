#include "helpviewer.h"
#include "ui_helpviewer.h"

helpViewer::helpViewer(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::helpViewer)
{
    ui->setupUi(this);
    connect( this->ui->webView->page()->profile() , SIGNAL(downloadRequested(QWebEngineDownloadItem*)),
                this, SLOT(downloader(QWebEngineDownloadItem*)));  // enables download request's detections
    this->setAttribute( Qt::WA_DeleteOnClose );
}

helpViewer::~helpViewer()
{
    disconnect( ui->webView->page()->profile() , SIGNAL(downloadRequested(QWebEngineDownloadItem*)),
                this, SLOT(downloader(QWebEngineDownloadItem*)));  // enables download request's detections
    delete ui;
}

/**
 * @brief helpViewer::show
 * Display the site with url encoded in a standard library string
 * @param url a standard url written in a std::string
 */
void helpViewer::show(std::string url){
    ui->webView->load(QUrl(url.data()));
}

/**
 * @brief helpViewer::downloadRequested
 * manages the download requests by asking the user where to put the downloadable file.
 * @param download downloadable content
 */
void helpViewer::downloader(QWebEngineDownloadItem* download) {

    if(!isActiveWindow()){return;}  // manage download request only for the active window

    QString selected = QFileDialog::getSaveFileName(this, tr("Choose download location:"), download->path() );
    if(selected.length() != 0){
        download->setPath(selected);
        qDebug() << "Format: " <<  download->savePageFormat();  // donwload settings
        qDebug() << "Path: " << download->path();
        download->accept();
    } else {
        download->cancel();
    }

}
