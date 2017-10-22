#ifndef HELPVIEWER_H
#define HELPVIEWER_H

#include <QWidget>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QWebEngineDownloadItem>
#include <QUrl>

namespace Ui {
class helpViewer;
}

class helpViewer : public QWidget
{
    Q_OBJECT
public slots:
    void downloader(QWebEngineDownloadItem* download);

public:
    explicit helpViewer(QWidget *parent = 0);
    void show(std::string);
    ~helpViewer();

private:
    Ui::helpViewer *ui;
    QString lastServed;

};

#endif // HELPVIEWER_H
