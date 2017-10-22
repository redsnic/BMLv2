#ifndef MAFTRANSLATOR_H
#define MAFTRANSLATOR_H

#include <QWidget>
#include <QFileDialog>
#include "MAFtoMatrixTranslator.h"
#include <QRegExp>
#include <QMessageBox>

namespace Ui {
class MAFtranslator;
}

class MAFtranslator : public QWidget
{
    Q_OBJECT

public slots:
    void getInput();
    void setOutput();
    void translate();
    void setStandard(bool);
    void setAll(bool);
    void setAllGenes();



public:
    explicit MAFtranslator(QWidget *parent = 0);
    ~MAFtranslator();

signals:
    void closeEvent();

private:
    Ui::MAFtranslator *ui;
    void closeEvent(QCloseEvent *event);
};

#endif // MAFTRANSLATOR_H
