#include "mainmenu.h"
#include <QApplication>


int main(int argc, char *argv[])
{
    if(!QDir("output").exists()){  // creates output dir if necessary
        QDir().mkdir("output");
    }
    QApplication a(argc, argv);
    MainMenu w;
    w.show();

    return a.exec();
}
