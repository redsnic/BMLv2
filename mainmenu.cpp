#include "mainmenu.h"
#include "ui_mainmenu.h"

MainMenu::MainMenu(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainMenu)
{
    ui->setupUi(this);
    /* get the threads' number of the host machine */
    ui->numberOfThreadsLine->setText(QString::number(QThread::idealThreadCount()));
    QtWebEngine::initialize();
}

MainMenu::~MainMenu()
{
    delete ui;
}

/**
 * @brief MainMenu::runMAF
 *
 * this porcedure creates the window for the MAF preprocessor
 */
void MainMenu::runMAF(){
    ui->mafButton->setText("Operation in progress...");
    ui->mafButton->setEnabled(false);
    MAFtranslator* tr = new MAFtranslator;
    connect(tr, SIGNAL(closeEvent()), SLOT(MAFtranslationComplete()));
    tr->setVisible(true);
}

/**
 * @brief MainMenu::MAFtranslationComplete
 *
 * restores the button to call the MAF preprocessor when the previously opened
 * window is closed.
 */
void MainMenu::MAFtranslationComplete(){
    ui->mafButton->setText("Run MAF file preprocessor utility");
    ui->mafButton->setEnabled(true);
}

/**
 * @brief MainMenu::openHelp
 * opens a window to display the help informations from a local or
 * remote web page
 * @param url        url of the help page
 */
void MainMenu::openHelp(std::string url){
    helpViewer* help = new helpViewer();
    help->setVisible(true);
    help->show(url);
}

/* Help pages */

void MainMenu::openHelp_About_BML(){
    openHelp("http://bml.molgen.mpg.de/");
}

void MainMenu::openHelp_BML_publication(){
    openHelp("https://academic.oup.com/bioinformatics/article/30/17/2456/2748189");
}

void MainMenu::openHelp_Graphviz(){
    openHelp("https://dreampuf.github.io/GraphvizOnline/");
}

void MainMenu::openHelp_cBioPortal(){
    openHelp("http://www.cbioportal.org/data_sets.jsp");
}

void MainMenu::openHelp_Guide(){
    openHelp(QUrl::fromLocalFile( QFileInfo("BML_Guide.html").absoluteFilePath() ).toString().toStdString());
}

void MainMenu::openHelp_Changelog(){
    openHelp(QUrl::fromLocalFile( QFileInfo("Changelog.html").absoluteFilePath() ).toString().toStdString());
}

/**
 * @brief MainMenu::browseMatrixPressed
 * opens dialog to select an input data matrix
 */
void MainMenu::browseMatrixPressed(){
    QString selected = QFileDialog::getOpenFileName(this, tr("Choose a data matrix file:"));
    ui->matrixPathLine->setText(selected);
}

/**
 * @brief MainMenu::cutoffAutoSet
 * manages text on cutoff textLine when auto mode is set
 * @param value true means auto is now checked, false that it is now unchecked
 */
void MainMenu::cutoffAutoSet(bool value){
    if(value){
        ui->CutoffLine->setText("Auto");
        ui->CutoffLine->setEnabled(false);
    }else{
        ui->CutoffLine->setText("");
        ui->CutoffLine->setEnabled(true);
    }
}

/**
 * @brief MainMenu::bootstrapEnabledSet
 * manages text on bootstrap replicates textLine
 * @param value true means bootstrap is now enabled, false otherwise
 */
void MainMenu::bootstrapEnabledSet(bool value){
    if(value){
        ui->replicatesLine->setText("");
        ui->replicatesLine->setEnabled(true);
    }else{
        ui->replicatesLine->setText("Disabled");
        ui->replicatesLine->setEnabled(false);
    }
}

/**
 * @brief MainMenu::log
 * write a string on log textArea
 * @param str text to be written
 */
void MainMenu::log(std::string str){
    ui->logText->moveCursor (QTextCursor::End);
    ui->logText->insertPlainText(str.data());
}

/**
 * @brief MainMenu::runBML
 * Launches execution of BML,
 * execution reports are written on log textArea.
 * Input informations are taken by the GUI.
 */
void MainMenu::runBML(){

    this->setEnabled(false);
    log("Validating input ...\n");
    if(inputValuesPreliminaryCheck()){
        log("Validation complete!\n");

        std::string job = ui->jobLine->text().toStdString();
        log("Prepared new job: ");
        log(job);
        log("\n");
        std::string inputFile = ui->matrixPathLine->text().toStdString();

        int numberOfThreads = ui->numberOfThreadsLine->text().toInt();
        int nTrees = ui->reseedsLine->text().toInt();
        double pthres = ui->thresholdLine->text().toFloat();
        bool performBootstrap = ui->checkBootstrap->isChecked();
        int nRep = ui->replicatesLine->text().toInt();

        bool autoCutoff = ui->autoCutoffCheckBox->isChecked();
        int cutoff = ui->CutoffLine->text().toInt();

        this->setEnabled(true);
        this->ui->runButton->setEnabled(false);
        this->ui->runButton->setText("Executing...");
        log("Execution in progress, plase wait (it might take a long time depending on your input) ...\n");
        this->repaint();
        qApp->processEvents();

        QFutureWatcher<int>* watcher = new QFutureWatcher<int>();
        connect(watcher, SIGNAL(finished()), this, SLOT(executeFinished()), Qt::QueuedConnection);

        struct Input* input = new struct Input();
        input->autoCutoff = autoCutoff;
        input->nRep = nRep;
        input->nTrees = nTrees;
        input->numberOfThreads = numberOfThreads;
        input->pthres = pthres;
        input->performBootstrap = performBootstrap;
        input->job = job;
        input->inputFile = inputFile;
        input->cutoff = cutoff;

        QFuture<int> future = QtConcurrent::run(execute,  input); // run BML on a separate thread

        watcher->setFuture(future);

    }else{
        log("Failed! Check input fields\n");
    }
    this->setEnabled(true);

}

/**
 * @brief MainMenu::executeFinished
 * manages the end of the execution in execute's thread
 */
void MainMenu::executeFinished(){

    QFutureWatcher<int>* watcher = static_cast< QFutureWatcher<int>* >(sender());

    if(watcher->result() == -1){
        log("Warning: An error occurred during the execution\n");
    } else {
        log("Task completed successfully!\n");
        log("if you want to see your output you can copy the content of the \n");
        log("output tree file (the one in .dot format) and paste it in the \n");
        log("help->Graphviz Online window. Note that a internet connection is needed.\n");
    }


    log("\njob completed!\n-----------------------------------------------------------------------------------\n\n");

    disconnect(watcher, SIGNAL(finished()), this, SLOT(executeFinished()));

    if(watcher->result() == 0){
         QDesktopServices::openUrl(QUrl::fromLocalFile(QFileInfo("output").absoluteFilePath() ));  // open output folder
    }

    delete(watcher);

    this->ui->runButton->setText("Execute BML analysis");
    this->ui->runButton->setEnabled(true);
}

/**
 * @brief wlog
 * prints str on the log area
 *
 * @param str    text to be written
 */
void MainMenu::wlog(string str){
    this->ui->logText->moveCursor (QTextCursor::End);
    ui->logText->insertPlainText(str.data());
}

/**
 * @brief MainMenu::file_exist
 * checks if a file exists
 * @param fileName file to be checked
 *
 * @return true if it exists, false otherwise
 */
bool MainMenu::file_exist(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

/**
 * @brief MainMenu::inputValuesPreliminaryCheck
 * Checks if the input values preset at the moment in the GUI are valid
 * for an execution of BML. If there are problems, a message is written on the
 * log textArea
 *
 * @return true if check is successfull, false otherwise
 */
bool MainMenu::inputValuesPreliminaryCheck(){

    if(ui->jobLine->text().length() <= 0){                    // jobID
        log("Error: no job ID selected\n");
        return false;
    }

    if(ui->matrixPathLine->text().length() <= 0){             // inputMatrix
        log("Error: no input matrix selected\n");
        return false;
    }

    if(!file_exist(ui->matrixPathLine->text().toStdString())){
        log("Error: invalid input matrix selected, check your path\n");
        return false;
    }

    if(ui->numberOfThreadsLine->text().length() <= 0){       // threads
        log("Error: number of threads was not selected\n");
        return false;
    }
    int nt = ui->numberOfThreadsLine->text().toInt();
    if(nt < 1){
        log("Error: number of threads must be an integer grater than 0\n");
        return false;
    }

    if(!ui->autoCutoffCheckBox->isChecked()){
        if(ui->CutoffLine->text().length() <= 0){            // cutoff
            log("Error: cutoff not set\n");
            return false;
        }
        int nt = ui->CutoffLine->text().toInt();
        if(nt <= 0){
            log("Error: cutoff must be an integer grater than 0\n");
            return false;
        }
    }


    if(ui->reseedsLine->text().length() <= 0){               // number of tree reseeds
        log("Error: number of tree reseeds was not selected\n");
        return false;
    }
    nt = ui->reseedsLine->text().toInt();
    if(nt < 1){
        log("Error: number of tree reseeds must be an integer grater than 0\n");
        return false;
    }

    if(ui->thresholdLine->text().length() <= 0){              // threshold
        log("Error: threshold was not selected\n");
        return false;
    }
    float ft = ui->thresholdLine->text().toFloat();
    if(ft <= 0 || ft >= 1){
        log("Error: threshold must be a float grater than 0 and lower than 1\n");
        return false;
    }

    if(ui->checkBootstrap->isChecked()){

        if(ui->replicatesLine->text().length() <= 0){         // number of bootstrap replicates
            log("Error: number of bootstrap replicates not set\n");
            return false;
        }
        int nt = ui->replicatesLine->text().toInt();
        if(nt < 0){
            log("Error: number of bootstrap replicates must be a positive integer\n");
            return false;
        }

    }

    return true;
}
