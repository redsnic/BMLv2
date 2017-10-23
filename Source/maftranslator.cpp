#include "maftranslator.h"
#include "ui_maftranslator.h"

MAFtranslator::MAFtranslator(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MAFtranslator)
{
    ui->setupUi(this);
    this->setAttribute( Qt::WA_DeleteOnClose );
}

MAFtranslator::~MAFtranslator()
{
    delete ui;
}

/**
 * @brief MAFtranslator::closeEvent
 * standard procedure to alert the parent that a children window is closed,
 * in this case it is used to re-enable the button to open the MAF preprocessor.
 * @param event
 */
void MAFtranslator::closeEvent(QCloseEvent *event){
    emit closeEvent();
}

/**
 * @brief MAFtranslator::getInput
 * MAF file user request
 */
void MAFtranslator::getInput(){
    QString selected = QFileDialog::getOpenFileName(this, tr("Choose a MAF file:"));
    ui->lineMAF->setText(selected);
}

/**
 * @brief MAFtranslator::setOutput
 * Output location user request
 */
void MAFtranslator::setOutput(){
    QString selected = QFileDialog::getSaveFileName(this, tr("Choose output location for data matrix:"));
    ui->lineMatrix->setText(selected);
}

/**
 * @brief MAFtranslator::translate
 * uses a MAFtoMatrixTranslator to elaborate a MAF file to a BML data matrix.
 * It also gets and prepares all the needed informations from the GUI.
 */
void MAFtranslator::translate(){
    MAFtoMatrixTranslator* maf = new MAFtoMatrixTranslator;
    this->setEnabled(false);

    std::string inPath = ui->lineMAF->text().toStdString();      //I/O locations
    std::string outPath = ui->lineMatrix->text().toStdString();

    bool all = ui->allCheckBox->isChecked();                     // parameters
    bool standard = ui->standardCheckBox->isChecked();

    bool allGenes = ui->CheckAllGenes->isChecked();
    bool useCustomGenes = !allGenes;

    std::vector<std::string> genes(0);

    std::vector<std::string> mutations(0);

    if(inPath.length() <= 0 || outPath.length() <= 0){           // no input or output location choosen
        QMessageBox fail;
        fail.setIcon(QMessageBox::Critical);
        fail.setText("Please insert input and output informations.");
        fail.exec();
        this->setEnabled(true);
        return;
    }

    if(!all && !standard){                                       // prepare selected mutations' list
        QString mutationsText = ui->mutationText->toPlainText();
        QStringList list = mutationsText.split(QRegExp("\\s+"));
        for(int i = 0; i<list.size(); i++){
            mutations.push_back(list.at(i).toUpper().toStdString());
        }
    }

    if(useCustomGenes){                                          // prepare selected genes' list
        QString geneText = ui->textCustomGenes->toPlainText();
        QStringList list = geneText.split(QRegExp("\\s+"));
        for(int i = 0; i<list.size(); i++){
            if(list.at(i) != ""){
                genes.push_back(list.at(i).toUpper().toStdString());
            }
        }
        std::sort(genes.begin(), genes.end(), [](std::string a, std::string b){ return b < a; }); // sort in lex order
        if(genes.size() == 0){
            QMessageBox fail;
            fail.setIcon(QMessageBox::Critical);
            fail.setText("No gene set to be used in custom genes mode, enter the genes'names in the text area.");
            fail.exec();
            this->setEnabled(true);
            return;
        }
    }

    try{                                                         // tries to execute the translation, intercepts I/O errors
        int val = maf->read(inPath, all, standard, mutations, useCustomGenes, genes);
        val -= maf->save(outPath);
        if (val<0){
            QMessageBox fail;
            fail.setIcon(QMessageBox::Critical);
            fail.setText("Error in input or output paths, please check.");
            fail.exec();
            this->setEnabled(true);
            return;
        }
        QMessageBox done;
        done.setIcon(QMessageBox::Information);
        done.setText("Data matrix created successfully.");
        done.exec();
    }catch(...){
        QMessageBox fail;
        fail.setIcon(QMessageBox::Critical);
        fail.setText("Error reading MAF file, check your path.");
        fail.exec();
    }

    this->setEnabled(true);
    delete(maf);
}

/**
 * @brief MAFtranslator::setStandard
 * standard checkbox management, so that it can't be selected with all chckbox.
 * it also disables custom mutations set textArea
 * @param value chekbox value
 */
void MAFtranslator::setStandard(bool value){
    if(value){
        ui->mutationText->setEnabled(false);
        ui->allCheckBox->setChecked(false);
    }else{
        ui->mutationText->setEnabled(true);
    }
}

/**
 * @brief MAFtranslator::setAll
 * all checkbox management, so that it can't be selected with standard chckbox.
 * it also disables custom mutations set textArea
 * @param value chekbox value
 */
void MAFtranslator::setAll(bool value){
    if(value){
        ui->mutationText->setEnabled(false);
        ui->standardCheckBox->setChecked(false);
    }else{
        ui->mutationText->setEnabled(true);
    }
}


/**
 * @brief MAFtranslator::setAllGenes
 * Reflects the choice of operating over all genes
 * by unchecking some checkbox and disabling some text fields
 */
void MAFtranslator::setAllGenes(){
    if(!ui->CheckAllGenes->isChecked()){
        ui->textCustomGenes->setEnabled(true);
    }else{
        ui->textCustomGenes->setEnabled(false);
    }
}




