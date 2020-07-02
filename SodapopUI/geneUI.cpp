#include "geneUI.h"
#include "ui_geneUI.h"

geneUI::geneUI(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::geneUI)
{
    ui->setupUi(this);

    QRegularExpression regExp("\\d+");
    QRegularExpression regExpFloat("\\d+.?\\d+");

    ui->numID->setValidator(new QRegularExpressionValidator(regExp));
    ui->abundance->setValidator(new QRegularExpressionValidator(regExp));
    ui->stability->setValidator(new QRegularExpressionValidator(regExpFloat));

}

geneUI::~geneUI()
{
    delete ui;
}

void geneUI::on_buttonBox_accepted()
{


    QString NumID = ui->numID->text();
    QString GeneID = ui->geneID->text();
    QString aaSeq = ui->aaSeq->text();
    QString nucSeq = ui->nucSeq->text();
    QString stability = ui->stability->text();
    QString concentration = ui->abundance->text();
    QString essential = ui->essential->isChecked() ? "1":"0";



    //verification of inputs

    if(ui->dir->text().isEmpty()){
        QMessageBox error;
        error.setText("You must choose a directory");
        error.exec();
        return;
    }
    else{


        QFile file(ui->dir->text()+NumID+".gene");
        if (file.exists()){
            QMessageBox overwrite;
            //no need to store pointer to yes button
            overwrite.addButton(tr("Yes"), QMessageBox::ActionRole);
            QAbstractButton *no = overwrite.addButton(tr("No"), QMessageBox::ActionRole);
            overwrite.setText("This file already exists, would you like to overwrite it?");
            overwrite.exec();
            if(overwrite.clickedButton() == no || overwrite.clickedButton() == overwrite.escapeButton()){
                return;
            }

        }
        if(!file.open(QIODevice::WriteOnly | QIODevice::Text)){
            return;
        }
        QTextStream out(&file);
        out<<"Gene_NUM\t"<<NumID+"\n"<<"Gene_ID\t"<<GeneID+"\n"<<"A_Seq\t"<<aaSeq+"\n"<<"N_Seq\t"<<nucSeq+"\n"<<"E\t"<<essential+"\n"<<"DG\t"<<stability+"\n"<<"CONC\t"<<concentration;
    }

}


void geneUI::on_buttonBox_rejected()
{
    this->~geneUI();
}

void geneUI::on_chooseDir_clicked()
{
    QString directory = QFileDialog::getExistingDirectory(this,
                                                         tr("Select the target directory"),
                                                         "/home/Documents",
                                                         QFileDialog::ShowDirsOnly);
    ui->dir->setText(directory+"/");
}

void geneUI::on_pushButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                                     tr("Select a FASTA file"),
                                                     "/home/Documents");
    ui->lineEdit->setText(fileName);

}
