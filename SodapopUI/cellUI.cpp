#include "cellUI.h"
#include "ui_cellUI.h"

cellUI::cellUI(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::cellUI)
{
    ui->setupUi(this);
    QHeaderView *header = ui->tableWidget->horizontalHeader();
    header->setStretchLastSection(true);
}

cellUI::~cellUI()
{
    delete ui;
}
QMap<QString, QString> geneMap;


void cellUI::on_buttonBox_accepted()
{
    if(ui->cellID->text().isEmpty()){
        QMessageBox error;
        error.setText("Please name your cell");
        error.exec();
        return;
    }

    QString filePath = ui->geneListPath->text().section('/', 0, -2);
    QFile file(filePath+"/"+ui->cellID->text()+".cell");

    if(!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        QMessageBox error;
        error.setText("Could not create file at specified location");
        error.exec();
        return;
    }
    QTextStream out(&file);
    out<<"org_id\t"+ui->cellID->text()+'\n'+"mrate\t"+ui->mutationRate->text()+'\n';

    //There are better ways to do this.
    //This only works because the only instances are those of interest


    for (auto button : findChildren<QCheckBox*>())if(button->isChecked()) {
        QString key = button->text();
        out<<"G\t"+geneMap[key].split('.')[0]+'\n';
    }
    this->~cellUI();
}

void cellUI::on_buttonBox_rejected()
{
    this->~cellUI();
}

void cellUI::on_geneListButton_clicked()
{
    ui->tableWidget->clearContents();
    QString filePath = QFileDialog::getOpenFileName(this,
                                                    tr("Select gene list file"),
                                                    "/home/",
                                                    "Gene list files (*.dat)");
    ui->geneListPath->setText(filePath);
    QFile file(filePath);

    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        QMessageBox error;
        error.setText("Could not open gene list file");
        error.exec();
        return;
    }
    int count=0;
    while (!file.atEnd()) {
        QString line = QString(file.readLine());
        if(line.contains("G\t")){
            ui->tableWidget->insertRow(count);
            ui->tableWidget->setCellWidget(count, 0 , new QCheckBox(line.section('\t', -1), this));
            geneMap.insert(line.section('\t', -1), line.section('\t',-2));
            count++;
        }
    }

}
