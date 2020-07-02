#include "geneListUI.h"
#include "ui_geneListUI.h"

geneListUI::geneListUI(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::geneListUI)
{
    ui->setupUi(this);
}

geneListUI::~geneListUI()
{
    delete ui;
}

QStringList geneFiles;
QStringList geneIDs;

QString geneListUI::parseGeneFile(QString fileName){
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        QMessageBox error;
        error.setText("Could not open gene file: "+fileName);
        error.exec();
    }
    while (!file.atEnd()) {
        QString line = QString(file.readLine());
        if(line.contains("Gene_ID")){
            return line.section('\t', -1);
        }
    }
    QMessageBox error;
    error.setText("Invalid gene file: " + fileName);
    error.exec();

}

void geneListUI::on_buttonBox_accepted()
{
    QFile file(ui->genePath->text()+ui->listName->text()+".dat");
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        return;
    }
    QTextStream out(&file);
    out<<"Gene_Count\t" + QString::number(geneFiles.size())+"\n";
    for(int i=0; i<geneFiles.size() ; i++){
        QString gene = geneFiles[i];
        out<<"G\t"+gene.section('/' , -1)+"\t"+geneIDs[i];
    }

    this->~geneListUI();
}

void geneListUI::on_buttonBox_rejected()
{
    this->~geneListUI();
}

void geneListUI::on_dirButton_clicked()
{
    geneFiles = QFileDialog::getOpenFileNames(this,
                                              tr("Select the genes you would like to include"),
                                              "/home/Documents",
                                              "Gene files (*.gene)");
    //flushing any previous content
    geneIDs = QStringList();

    for(const auto& fileName : geneFiles){
        geneIDs.append(parseGeneFile(fileName));
    }
    ui->genePath->setText(geneFiles[0].section('/', 0, -2)+"/");
    ui->listWidget->addItems(geneIDs);
}
