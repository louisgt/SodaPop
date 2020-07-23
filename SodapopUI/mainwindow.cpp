#include "mainwindow.h"
#include "cellUI.h"
#include "geneListUI.h"
#include "geneUI.h"
#include "summaryUI.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QRegularExpression num("\\d+");

    ui->sizeEdit->setValidator(new QRegularExpressionValidator(num));
    ui->genEdit->setValidator(new QRegularExpressionValidator(num));
    ui->snapEdit->setValidator(new QRegularExpressionValidator(num));

    funcType.insert("Folding", "1");
    funcType.insert("Misfolding", "2");
    funcType.insert("Additive", "3");
    funcType.insert("Multiplicative", "4");
    funcType.insert("Neutral", "5");
    funcType.insert("Growth rate", "6");

}

MainWindow::~MainWindow()
{
    delete ui;
}
QMap<QString, QString> funcType;


const QStringList biologicalFit = { "Folding", "Misfolding", "Growth rate" };
const QStringList theoreticalFit = { "Additive", "Multiplicative", "Neutral" };





void MainWindow::on_actionSummary_triggered()
{
    summaryUI *summWindow = new summaryUI();
    summWindow->show();
}

void MainWindow::on_actionCell_triggered()
{
    cellUI *cellWindow = new cellUI();
    cellWindow->show();
}

void MainWindow::on_actionGene_list_triggered()
{
    geneListUI *genelistWindow = new geneListUI();
    genelistWindow->show();
}

void MainWindow::on_actionGene_triggered()
{
    geneUI *geneWindow = new geneUI();
    geneWindow->show();
}

void MainWindow::on_actionSnapshot_triggered()
{
    /*runs sodasumm*/
    qDebug() << "yo wassup";
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Select population summary"),
                                                    "/home",
                                                    tr("Population summary files (*.dat)")
                                                    );
    QString program = "/media/tacito/Data/SerohijosLab/SodaPop/SodaPop-master/sodasumm";
    process.setWorkingDirectory("/home/tacito/Documents/gene");
    QStringList arg;
    QStringList error;

    arg << filename << "0";
    process.start(program, arg);
    process.setReadChannel(process.StandardError);
     do{
        qDebug() << "Entered the loop";
        process.waitForReadyRead();
        qDebug() << process.readLine();
        qDebug() << "Left the loop";
    }while (!process.waitForFinished());
    qDebug()<< process.readAllStandardError();
    qDebug()<<"done";


    /*QMessageBox wip;
    wip.setText("Work in progress");
    wip.exec();*/

}



void MainWindow::on_theoretical_clicked(bool checked)
{

    ui->fitLandscape->setEnabled(false);
    ui->DFE->setEnabled(true);

    ui->functionType->clear();
    ui->functionType->addItems(theoreticalFit);
}



void MainWindow::on_experimental_clicked(bool checked)
{
    ui->fitLandscape->setEnabled(true);
    ui->DFE->setEnabled(false);

    ui->functionType->clear();
    ui->functionType->addItems(biologicalFit);
}

void MainWindow::on_summButton_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Select population snapshot"),
                                                    "/home/",
                                                    tr("Snapshot (*.snap)"));
    ui->summPath->setText(filename);
}

void MainWindow::on_geneListButton_clicked()
{
    QString filename =QFileDialog::getOpenFileName(this,
                                                   tr("Select gene list"),
                                                   "/home/",
                                                   tr("Gene list (*.dat)"));
    ui->genelistPath->setText(filename);
}

void MainWindow::on_pushButton_clicked()
{

}

void MainWindow::on_landscapeA_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Select matrix"),
                                                    "/home/",
                                                    tr("Substitution matrix (*.matrix)"));
    ui->landscapeAEdit->setText(filename);
}

void MainWindow::on_landscapeB_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Select matrix"),
                                                    "/home/",
                                                    tr("Substitution matrix (*.matrix)"));
    ui->landscapeBEdit->setText(filename);
}
QStringList MainWindow::buildArgs()

{


    QStringList arguments;

    arguments << getSimType() << "-f" << funcType[ui->functionType->currentText()] << "-p" << ui->summPath->text() <<
                 "-g" << ui->genelistPath->text() << "-o test" << "-n " << ui->sizeEdit->text() <<
                 "-m" << ui->genEdit->text() << "-t" << ui->snapEdit->text();


    return arguments;
}

QStringList MainWindow::getSimType()
{
    QStringList simArgs;
    QString choice = ui->functionTypeGroup->checkedButton()->text();
    if(choice=="Theoretical")
    {
        simArgs << "--sim-type s";
        if(ui->distributionType->currentText()=="Normal"){
            simArgs << "--normal" << "--alpha" << ui->shapeEdit->text() << "--beta" << ui->scaleEdit->text();
        }
        else if(ui->distributionType->currentText() == "Gamma"){
            simArgs << "--gamma" << "--alpha" << ui->shapeEdit->text() << "--beta" << ui->scaleEdit->text();
        }

    }
    else if(choice=="Experimental")
    {
        simArgs<< "--sim-type stability" << "-i" << ui->landscapeAEdit->text();

    }
    return simArgs;
}

void MainWindow::on_pushButton_released(){
    process.setWorkingDirectory(ui->workingDir->text());
    QString program("./sodapop");

    qDebug() << buildArgs();

    process.start(program, buildArgs());
    process.waitForFinished();

    qDebug()<< process.readAllStandardError();
    return;

}

void MainWindow::on_chooseWorkDir_clicked()
{
    QString dir = QFileDialog::getExistingDirectory(this,
                                                    "Select working directory");
    ui->workingDir->setText(dir);
    process.setWorkingDirectory(dir);
}
