#include "mainwindow.h"
#include "cellUI.h"
#include "geneListUI.h"
#include "geneUI.h"
#include "summaryUI.h"
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

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
                                                    tr("Select population summary"),
                                                    "/home/Documents",
                                                    tr("Summary (*.dat)"));
    ui->summPath->setText(filename);
}

void MainWindow::on_geneListButton_clicked()
{
    QString filename =QFileDialog::getOpenFileName(this,
                                                   tr("Select gene list"),
                                                   "/home/Documents",
                                                   tr("Gene list (*.dat)"));
    ui->genelistPath->setText(filename);
}
