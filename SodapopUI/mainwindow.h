#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>
#include "ui_mainwindow.h"
#include <QFileDialog>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    int foo;

private slots:
    void on_actionSummary_triggered();

    void on_actionCell_triggered();

    void on_actionGene_list_triggered();

    void on_actionGene_triggered();

    void on_actionSnapshot_triggered();

    void on_theoretical_clicked(bool checked);

    void on_experimental_clicked(bool checked);

    void on_summButton_clicked();

    void on_geneListButton_clicked();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
