#ifndef CELLUI_H
#define CELLUI_H

#include <QDialog>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QMap>
#include <QCheckBox>


namespace Ui {
class cellUI;
}

class cellUI : public QDialog
{
    Q_OBJECT

public:
    explicit cellUI(QWidget *parent = nullptr);
    ~cellUI();

private slots:

    void on_buttonBox_accepted();

    void on_buttonBox_rejected();

    void on_geneListButton_clicked();

private:
    Ui::cellUI *ui;
};

#endif // CELLUI_H
