#ifndef GENELISTUI_H
#define GENELISTUI_H

#include <QDialog>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>

namespace Ui {
class geneListUI;
}

class geneListUI : public QDialog
{
    Q_OBJECT

public:
    explicit geneListUI(QWidget *parent = nullptr);
    QString parseGeneFile(QString fileName);

    ~geneListUI();


private slots:
    void on_buttonBox_accepted();

    void on_buttonBox_rejected();

    void on_dirButton_clicked();

private:
    Ui::geneListUI *ui;

};

#endif // GENELISTUI_H
