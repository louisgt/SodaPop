#ifndef GENEUI_H
#define GENEUI_H

#include <QDialog>
#include <QFile>
#include <QFileDialog>
#include <QTextStream>
#include <QMessageBox>
//#include <seqan/sequence.h>
//#include <seqan/seq_io.h>

namespace Ui {
class geneUI;
}

class geneUI : public QDialog
{
    Q_OBJECT

public:
    explicit geneUI(QWidget *parent = nullptr);
    ~geneUI();

private slots:
    void on_buttonBox_accepted();

    void on_buttonBox_rejected();

    void on_chooseDir_clicked();

    void on_pushButton_clicked();

private:
    Ui::geneUI *ui;
};

#endif // geneUI_H
