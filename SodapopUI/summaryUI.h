#ifndef SUMMARYUI_H
#define SUMMARYUI_H

#include <QDialog>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QCheckBox>


namespace Ui {
class summaryUI;
}

class summaryUI : public QDialog
{
    Q_OBJECT

public:
    explicit summaryUI(QWidget *parent = nullptr);
    ~summaryUI();

private slots:
    void on_cellBrowse_clicked();

    void on_buttonBox_accepted();

    void on_buttonBox_rejected();

private:
    Ui::summaryUI *ui;
    QRegExp countRe;
};

#endif // summaryUI_H
