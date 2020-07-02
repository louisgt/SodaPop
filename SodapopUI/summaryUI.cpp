#include "summaryUI.h"
#include "ui_summaryUI.h"

summaryUI::summaryUI(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::summaryUI)
{
    ui->setupUi(this);
}

summaryUI::~summaryUI()
{
    delete ui;
}
QStringList cells;

void summaryUI::on_cellBrowse_clicked()
{
//    int column = 0;
    cells = QFileDialog::getOpenFileNames(this,
                                          tr("Select the cells you'd like to include in your summary"),
                                          "/home/Documents/",
                                          "Cell files (*.cell)");
    ui->cellsDir->setText(cells[1].section('/',0,-1)+'/');
//    for(const auto& cell : cells){
//       QString cellName = cell.split('/')[-1];
//       ui->gridLayout->addWidget(new QLabel(cellName), 0, column);
//       ui->gridLayout->addWidget(new)
//    }


}

void summaryUI::on_buttonBox_accepted()
{
    QFile file(ui->cellsDir->text()+ui->summaryName->text()+".dat");
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        QMessageBox error;
        error.setText("Could not create file at the requested location");
        error.exec();
    }
    QTextStream out(&file);
    out<<"//Population summary\n//Count\tFile\tComment\n";

}

void summaryUI::on_buttonBox_rejected()
{
    this->~summaryUI();
}
