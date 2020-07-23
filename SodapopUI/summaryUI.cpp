#include "summaryUI.h"
#include "ui_summaryUI.h"

summaryUI::summaryUI(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::summaryUI)
{
    ui->setupUi(this);
    QHeaderView *header = ui->tableWidget->horizontalHeader();
    header->setStretchLastSection(true);
    countRe.setPattern("\\d+");
}

summaryUI::~summaryUI()
{
    delete ui;
}
QStringList cells;

void summaryUI::on_cellBrowse_clicked()
{

    ui->tableWidget->clearContents();
    ui->cellsDir->clear();
    QStringList cellIds;
    cells = QFileDialog::getOpenFileNames(this,
                                          tr("Select the cells you would like to include in your summary"),
                                          "/home/",
                                          "Cell files (*.cell)");
    if(cells.isEmpty()){
        return;
    }
    for(const auto& filename:cells){
        QString id = filename.section('/', -1).section('.', 0, 1);
        cellIds.append(id);
    }
    ui->cellsDir->setText(cells[0].section('/', 0, -2)+'/');
    ui->tableWidget->setRowCount(cells.size());
    for (int i=0; i<cellIds.size(); i++) {
        QLineEdit *countEdit = new QLineEdit();
        countEdit->setValidator(new QRegExpValidator(countRe));

        ui->tableWidget->setCellWidget(i, 0, new QCheckBox(cellIds[i], this));
        ui->tableWidget->setCellWidget(i, 1, countEdit);
        ui->tableWidget->setCellWidget(i, 2, new QLineEdit());

    }



}

void summaryUI::on_buttonBox_accepted()
{
    //QRegularExpressionValidator *countVal = new QRegularExpressionValidator(QRegularExpression("//d+"), this);
    //Validator does NOT like temporaries
    //int idx = 0;

    QFile file(ui->cellsDir->text()+ui->summaryName->text()+".dat");
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        QMessageBox error;
        error.setText("Could not create file at the requested location");
        error.exec();
        return;
    }
    QTextStream out(&file);
    out<<"//Population summary\n//Count\tFile\tComment\n";
    for (int i=0; i<ui->tableWidget->rowCount(); i++) {
        QCheckBox *cell = (QCheckBox*)ui->tableWidget->cellWidget(i, 0);
        if(!cell->isChecked()){
            continue;
        }

        //typecasting is the devil, but it works
        QLineEdit *countEdit = (QLineEdit*)ui->tableWidget->cellWidget(i, 1);
        QLineEdit *commentEdit = (QLineEdit*)ui->tableWidget->cellWidget(i, 2);

        QString count = countEdit->text();
        QString comment = commentEdit->text();

#
        /*if(!countVal->validate(count, idx))
        {
            file.remove();
            QMessageBox error;
            error.setText("Count must be an integer");
            error.exec();
            return;
        }*/

        out<<count+'\t'+cells[i]+'\t'+comment+'\n';
    }
    this->~summaryUI();
}

void summaryUI::on_buttonBox_rejected()
{
    this->~summaryUI();
}
