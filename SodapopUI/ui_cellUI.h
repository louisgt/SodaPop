/********************************************************************************
** Form generated from reading UI file 'cellUI.ui'
**
** Created by: Qt User Interface Compiler version 5.9.9
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CELLUI_H
#define UI_CELLUI_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_cellUI
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout_3;
    QPushButton *geneListButton;
    QLineEdit *geneListPath;
    QWidget *gridLayoutWidget;
    QGridLayout *gridLayout;
    QWidget *horizontalLayoutWidget_2;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QLineEdit *cellID;
    QDialogButtonBox *buttonBox;
    QWidget *horizontalLayoutWidget_3;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QLineEdit *mutationRate;

    void setupUi(QDialog *cellUI)
    {
        if (cellUI->objectName().isEmpty())
            cellUI->setObjectName(QStringLiteral("cellUI"));
        cellUI->resize(600, 387);
        horizontalLayoutWidget = new QWidget(cellUI);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(9, 9, 582, 27));
        horizontalLayout_3 = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        geneListButton = new QPushButton(horizontalLayoutWidget);
        geneListButton->setObjectName(QStringLiteral("geneListButton"));

        horizontalLayout_3->addWidget(geneListButton);

        geneListPath = new QLineEdit(horizontalLayoutWidget);
        geneListPath->setObjectName(QStringLiteral("geneListPath"));

        horizontalLayout_3->addWidget(geneListPath);

        gridLayoutWidget = new QWidget(cellUI);
        gridLayoutWidget->setObjectName(QStringLiteral("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(9, 39, 581, 281));
        gridLayout = new QGridLayout(gridLayoutWidget);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayoutWidget_2 = new QWidget(cellUI);
        horizontalLayoutWidget_2->setObjectName(QStringLiteral("horizontalLayoutWidget_2"));
        horizontalLayoutWidget_2->setGeometry(QRect(10, 350, 581, 31));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget_2);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(horizontalLayoutWidget_2);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout->addWidget(label);

        cellID = new QLineEdit(horizontalLayoutWidget_2);
        cellID->setObjectName(QStringLiteral("cellID"));

        horizontalLayout->addWidget(cellID);

        buttonBox = new QDialogButtonBox(horizontalLayoutWidget_2);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        horizontalLayout->addWidget(buttonBox);

        horizontalLayoutWidget_3 = new QWidget(cellUI);
        horizontalLayoutWidget_3->setObjectName(QStringLiteral("horizontalLayoutWidget_3"));
        horizontalLayoutWidget_3->setGeometry(QRect(10, 320, 581, 31));
        horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget_3);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(horizontalLayoutWidget_3);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout_2->addWidget(label_2);

        mutationRate = new QLineEdit(horizontalLayoutWidget_3);
        mutationRate->setObjectName(QStringLiteral("mutationRate"));

        horizontalLayout_2->addWidget(mutationRate);


        retranslateUi(cellUI);

        QMetaObject::connectSlotsByName(cellUI);
    } // setupUi

    void retranslateUi(QDialog *cellUI)
    {
        cellUI->setWindowTitle(QApplication::translate("cellUI", "Select the genes to include in the cell's genome", Q_NULLPTR));
        geneListButton->setText(QApplication::translate("cellUI", "Select list of genes", Q_NULLPTR));
        geneListPath->setText(QApplication::translate("cellUI", "Gene list path", Q_NULLPTR));
        label->setText(QApplication::translate("cellUI", "Cell ID", Q_NULLPTR));
        label_2->setText(QApplication::translate("cellUI", "Mutation rate", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class cellUI: public Ui_cellUI {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CELLUI_H
