/********************************************************************************
** Form generated from reading UI file 'summaryUI.ui'
**
** Created by: Qt User Interface Compiler version 5.9.9
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SUMMARYUI_H
#define UI_SUMMARYUI_H

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

class Ui_summaryUI
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QPushButton *cellBrowse;
    QLineEdit *cellsDir;
    QWidget *horizontalLayoutWidget_2;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label;
    QLineEdit *summaryName;
    QDialogButtonBox *buttonBox;
    QWidget *gridLayoutWidget;
    QGridLayout *gridLayout;

    void setupUi(QDialog *summaryUI)
    {
        if (summaryUI->objectName().isEmpty())
            summaryUI->setObjectName(QStringLiteral("summaryUI"));
        summaryUI->resize(700, 400);
        horizontalLayoutWidget = new QWidget(summaryUI);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(10, 10, 681, 31));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        cellBrowse = new QPushButton(horizontalLayoutWidget);
        cellBrowse->setObjectName(QStringLiteral("cellBrowse"));

        horizontalLayout->addWidget(cellBrowse);

        cellsDir = new QLineEdit(horizontalLayoutWidget);
        cellsDir->setObjectName(QStringLiteral("cellsDir"));

        horizontalLayout->addWidget(cellsDir);

        horizontalLayoutWidget_2 = new QWidget(summaryUI);
        horizontalLayoutWidget_2->setObjectName(QStringLiteral("horizontalLayoutWidget_2"));
        horizontalLayoutWidget_2->setGeometry(QRect(20, 370, 661, 27));
        horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget_2);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(horizontalLayoutWidget_2);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout_2->addWidget(label);

        summaryName = new QLineEdit(horizontalLayoutWidget_2);
        summaryName->setObjectName(QStringLiteral("summaryName"));

        horizontalLayout_2->addWidget(summaryName);

        buttonBox = new QDialogButtonBox(horizontalLayoutWidget_2);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        horizontalLayout_2->addWidget(buttonBox);

        gridLayoutWidget = new QWidget(summaryUI);
        gridLayoutWidget->setObjectName(QStringLiteral("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(9, 49, 681, 311));
        gridLayout = new QGridLayout(gridLayoutWidget);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);

        retranslateUi(summaryUI);

        QMetaObject::connectSlotsByName(summaryUI);
    } // setupUi

    void retranslateUi(QDialog *summaryUI)
    {
        summaryUI->setWindowTitle(QApplication::translate("summaryUI", "Create a population summary", Q_NULLPTR));
        cellBrowse->setText(QApplication::translate("summaryUI", "Choose directory", Q_NULLPTR));
        label->setText(QApplication::translate("summaryUI", "Name your summary", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class summaryUI: public Ui_summaryUI {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SUMMARYUI_H
