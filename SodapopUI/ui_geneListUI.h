/********************************************************************************
** Form generated from reading UI file 'geneListUI.ui'
**
** Created by: Qt User Interface Compiler version 5.9.9
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GENELISTUI_H
#define UI_GENELISTUI_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_geneListUI
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout_4;
    QPushButton *dirButton;
    QLineEdit *genePath;
    QWidget *gridLayoutWidget;
    QGridLayout *gridLayout_2;
    QListWidget *listWidget;
    QWidget *horizontalLayoutWidget_2;
    QHBoxLayout *horizontalLayout;
    QLineEdit *listName;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *geneListUI)
    {
        if (geneListUI->objectName().isEmpty())
            geneListUI->setObjectName(QStringLiteral("geneListUI"));
        geneListUI->resize(600, 354);
        horizontalLayoutWidget = new QWidget(geneListUI);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(10, 10, 582, 27));
        horizontalLayout_4 = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, 0, 0);
        dirButton = new QPushButton(horizontalLayoutWidget);
        dirButton->setObjectName(QStringLiteral("dirButton"));

        horizontalLayout_4->addWidget(dirButton);

        genePath = new QLineEdit(horizontalLayoutWidget);
        genePath->setObjectName(QStringLiteral("genePath"));

        horizontalLayout_4->addWidget(genePath);

        gridLayoutWidget = new QWidget(geneListUI);
        gridLayoutWidget->setObjectName(QStringLiteral("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(10, 40, 581, 271));
        gridLayout_2 = new QGridLayout(gridLayoutWidget);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        gridLayout_2->setContentsMargins(0, 0, 0, 0);
        listWidget = new QListWidget(gridLayoutWidget);
        listWidget->setObjectName(QStringLiteral("listWidget"));

        gridLayout_2->addWidget(listWidget, 0, 0, 1, 1);

        horizontalLayoutWidget_2 = new QWidget(geneListUI);
        horizontalLayoutWidget_2->setObjectName(QStringLiteral("horizontalLayoutWidget_2"));
        horizontalLayoutWidget_2->setGeometry(QRect(10, 320, 581, 31));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget_2);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        listName = new QLineEdit(horizontalLayoutWidget_2);
        listName->setObjectName(QStringLiteral("listName"));

        horizontalLayout->addWidget(listName);

        buttonBox = new QDialogButtonBox(horizontalLayoutWidget_2);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        horizontalLayout->addWidget(buttonBox);


        retranslateUi(geneListUI);
        QObject::connect(buttonBox, SIGNAL(accepted()), geneListUI, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), geneListUI, SLOT(reject()));

        QMetaObject::connectSlotsByName(geneListUI);
    } // setupUi

    void retranslateUi(QDialog *geneListUI)
    {
        geneListUI->setWindowTitle(QApplication::translate("geneListUI", "Select genes for gene list", Q_NULLPTR));
        dirButton->setText(QApplication::translate("geneListUI", "Select gene directory", Q_NULLPTR));
        genePath->setText(QApplication::translate("geneListUI", "Gene path", Q_NULLPTR));
        listName->setText(QApplication::translate("geneListUI", "Name of gene list", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class geneListUI: public Ui_geneListUI {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GENELISTUI_H
