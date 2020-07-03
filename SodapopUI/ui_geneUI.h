/********************************************************************************
** Form generated from reading UI file 'geneUI.ui'
**
** Created by: Qt User Interface Compiler version 5.9.9
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GENEUI_H
#define UI_GENEUI_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
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

class Ui_geneUI
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QPushButton *chooseDir;
    QLineEdit *dir;
    QDialogButtonBox *buttonBox;
    QWidget *horizontalLayoutWidget_2;
    QHBoxLayout *horizontalLayout_3;
    QPushButton *pushButton;
    QLineEdit *lineEdit;
    QWidget *gridLayoutWidget;
    QGridLayout *gridLayout;
    QLabel *label_5;
    QCheckBox *essential;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_4;
    QLabel *label_3;
    QLabel *label_6;
    QLineEdit *numID;
    QLineEdit *geneID;
    QLineEdit *aaSeq;
    QLineEdit *nucSeq;
    QLineEdit *stability;
    QLineEdit *abundance;

    void setupUi(QDialog *geneUI)
    {
        if (geneUI->objectName().isEmpty())
            geneUI->setObjectName(QStringLiteral("geneUI"));
        geneUI->resize(600, 436);
        horizontalLayoutWidget = new QWidget(geneUI);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(10, 400, 581, 31));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        chooseDir = new QPushButton(horizontalLayoutWidget);
        chooseDir->setObjectName(QStringLiteral("chooseDir"));

        horizontalLayout->addWidget(chooseDir);

        dir = new QLineEdit(horizontalLayoutWidget);
        dir->setObjectName(QStringLiteral("dir"));
        dir->setEnabled(true);

        horizontalLayout->addWidget(dir);

        buttonBox = new QDialogButtonBox(horizontalLayoutWidget);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        horizontalLayout->addWidget(buttonBox);

        horizontalLayoutWidget_2 = new QWidget(geneUI);
        horizontalLayoutWidget_2->setObjectName(QStringLiteral("horizontalLayoutWidget_2"));
        horizontalLayoutWidget_2->setGeometry(QRect(9, 10, 581, 31));
        horizontalLayout_3 = new QHBoxLayout(horizontalLayoutWidget_2);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        pushButton = new QPushButton(horizontalLayoutWidget_2);
        pushButton->setObjectName(QStringLiteral("pushButton"));

        horizontalLayout_3->addWidget(pushButton);

        lineEdit = new QLineEdit(horizontalLayoutWidget_2);
        lineEdit->setObjectName(QStringLiteral("lineEdit"));

        horizontalLayout_3->addWidget(lineEdit);

        gridLayoutWidget = new QWidget(geneUI);
        gridLayoutWidget->setObjectName(QStringLiteral("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(9, 59, 581, 331));
        gridLayout = new QGridLayout(gridLayoutWidget);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        label_5 = new QLabel(gridLayoutWidget);
        label_5->setObjectName(QStringLiteral("label_5"));

        gridLayout->addWidget(label_5, 4, 0, 1, 1);

        essential = new QCheckBox(gridLayoutWidget);
        essential->setObjectName(QStringLiteral("essential"));

        gridLayout->addWidget(essential, 6, 0, 1, 1);

        label = new QLabel(gridLayoutWidget);
        label->setObjectName(QStringLiteral("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        label_2 = new QLabel(gridLayoutWidget);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        label_4 = new QLabel(gridLayoutWidget);
        label_4->setObjectName(QStringLiteral("label_4"));

        gridLayout->addWidget(label_4, 3, 0, 1, 1);

        label_3 = new QLabel(gridLayoutWidget);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout->addWidget(label_3, 2, 0, 1, 1);

        label_6 = new QLabel(gridLayoutWidget);
        label_6->setObjectName(QStringLiteral("label_6"));

        gridLayout->addWidget(label_6, 5, 0, 1, 1);

        numID = new QLineEdit(gridLayoutWidget);
        numID->setObjectName(QStringLiteral("numID"));

        gridLayout->addWidget(numID, 0, 1, 1, 1);

        geneID = new QLineEdit(gridLayoutWidget);
        geneID->setObjectName(QStringLiteral("geneID"));

        gridLayout->addWidget(geneID, 1, 1, 1, 1);

        aaSeq = new QLineEdit(gridLayoutWidget);
        aaSeq->setObjectName(QStringLiteral("aaSeq"));

        gridLayout->addWidget(aaSeq, 2, 1, 1, 1);

        nucSeq = new QLineEdit(gridLayoutWidget);
        nucSeq->setObjectName(QStringLiteral("nucSeq"));

        gridLayout->addWidget(nucSeq, 3, 1, 1, 1);

        stability = new QLineEdit(gridLayoutWidget);
        stability->setObjectName(QStringLiteral("stability"));

        gridLayout->addWidget(stability, 4, 1, 1, 1);

        abundance = new QLineEdit(gridLayoutWidget);
        abundance->setObjectName(QStringLiteral("abundance"));

        gridLayout->addWidget(abundance, 5, 1, 1, 1);


        retranslateUi(geneUI);
        QObject::connect(buttonBox, SIGNAL(accepted()), geneUI, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), geneUI, SLOT(reject()));

        QMetaObject::connectSlotsByName(geneUI);
    } // setupUi

    void retranslateUi(QDialog *geneUI)
    {
        geneUI->setWindowTitle(QApplication::translate("geneUI", "Create a gene file", Q_NULLPTR));
        chooseDir->setText(QApplication::translate("geneUI", "Choose directory", Q_NULLPTR));
        pushButton->setText(QApplication::translate("geneUI", "Choose file", Q_NULLPTR));
        lineEdit->setText(QApplication::translate("geneUI", "(Optional, uses fasta to autocomplete)", Q_NULLPTR));
        label_5->setText(QApplication::translate("geneUI", "Stability (kcal/mol)", Q_NULLPTR));
        essential->setText(QApplication::translate("geneUI", "Essential", Q_NULLPTR));
        label->setText(QApplication::translate("geneUI", "Numeric ID", Q_NULLPTR));
        label_2->setText(QApplication::translate("geneUI", "Gene ID", Q_NULLPTR));
        label_4->setText(QApplication::translate("geneUI", "Nucleotide sequence", Q_NULLPTR));
        label_3->setText(QApplication::translate("geneUI", "Amino acid sequence", Q_NULLPTR));
        label_6->setText(QApplication::translate("geneUI", "Abundance", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class geneUI: public Ui_geneUI {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GENEUI_H
