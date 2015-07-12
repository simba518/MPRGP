/********************************************************************************
** Form generated from reading UI file 'record_replay.ui'
**
** Created: Fri Apr 4 13:26:46 2014
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RECORD_REPLAY_H
#define UI_RECORD_REPLAY_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_record_replay_dialog
{
public:
    QVBoxLayout *verticalLayout_2;
    QGridLayout *gridLayout_2;
    QGridLayout *gridLayout;
    QPushButton *pushButton_saveAs;
    QPushButton *pushButton_clear;
    QPushButton *pushButton_cut;
    QPushButton *pushButton_load;
    QPushButton *pushButton_save;
    QVBoxLayout *verticalLayout_3;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout;
    QRadioButton *radioButton_pause;
    QRadioButton *radioButton_play;
    QRadioButton *radioButton_record;
    QCheckBox *checkBox_observer;
    QGridLayout *gridLayout_3;
    QProgressBar *progressBar;
    QLabel *current_total;
    QLabel *file_path;
    QLabel *label;

    void setupUi(QWidget *record_replay_dialog)
    {
        if (record_replay_dialog->objectName().isEmpty())
            record_replay_dialog->setObjectName(QString::fromUtf8("record_replay_dialog"));
        record_replay_dialog->resize(360, 240);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(record_replay_dialog->sizePolicy().hasHeightForWidth());
        record_replay_dialog->setSizePolicy(sizePolicy);
        record_replay_dialog->setMinimumSize(QSize(360, 240));
        record_replay_dialog->setMaximumSize(QSize(360, 240));
        verticalLayout_2 = new QVBoxLayout(record_replay_dialog);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        pushButton_saveAs = new QPushButton(record_replay_dialog);
        pushButton_saveAs->setObjectName(QString::fromUtf8("pushButton_saveAs"));

        gridLayout->addWidget(pushButton_saveAs, 2, 0, 1, 1);

        pushButton_clear = new QPushButton(record_replay_dialog);
        pushButton_clear->setObjectName(QString::fromUtf8("pushButton_clear"));

        gridLayout->addWidget(pushButton_clear, 6, 0, 1, 1);

        pushButton_cut = new QPushButton(record_replay_dialog);
        pushButton_cut->setObjectName(QString::fromUtf8("pushButton_cut"));

        gridLayout->addWidget(pushButton_cut, 3, 0, 1, 1);

        pushButton_load = new QPushButton(record_replay_dialog);
        pushButton_load->setObjectName(QString::fromUtf8("pushButton_load"));

        gridLayout->addWidget(pushButton_load, 0, 0, 1, 1);

        pushButton_save = new QPushButton(record_replay_dialog);
        pushButton_save->setObjectName(QString::fromUtf8("pushButton_save"));

        gridLayout->addWidget(pushButton_save, 1, 0, 1, 1);


        gridLayout_2->addLayout(gridLayout, 3, 0, 1, 1);

        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        groupBox = new QGroupBox(record_replay_dialog);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        verticalLayout = new QVBoxLayout(groupBox);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        radioButton_pause = new QRadioButton(groupBox);
        radioButton_pause->setObjectName(QString::fromUtf8("radioButton_pause"));
        radioButton_pause->setChecked(true);
        radioButton_pause->setAutoRepeat(false);

        verticalLayout->addWidget(radioButton_pause);

        radioButton_play = new QRadioButton(groupBox);
        radioButton_play->setObjectName(QString::fromUtf8("radioButton_play"));

        verticalLayout->addWidget(radioButton_play);

        radioButton_record = new QRadioButton(groupBox);
        radioButton_record->setObjectName(QString::fromUtf8("radioButton_record"));

        verticalLayout->addWidget(radioButton_record);


        verticalLayout_3->addWidget(groupBox);

        checkBox_observer = new QCheckBox(record_replay_dialog);
        checkBox_observer->setObjectName(QString::fromUtf8("checkBox_observer"));
        checkBox_observer->setChecked(false);
        checkBox_observer->setTristate(false);

        verticalLayout_3->addWidget(checkBox_observer);


        gridLayout_2->addLayout(verticalLayout_3, 3, 1, 1, 1);


        verticalLayout_2->addLayout(gridLayout_2);

        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        progressBar = new QProgressBar(record_replay_dialog);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setValue(24);

        gridLayout_3->addWidget(progressBar, 0, 0, 1, 1);

        current_total = new QLabel(record_replay_dialog);
        current_total->setObjectName(QString::fromUtf8("current_total"));

        gridLayout_3->addWidget(current_total, 0, 1, 1, 1);


        verticalLayout_2->addLayout(gridLayout_3);

        file_path = new QLabel(record_replay_dialog);
        file_path->setObjectName(QString::fromUtf8("file_path"));

        verticalLayout_2->addWidget(file_path);

        label = new QLabel(record_replay_dialog);
        label->setObjectName(QString::fromUtf8("label"));
        QFont font;
        font.setPointSize(8);
        font.setItalic(true);
        label->setFont(font);
        label->setTextFormat(Qt::RichText);
        label->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        verticalLayout_2->addWidget(label);


        retranslateUi(record_replay_dialog);

        QMetaObject::connectSlotsByName(record_replay_dialog);
    } // setupUi

    void retranslateUi(QWidget *record_replay_dialog)
    {
        record_replay_dialog->setWindowTitle(QApplication::translate("record_replay_dialog", "record&replay", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        record_replay_dialog->setToolTip(QString());
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        record_replay_dialog->setStatusTip(QString());
#endif // QT_NO_STATUSTIP
#ifndef QT_NO_WHATSTHIS
        record_replay_dialog->setWhatsThis(QString());
#endif // QT_NO_WHATSTHIS
        pushButton_saveAs->setText(QApplication::translate("record_replay_dialog", "saveAs", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_clear->setToolTip(QApplication::translate("record_replay_dialog", "remove all operations.", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        pushButton_clear->setStatusTip(QString());
#endif // QT_NO_STATUSTIP
#ifndef QT_NO_WHATSTHIS
        pushButton_clear->setWhatsThis(QString());
#endif // QT_NO_WHATSTHIS
        pushButton_clear->setText(QApplication::translate("record_replay_dialog", "clear", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_cut->setToolTip(QApplication::translate("record_replay_dialog", "remove the remainning operations.", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_WHATSTHIS
        pushButton_cut->setWhatsThis(QString());
#endif // QT_NO_WHATSTHIS
#ifndef QT_NO_ACCESSIBILITY
        pushButton_cut->setAccessibleDescription(QString());
#endif // QT_NO_ACCESSIBILITY
        pushButton_cut->setText(QApplication::translate("record_replay_dialog", "cut", 0, QApplication::UnicodeUTF8));
        pushButton_load->setText(QApplication::translate("record_replay_dialog", "load", 0, QApplication::UnicodeUTF8));
        pushButton_save->setText(QApplication::translate("record_replay_dialog", "save", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QString());
        radioButton_pause->setText(QApplication::translate("record_replay_dialog", "pause", 0, QApplication::UnicodeUTF8));
        radioButton_play->setText(QApplication::translate("record_replay_dialog", "play", 0, QApplication::UnicodeUTF8));
        radioButton_record->setText(QApplication::translate("record_replay_dialog", "record", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        checkBox_observer->setToolTip(QApplication::translate("record_replay_dialog", "send message to observer.", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        checkBox_observer->setStatusTip(QString());
#endif // QT_NO_STATUSTIP
#ifndef QT_NO_ACCESSIBILITY
        checkBox_observer->setAccessibleDescription(QString());
#endif // QT_NO_ACCESSIBILITY
        checkBox_observer->setText(QApplication::translate("record_replay_dialog", "observer", 0, QApplication::UnicodeUTF8));
        current_total->setText(QApplication::translate("record_replay_dialog", "current/total", 0, QApplication::UnicodeUTF8));
        file_path->setText(QApplication::translate("record_replay_dialog", "file:", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("record_replay_dialog", "author: Simba,  email: lisiwangcg@gmail.com", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class record_replay_dialog: public Ui_record_replay_dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RECORD_REPLAY_H
