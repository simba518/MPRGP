#include <iostream>
#include <sstream>
#include "QInputEventRecorderWidget.h"
using namespace std;
using namespace UTILITY;

QInputEventRecorderWidget::QInputEventRecorderWidget
(QObject *obj,QInputEventRecorderObserver*ob,const string fname,
 QWidget *parent,Qt::WFlags flags):QWidget(parent,flags), record_filename(fname){

  dialog.setupUi(this);
  createConnections();
  recorder.load(fname.c_str());
  updateTimeOffset();
  dialog.file_path->setText(record_filename.c_str());

  setWindowFlags(Qt::WindowMinimizeButtonHint);
  recorder.setObj(obj);
  recorder.setObserver(ob);
}

void QInputEventRecorderWidget::createConnections(){

  connect(dialog.pushButton_load, SIGNAL(clicked()),this,SLOT(load()));
  connect(dialog.pushButton_save, SIGNAL(clicked()),this,SLOT(save()));
  connect(dialog.pushButton_saveAs, SIGNAL(clicked()),this,SLOT(saveAs()));
  connect(dialog.pushButton_cut,SIGNAL(clicked()),&recorder,SLOT(removeTheUnplayedEvents()));
  connect(dialog.pushButton_clear,SIGNAL(clicked()),&recorder,SLOT(clear()));

  connect(dialog.radioButton_pause, SIGNAL(clicked()), &recorder, SLOT(stop()));
  connect(dialog.radioButton_play, SIGNAL(clicked()), &recorder, SLOT(replay()));
  connect(dialog.radioButton_record, SIGNAL(clicked()), &recorder, SLOT(record()));

  connect(dialog.checkBox_observer, SIGNAL(clicked()), &recorder, SLOT(toggleObserver()));

  connect(&recorder,SIGNAL(isPaused(bool)),dialog.radioButton_pause,SLOT(setChecked(bool)));
  connect(&recorder,SIGNAL(updateTimeOffset()),this,SLOT(updateTimeOffset()));
}

void QInputEventRecorderWidget::updateTimeOffset(){
  
  const int total = recorder.totalTimeOffset();
  const int current = recorder.currentTimeOffset();

  dialog.progressBar->setRange(0,total>0?total:1);
  dialog.progressBar->setValue(current);

  stringstream ss;
  ss.precision(5);
  ss << (current/1000.0f) << " / " << (total/1000.0f) << "s";
  const string ssd = ss.str(); 
  dialog.current_total->setText(ssd.c_str());
}

void QInputEventRecorderWidget::load(){

  recorder.stop();
  QString filename = QFileDialog::getOpenFileName(NULL);
  if (!filename.isEmpty()){
	recorder.load(filename);
	record_filename = filename.toStdString();
	dialog.file_path->setText(record_filename.c_str());
  }
}

void QInputEventRecorderWidget::save(){
  recorder.save(record_filename.c_str());
}

void QInputEventRecorderWidget::saveAs(){

  recorder.stop();
  const QString filename = QFileDialog::getSaveFileName(NULL);
  if (!filename.isEmpty()){
	recorder.save(filename);
	record_filename = filename.toStdString();
	dialog.file_path->setText(record_filename.c_str());
  }
}
