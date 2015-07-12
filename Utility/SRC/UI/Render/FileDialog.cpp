#include <QObject>
#include <QtGui/QMessageBox>
#include "FileDialog.h"

using namespace QGLVEXT;

string FileDialog::default_file_dialog_dir;

FileDialog::FileDialog(QMainWindow * _main_win){

  main_win = _main_win;  
  fileDialog = pQFileDialog (new QFileDialog(_main_win));
}

string FileDialog::load(const string appendix, const string title, const string default_file_name) {
  
  string temp_app;
  if(appendix.size() > 0){

	temp_app = "(*.";
	temp_app += appendix;
	temp_app += ");;All Files (*)";
  }else{
	temp_app = "All Files (*)";
  }

  QString q_title = QObject::tr(title.c_str());
  QString q_fname = QObject::tr((default_file_name).c_str());
  if(q_fname.size() <= 0){
	q_fname = last_opened_directory.c_str();
  }
  if(q_fname.size() <= 0){
	q_fname = default_file_dialog_dir.c_str();
  }

  QString q_appendix = QObject::tr(temp_app.c_str());
  QString filename = fileDialog->getOpenFileName(main_win,q_title,q_fname,q_appendix);

  last_opened_directory = filename.toStdString(); 
  return filename.toStdString();
}

string FileDialog::save(const string appendix, const string title, const string default_file_name) {
  
  string temp_app;
  if(appendix.size() > 0){

	temp_app = "(*.";
	temp_app += appendix;
	temp_app += ");;All Files (*)";
  }else{
	temp_app = "All Files (*)";
  }

  QString q_title = QObject::tr(title.c_str());
  QString q_fname = QObject::tr((default_file_name).c_str());
  if(q_fname.size() <= 0){
	q_fname = last_opened_directory.c_str();
  }

  if(q_fname.size() <= 0){
	q_fname = default_file_dialog_dir.c_str();
  }

  QString q_appendix = QObject::tr(temp_app.c_str());
  QString filename = fileDialog->getSaveFileName(main_win,q_title,q_fname,q_appendix);
  last_opened_directory = filename.toStdString(); 
  return filename.toStdString();
}

void FileDialog::warning(const bool succ, const string message) const{
  
  if(!succ){
	QMessageBox::warning(main_win,main_win->windowTitle(),message.c_str());
  }
}

void FileDialog::setDefaultDir(const string dir){
  
  default_file_dialog_dir = dir;
}
