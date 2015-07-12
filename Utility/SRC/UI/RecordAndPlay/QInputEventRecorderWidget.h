#ifndef _QINPUTEVENTRECORDERWIDGET_H_
#define _QINPUTEVENTRECORDERWIDGET_H_

#include <QFileDialog>
#include <QKeyEvent>
#include <QInputEventRecorder.h>
#include <string>
#include "ui_record_replay.h"
using namespace std;

namespace UTILITY{
  
  /**
   * @class QInputEventRecorderWidget provide the record/replay application with
   * a control panel.
   * @ ussage exmaple:
    
   #include <qapplication.h>
   #include <QGLViewerExt.h>
   #include <QInputEventRecorderWidget.h>
   using namespace QGLVEXT;
   using namespace UTILITY;
 
   int main(int argc, char *argv[]){

   QApplication application(argc,argv);
   QWidget window;
   QGLViewerExt viewer(&window);
   window.show();

   QInputEventRecorderWidget recorder(&window);
   recorder.show();

   return  application.exec();
   }

   */
  class QInputEventRecorderWidget:public QWidget{

	Q_OBJECT
	
  public:
	QInputEventRecorderWidget(QObject *obj=0,QInputEventRecorderObserver*ob=0,
							  const string fname="tempt_qt_event_record.b",
							  QWidget *parent=0,Qt::WFlags flags=0);

	void setObj(QObject* widget){
	  recorder.setObj(widget);
	}

  protected:
	void createConnections();
	
  protected slots:
	void load();
	void save();
	void saveAs();
	void updateTimeOffset();

  private:
	string record_filename;
	QInputEventRecorder recorder;
	Ui_record_replay_dialog dialog;
  };
  
}//end of namespace

#endif /*_QINPUTEVENTRECORDERWIDGET_H_*/
