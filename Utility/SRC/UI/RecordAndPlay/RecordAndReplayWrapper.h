#ifndef RECORDANDREPLAYWRAPPER_H
#define RECORDANDREPLAYWRAPPER_H

#include <QInputEventRecorder.h>
#include <QObject>
#include <QtGui/QMainWindow>
#include <QWidget>
#include <boost/shared_ptr.hpp>

namespace UTILITY
{
	class RecordAndReplayWrapper: public QObject
	{
		Q_OBJECT

	public:

		RecordAndReplayWrapper(QWidget* widget);
		~RecordAndReplayWrapper(){}

		public slots:
			void record();
			void replay();
			void stop();
			void save();
			void load();

	protected:
	private:
		QWidget* widget; 
		QInputEventRecorder m_recorder;
	};
	typedef boost::shared_ptr<RecordAndReplayWrapper> pRecordAndReplayWrapper;
}
#endif