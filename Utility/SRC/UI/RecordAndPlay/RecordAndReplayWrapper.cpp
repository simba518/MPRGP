#include <RecordAndReplayWrapper.h>
#include <QFileDialog>
#include <QString>
using namespace UTILITY;

RecordAndReplayWrapper::RecordAndReplayWrapper(QWidget* _widget)
{
	widget = _widget;
	m_recorder.setObj(widget);
}

void RecordAndReplayWrapper::record()
{
	m_recorder.record();
}

void RecordAndReplayWrapper::replay()
{
	m_recorder.replayScaled(1.0);
}

void RecordAndReplayWrapper::save()
{
	QString filename = QFileDialog::getSaveFileName(widget);
	if (!filename.isEmpty())
	{
		m_recorder.save(filename);
	}
}

void RecordAndReplayWrapper::stop()
{
	m_recorder.stop();
}

void RecordAndReplayWrapper::load()
{
	stop();
	QString filename = QFileDialog::getOpenFileName(widget);
	if (!filename.isEmpty())
	{
		m_recorder.load(filename);
	}
}
