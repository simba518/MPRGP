#ifndef _FILEDIALOG_H_
#define _FILEDIALOG_H_

#include <string>
#include <boost/shared_ptr.hpp>
#include <QtGui/QMainWindow>
#include <QtGui/QFileDialog>

using std::string;

namespace QGLVEXT{

  typedef boost::shared_ptr<QFileDialog> pQFileDialog;
  
  /**
   * @class FileDialog Qt dialog for loading and saving files.
   * 
   */
  class FileDialog{
	
  public:
	FileDialog(QMainWindow * main_win);
	
	string load(const string appendix="", const string title = "", const string default_file_name = ""); 
	string save(const string appendix="", const string title = "", const string default_file_name = "");
	void warning(const bool succ = false, const string message = "failed!") const;

	static void setDefaultDir(const string dir);
	
  protected:
	QMainWindow *main_win;
	pQFileDialog fileDialog;
	string last_opened_directory;
	static string default_file_dialog_dir;
  };
  
  typedef boost::shared_ptr<FileDialog> pFileDialog;
  
}//end of namespace

#endif /*_FILEDIALOG_H_*/
