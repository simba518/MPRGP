#include <qapplication.h>
#include <QGLViewerExt.h>
#include <QInputEventRecorderCmd.h>
#include "QtEventRecorder.h"
#include "SimpleMeshRender.h"
#include "DefManipulator.h"
#include <QPushButton>
#include <QInputEventRecorderWidget.h>
using namespace QGLVEXT;
using namespace UTILITY;
using namespace Eigen;

#ifdef WIN32
#define TEST_DATA_DIR "D:\\lsw\\Utility\\TestCase\\TestData\\"
#endif

int main(int argc, char** argv){

  QApplication application(argc,argv);  
  QWidget window;
  QGLViewerExt viewer(&window);
  QInputEventRecorderWidget recorder;
  recorder.setObj(&window);
  recorder.show();

  // QInputEventRecorderCmd recorder(&window);

  // { // load meshses
  // 	// const string fn = string(TEST_DATA_DIR)+"beam.obj";
  // 	// Objmesh obj;
  // 	// obj.load("D:\\lsw\\mesh1.obj");
  // 	// viewer.addSelfRenderEle(pSelfRenderEle(new ObjRender(obj)));

  // 	// const string fn2 = string(TEST_DATA_DIR)+"dino.abq";
  // 	// TetMesh tet;
  // 	// tet.load(fn2);
  // 	// const VectorXd u = VectorXd::Ones(tet.nodes().size()*3)*5.0f;
  // 	// viewer.addSelfRenderEle(pSelfRenderEle(new TetRender(tet,u)));
  // 	// viewer.addSelfRenderEle(pSelfRenderEle(new TetRender(tet)));
  // }

  // { // parse record command
  // 	// if (2 == argc)
  // 	//   recorder.setCmd(argv[argc-1]);
  // 	// else if(3 == argc)
  // 	//   recorder.setCmd(argv[argc-2],argv[argc-1]);
  // }

  // // test manipulation

  pLocalframeManipulatoion mani = pLocalframeManipulatoion(new DefManipulatorExt(&viewer));
  LocalframeManipulatoionCtrl mani_ctrl(&viewer, mani);

  viewer.show3DGrid();
  viewer.resize(600,600);
  viewer.setWindowTitle("simpleViewer");
  window.show();

  return application.exec();
}
