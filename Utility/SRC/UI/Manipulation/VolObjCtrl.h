#ifndef _VOLOBJCTRL_H_
#define _VOLOBJCTRL_H_

#include <boost/shared_ptr.hpp>
#include <QObject>
#include <QtGui/QMainWindow>
#include <TetMeshEmbeding.h>
#include <FileDialog.h>
#include <JsonFilePaser.h>
using namespace QGLVEXT;

namespace UTILITY{
  
  /**
   * @class VolObjCtrl manager the VolObjMesh object: load obj mesh, volume
   * mesh, interpolate weights, and precompute interpolate weights.
   * 
   */
  class VolObjCtrl: public QObject{

	Q_OBJECT

  public: 
	VolObjCtrl(QMainWindow *main_win, pTetMeshEmbeding volobj):volobjmesh(volobj){
  	  file_dialog = pFileDialog(new FileDialog(main_win));
	}
	bool initialize(const string filename){
	  if (!volobjmesh)
	  	volobjmesh = pTetMeshEmbeding(new TetMeshEmbeding());
	  const bool succ = init(filename);
	  sendMsgResetScene();
	  return succ;
	}
	pTetMeshEmbeding getVolObjMesh(){
	  return volobjmesh;
	}
	pTetMeshEmbeding_const getVolObjMesh()const{
	  return volobjmesh;
	}
	pTetMesh_const getVolMesh()const{
	  assert (volobjmesh != NULL);
	  return volobjmesh->getTetMesh();
	}
	pObjmesh_const getObjMesh()const{
	  assert (volobjmesh != NULL);
	  return volobjmesh->getObjMesh();
	}
	
  public slots:
	void loadObjMesh(){
	  const string filename = file_dialog->load();
	  if(filename.size() >0){
		file_dialog->warning(volobjmesh != NULL && volobjmesh->loadObjMesh(filename));
		sendMsgResetScene();
	  }
	}
	void loadVolMesh(){
	  const string filename = file_dialog->load();
	  if(filename.size() >0){
		file_dialog->warning(volobjmesh != NULL && volobjmesh->loadTetMesh(filename));
		sendMsgResetScene();
	  }
	}
	void loadInterpWeights(){
	  const string filename = file_dialog->load();
	  if(filename.size() >0){
		file_dialog->warning(volobjmesh != NULL && volobjmesh->loadWeights(filename));
	  }
	}
	void saveInterpWeights()const{
	  const string filename = file_dialog->save();
	  if(filename.size() >0){
		file_dialog->warning(volobjmesh != NULL && volobjmesh->writeWeights(filename));
	  }
	}

	void precomputeInterpWeights(){
	  if(volobjmesh != NULL){
		volobjmesh->buildInterpWeights();
	  }
	}
	void sendMsgResetScene()const{
	  if(volobjmesh){
		double min[3],max[3];
		volobjmesh->getBBox(min,max);
		emit resetSceneMsg(min[0],min[1],min[2],max[0],max[1],max[2]);
	  }
	}

  protected:
	bool init(const string init_filename){

	  UTILITY::JsonFilePaser inf;
	  bool succ = false;
	  if ( inf.open(init_filename) ){
		string filename;
		if ( inf.readFilePath("obj_file",filename,true)){
		  succ = volobjmesh->loadObjMesh(filename);
		}
		if ( inf.readFilePath("vol_file",filename,true) ){
		  succ &= volobjmesh->loadTetMesh(filename);
		}
		if ( inf.readFilePath("vol2obj_weights",filename,false)){
		  if(!volobjmesh->loadWeights(filename)&&succ){
			volobjmesh->buildInterpWeights();
			volobjmesh->writeWeights(filename);
		  }
		}else if(succ){
		  volobjmesh->buildInterpWeights();
		}
	  }
	  return succ;
	}
 
  signals:
	// (x0,y0,z0) is the corner with the smallest coordinates of the bounding box.
	// (x1,y1,z1) is the corner with the largest coordinates of the bounding box.
	void resetSceneMsg(double x0,double y0,double z0,double x1,double y1,double z1)const;
	void update();
	
  private:
	pFileDialog file_dialog;
	pTetMeshEmbeding volobjmesh;
  };
  
  typedef boost::shared_ptr<VolObjCtrl> pVolObjCtrl;
  
}//end of namespace

#endif /* _VOLOBJCTRL_H_ */
