#include "mainConfig.h"
#ifdef MAIN_GUI

#include "PVSMScene.h"
#include "FPSCamera.h"
#include "StaticDraw.h"
#include "DeformEditor.h"
#include "FEMEmbeddedMesh.h"
#include <GL/freeglut.h>

USE_PRJ_NAMESPACE
	
boost::shared_ptr<PVSMScene> scene;
boost::shared_ptr<FEMSolver> solver;

void display(){scene->render();}
void idle(){scene->idle();}
void reshape(int w,int h){scene->reshape(w,h);}
void timer(int id)
{
	scene->frame();
	glutTimerFunc(1000/scene->getFPS(),timer,0);
}
void mouse(int button,int state,int x,int y){scene->mouse(button,state,x,y);}
void motion(int x,int y){scene->motion(x,y);}
void passiveMotion(int x,int y){scene->passiveMotion(x,y);}
void mouseWheel(int wheel,int dir,int x,int y){scene->mouseWheel(wheel,dir,x,y);}
void keyboard(unsigned char key,int x,int y){scene->keyboard(key,x,y);}
void keyboardUp(unsigned char key,int x,int y){scene->keyboardUp(key,x,y);}

void init()
{
	//setup renderer
	//scene->setRenderer(boost::shared_ptr<DefaultRenderer>(new DefaultRenderer));
	//scene->setRenderer(boost::shared_ptr<DefaultRenderer>(new PerPixelRenderer));
	scene->setRenderer(boost::shared_ptr<DefaultRenderer>(new VSMPerPixelRenderer));
	scene->getRenderer().setShowLights(true,0.01f);

	boost::shared_ptr<PVSMScene::PVSMAction> CAM(new FPSCamera());
	scene->addAction(CAM);
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new StaticObjMesh()));
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new StaticBox()));
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new StaticPlane()));
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new StaticSphere()));
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new StaticLine()));
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new StaticPoint()));
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new Light()));

	boost::shared_ptr<PVSMScene::PVSMAction> FEM(new DeformEditor(*solver));
	scene->addAction(FEM);
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new FEMEmbeddedMesh(solver->getMesh())));
	scene->addAction(boost::shared_ptr<PVSMScene::PVSMAction>(new FEMComputeMesh(solver->getMesh())));
	scene->addSceneAction(boost::shared_ptr<PVSMSwitch>(new PVSMSwitch(CAM,FEM)));
}
int main(int argc,char** argv)
{
	glutInit(&argc,argv);
	glewInit();
	scene.reset(new PVSMScene);
	solver.reset(new ReducedFEMSolver(3));

	ASSERT_MSG(argc >= 2,"Usage: exe [SceneFile].xml")
	init();
	ASSERT_MSG(scene->read(std::string(argv[1])),"Error Reading File!")
	{
		boost::filesystem::path path(argv[1]);
		path.replace_extension(".xml");
		scene->clearUnused();
		scene->write(path.string());
	}
	PVSMScene::GlobalSetting settings=scene->getGlobalSetting();

	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutInitWindowSize(settings._size[0],settings._size[1]);
	glutInitWindowPosition(settings._pos[0],settings._pos[1]);

	glutCreateWindow(settings._title._value.c_str());
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutReshapeFunc(reshape);
	glutTimerFunc(0,timer,0);

	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutPassiveMotionFunc(passiveMotion);
	glutMouseWheelFunc(mouseWheel);

	glutKeyboardFunc(keyboard);
	glutKeyboardUpFunc(keyboardUp);
	
	scene->init();
	glutMainLoop();
	return 0;
}

#endif