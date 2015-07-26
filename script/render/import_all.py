import maya.cmds as cmds

def selectDir():
	basicFilter = "*.obj"
	filename=cmds.fileDialog2(fileMode=4,ds=1,caption="Select Sequence",fileFilter=basicFilter)
	uniToStr=str(filename)
	global objList
	filePath=uniToStr.split("'")
	search=".obj"
	objList=[mesh for mesh in filePath if search in mesh]
	print objList
	cmds.textScrollList(path, edit=True, append=objList)

def actualImporter(myobjList, group_name):
    mycurrent_scene=cmds.ls(assemblies=True)
    myobjList.sort()
    for each in myobjList:
        print each
        loaded = cmds.file(each, i=True, dns=True)
    return getMeshes(group_name, mycurrent_scene)

def importer():
    actualImporter(objList, "GROUP")

def loadWin():
	if cmds.window('objwindow',exists=True):
		cmds.deleteUI('objwindow',window=True)
	win=cmds.window('objwindow',title="Load Obj sequence",widthHeight=(250,450))
	cmds.columnLayout()
	cmds.text('Path to files')
	global path
	path=cmds.textScrollList(w=250)
	browseButton=cmds.button(label='browse',w=250,c='selectDir()')
	importButton=cmds.button(label='import',w=250,c='importer()')
	autoImportButton=cmds.button(label='import all',w=250,c='loadAllMeshes()')
	cmds.showWindow(win)

def getMeshes(group_name, mycurrent_scene):
	newScene=cmds.ls(assemblies=True)
	meshList=[mesh for mesh in newScene if mesh not in mycurrent_scene]
	for i in range(len(meshList)):
		cmds.setKeyframe(meshList[i],attribute='visibility',value=0,t=0)
		cmds.setKeyframe(meshList[i],attribute='visibility',value=1,t=i+1)
		cmds.setKeyframe(meshList[i],attribute='visibility',value=0,t=i+2)
	return cmds.group(meshList,n=group_name)
	

def indexToStr(index):
    if index < 10:
        return '00'+str(index)
    elif index < 100:
        return '0'+str(index)
    return str(index)

def getObjList(bath_path, begin, end):
    objfileList = []
    for index in range(begin, end+1):
        objfileList.append(bath_path+indexToStr(index)+".obj")
    return objfileList

def loadAllMeshes():
    base_dir = 'E:\\\\lsw\\\\objfiles\\\\'
    my_shaders = ['lambert21','lambert22','lambert11','lambert21','lambert19','lambert18','lambert20','lambert17']
    for i in range (0,8):
        myobjList = getObjList(base_dir+'obj_'+str(i)+'\\\\obj_'+str(i)+'_frame_',399,399)
        cmds.textScrollList(path, edit=True, append=myobjList)
        group_i = actualImporter(myobjList, "GROUP")
        cmds.select(group_i)
        cmds.hyperShade(assign=my_shaders[i])

loadWin()
