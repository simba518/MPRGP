#ifndef _BUFFERRENDER_H
#define _BUFFERRENDER_H

namespace UTILITY{

class BufferRender
{
public:
	//add by qnn
	void BuildBuffer();

	void UpdateBuffer();

	void DrawBuffer();

	void SetObjMesh(pObjmesh _objMesh)
	{
		objMesh = _objMesh;
	}

protected:
private:
	
	bool m_bufferInit; // VBO init

	unsigned int m_vertexBuffer;// vertex buffer

	unsigned int m_indexBuffer;// index buffer 

	double * m_points;

	double * m_normals;

	unsigned int * m_indices;

	pObjmesh objMesh;
};

}

#endif 
