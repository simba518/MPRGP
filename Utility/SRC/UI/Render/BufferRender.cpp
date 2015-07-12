#include <BufferRender.h>
#include <gl/glew.h>
#include <gl/GL.h>
#include <Objmesh.h>
#include <TetMesh.h>

BufferRender::BufferRender()
{
	m_bufferInit = false;
	m_points = NULL;
	m_normals = NULL;
	m_indices = NULL;
	objMesh = NULL;
}

void BufferRender::BuildBuffer()
{
	if (m_bufferInit)
	{
		delete[] m_points;
		delete[] m_normals;
		delete[] m_indices;
		glDeleteBuffersARB(1, &m_vertexBuffer);
		glDeleteBuffersARB(1, &m_indexBuffer);
	}
	
	glGenBuffersARB(1, &m_vertexBuffer);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexBuffer);

	const Eigen::VectorXd &verts = objMesh->getVerts();
	const Eigen::VectorXd &norms = objMesh->getVertNormal();
	const Eigen::VectorXi &faces = objMesh->getFaces();
	const Eigen::VectorXi &normIndex = objMesh->getNormalIndex();

	m_points = new double[verts.size()];
	m_normals = new double[verts.size()];
	m_indices = &faces;

	for (int f = 0; f < faces.size(); ++f){
		const int v3 = faces[f]*3;
		const int n3 = normIndex[f]*3;
		m_points[v3+0] = verts[v3+0]; m_points[v3+1] = verts[v3+1]; m_points[v3+2] = verts[v3+2];
		m_normals[v3+0] = norms[n3+0];  m_normals[v3+1] = norms[n3+1]; m_normals[v3+2] = norms[n3+2];
	}

	glBindBufferARB(GL_ARRAY_BUFFER_ARB, sizeof(GLdouble)*verts.size()*2, 0, GL_DYNAMIC_DRAW_ARB);
	glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, sizeof(GLdouble)*verts.size(), m_points);
	glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLdouble)*verts.size(), sizeof(GLdouble)*verts.size(), m_normals);

	glGenBuffersARB(1, &m_indexBuffer);
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, m_indexBuffer);
	glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, sizeof(GLdouble)*verts.size(), m_indices, GL_DYNAMIC_DRAW_ARB);

	delete[] m_points;
	m_points = NULL;
	delete[] m_normals;
	m_normals = NULL;
	delete[] m_indices;
	m_indices = NULL;
	m_bufferInit = true;

}

void BufferRender::UpdateBuffer()
{
	if (!m_bufferInit)
	{
		return;
	}

}

void BufferRender::DrawBuffer()
{

}
