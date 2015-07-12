#include <WindowsHeaders.h>
#include <MeshRender.h>
#include <GL/gl.h>

#define GL_SET_MTL   glEnable(GL_SMOOTH);				\
  glEnable(GL_LIGHTING);								\
  glDisable(GL_COLOR_MATERIAL);							\
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mtl.diffuse);		\
  glMaterialfv(GL_FRONT, GL_AMBIENT, mtl.ambient);		\
  glMaterialfv(GL_FRONT, GL_SPECULAR, mtl.specular);	\
  glMaterialf(GL_FRONT, GL_SHININESS, mtl.shininess);	\
  glMaterialfv(GL_FRONT, GL_EMISSION, mtl.emission);

/// @todo improve the drawing quality and speed.
void UTILITY::draw(const Objmesh& obj,const ObjMtl&mtl){
  
  const Eigen::VectorXd &verts = obj.getVerts();
  const Eigen::VectorXd &norms = obj.getVertNormal();
  const Eigen::VectorXi &faces = obj.getFaces();
  const Eigen::VectorXi &normIndex = obj.getNormalIndex();

  GL_SET_MTL;

  glBegin(GL_TRIANGLES);
  for (int f = 0; f < faces.size(); ++f){
	const int v3 = faces[f]*3;
	const int n3 = normIndex[f]*3;
	assert_in(v3,0,verts.size()-3);
	assert_in(n3,0,norms.size()-3);
	glNormal3d(norms[n3+0],norms[n3+1],norms[n3+2]);
	glVertex3d(verts[v3+0],verts[v3+1],verts[v3+2]);
  }
  glEnd();
  glDisable(GL_LIGHTING);
}

/// @todo improve the drawing quality and speed.
void UTILITY::draw(const TetMesh &tetmesh,const ObjMtl&mtl,const double *u){
  
  const VVec3d &verts = tetmesh.nodes();
  const VVec3i &faces = tetmesh.surface();
  const VVec3d &norms = tetmesh.normal();
  
  GL_SET_MTL;

  glPushAttrib(GL_POLYGON_BIT);
  glDisable(GL_LIGHTING);
  glColor3f(0.1,0.1,0.1);
  glLineWidth(1.0f);
  glPolygonMode(GL_FRONT_AND_BACK ,GL_LINE);
  glBegin(GL_TRIANGLES);

  if (u == NULL){
	  for (size_t f = 0; f < faces.size(); ++f){

		  const Vector3i &v = faces[f];
		  glNormal3d(norms[f][0],norms[f][1],norms[f][2]);
		  glVertex3d(verts[v[0]][0],verts[v[0]][1],verts[v[0]][2]);

		  glNormal3d(norms[f][0],norms[f][1],norms[f][2]);
		  glVertex3d(verts[v[1]][0],verts[v[1]][1],verts[v[1]][2]);

		  glNormal3d(norms[f][0],norms[f][1],norms[f][2]);
		  glVertex3d(verts[v[2]][0],verts[v[2]][1],verts[v[2]][2]);
	  }
  }else {
	  for (size_t f = 0; f < faces.size(); ++f){

		  const Vector3i &v = faces[f];

		  const double x0 = verts[v[0]][0]+u[v[0]*3+0];
		  const double y0 = verts[v[0]][1]+u[v[0]*3+1];
		  const double z0 = verts[v[0]][2]+u[v[0]*3+2];

		  const double x1 = verts[v[1]][0]+u[v[1]*3+0];
		  const double y1 = verts[v[1]][1]+u[v[1]*3+1];
		  const double z1 = verts[v[1]][2]+u[v[1]*3+2];

		  const double x2 = verts[v[2]][0]+u[v[2]*3+0];
		  const double y2 = verts[v[2]][1]+u[v[2]*3+1];
		  const double z2 = verts[v[2]][2]+u[v[2]*3+2];

		  glNormal3d(norms[f][0],norms[f][1],norms[f][2]);
		  glVertex3d(x0,y0,z0);

		  glNormal3d(norms[f][0],norms[f][1],norms[f][2]);
		  glVertex3d(x1,y1,z1);

		  glNormal3d(norms[f][0],norms[f][1],norms[f][2]);
		  glVertex3d(x2,y2,z2);
	  }
  }

  glEnd();
  glPopAttrib();
  
}
