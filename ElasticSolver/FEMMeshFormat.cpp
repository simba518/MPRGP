#include "FEMMeshFormat.h"
#include "FEMMesh.h"
#include "FEMUtils.h"
#include "FEMCollision.h"
#include "ImageOp.h"
#include <boost/filesystem/fstream.hpp>
#include <boost/tokenizer.hpp>

USE_PRJ_NAMESPACE

void FEMMeshFormat::meshToABQ(std::istream& is,std::ostream& os)
{
    Vec4i tets;
    Vec3d pos;
    sizeType nr,one;
    string line;

    getline(is,line);
    istringstream(line) >> nr;
    os << "*NODE" << std::endl;
    for(sizeType i=0; i<nr; i++) {
        getline(is,line);
        istringstream iss(line);
        iss >> pos[0] >> pos[1] >> pos[2];
        os << (i+1) << "," <<
           pos[0] << "," << pos[1] << "," << pos[2] << std::endl;
    }

    os << "*ELEMENT, type=C3D4, ELSET=PART1" << std::endl;
    getline(is,line);
    istringstream(line) >> nr;
    for(sizeType i=0; i<nr; i++) {
        getline(is,line);
        istringstream iss(line);
        iss >> one >> tets[0] >> tets[1] >> tets[2] >> tets[3];
        os << (i+1) << "," <<
           tets[0] << "," << tets[1] << "," <<
           tets[2] << "," << tets[3] << std::endl;
    }

    os << "*ELSET,ELSET=EALL,GENERATE" << std::endl;
    os << "1," << nr << std::endl;
}
void FEMMeshFormat::meshToABQ(const std::string& is,const std::string& os)
{
    boost::filesystem::ifstream inf(is);
    boost::filesystem::ofstream outf(os);
    meshToABQ(inf,outf);
}
void FEMMeshFormat::VTKToABQ(std::istream& is,std::ostream& os)
{
    int nr;
    std::string line;
    while(getline(is,line).good()) {
        if(beginsWith(line,"POINTS")) {
#ifdef _MSC_VER
            sscanf_s(line.c_str(),"POINTS %d",&nr);
#else
            sscanf(line.c_str(),"POINTS %d",&nr);
#endif
            vector<scalarD> poses(nr*3);
            os << "*NODE" << std::endl;
            for(sizeType i=0; i<nr*3;) {
                getline(is,line);
                boost::char_separator<char> sep(" ");
                boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
                boost::tokenizer<boost::char_separator<char> >::iterator tokIter;
                for(tokIter=tokens.begin(); tokIter!=tokens.end(); ++tokIter)
                    istringstream(*tokIter)>>poses[i++];
            }
            for(sizeType i=0; i<nr*3; i+=3)
                os << (i/3+1) << ',' <<
                   poses[i+0] << ',' <<
                   poses[i+1] << ',' <<
                   poses[i+2] << std::endl;
        } else if(beginsWith(line,"CELLS")) {
            sizeType tmp;
            Vec4i tet;
#ifdef _MSC_VER
            sscanf_s(line.c_str(),"CELLS %d",&nr);
#else
            sscanf(line.c_str(),"CELLS %d",&nr);
#endif
            os << "*ELEMENT, type=C3D4, ELSET=PART1" << std::endl;
            for(sizeType i=0; i<nr; i++) {
                getline(is,line);
                istringstream(line) >> tmp >> tet[0] >> tet[1] >> tet[2] >> tet[3];
                os << (i+1) << "," <<
                   (tet[0]+1) << "," << (tet[1]+1) << "," <<
                   (tet[2]+1) << "," << (tet[3]+1) << std::endl;
            }
        }
    }
    os << "*ELSET,ELSET=EALL,GENERATE" << std::endl;
    os << "1," << nr << std::endl;
}
void FEMMeshFormat::VTKToABQ(const std::string& is,const std::string& os)
{
    boost::filesystem::ifstream inf(is);
    boost::filesystem::ofstream outf(os);
    VTKToABQ(inf,outf);
}
void FEMMeshFormat::ABQToObj(const std::string& is,std::ostream& os)
{
    FEMMesh mesh(3,boost::shared_ptr<FEMCollision>(new FEMCollision));
    mesh.reset(is,0.0f);
    mesh.getB(0).writeObj(os);
}
void FEMMeshFormat::ABQToObj(const std::string& is,const std::string& os)
{
    boost::filesystem::ofstream outf(os);
    ABQToObj(is,outf);
}
void createCell
(vector<Vec4i,Eigen::aligned_allocator<Vec4i> >& tss,
 const boost::unordered_map<Vec3i,sizeType,Hash>& vMap,
 const Vec3i& base,sizeType GIVEN,sizeType II,sizeType IJ,sizeType I0)
{
#define GI(I) vMap.find(I+base)->second
#define ADDTWO	\
tss.push_back(Vec4i(ctr,cid,GI(a),GI((a+b)/2)));	\
tss.push_back(Vec4i(ctr,cid,GI(b),GI((a+b)/2)));
    sizeType ctr=GI(Vec3i::Ones());
    Vec3i a,b,c;
    //set fixed
    a[I0]=GIVEN;
    b[I0]=GIVEN;
    //set c
    c[I0]=GIVEN;
    c[II]=1;
    c[IJ]=1;
    sizeType cid=GI(c);
    //two tet
    a[II]=0;
    a[IJ]=0;
    b[II]=0;
    b[IJ]=2;
    ADDTWO
    //two tet
    a[II]=0;
    a[IJ]=2;
    b[II]=2;
    b[IJ]=2;
    ADDTWO
    //two tet
    a[II]=2;
    a[IJ]=2;
    b[II]=2;
    b[IJ]=0;
    ADDTWO
    //two tet
    a[II]=2;
    a[IJ]=0;
    b[II]=0;
    b[IJ]=0;
    ADDTWO
}
void findReplace(Vec3i& t,int in,int out)
{
    if(t[0] == in)t[0]=out;
    else if(t[1] == in)t[1]=out;
    else if(t[2] == in)t[2]=out;
}
void FEMMeshFormat::segmentObj(scalar deg,ObjMesh& mesh)
{
    typedef ObjMesh::Edge EDGE;
    typedef ObjMesh::EdgeMap EMAP;
    typedef map<pair<int,int>,EDGE,EMAP::LSS> EMAPM;

    //build edge, find face normal
    EMAP eMap;
    mesh.smooth();
    mesh.buildEdge(eMap);
    //joint patch
    DisjointSet<int> dset(mesh.getI().size());
    for(EMAPM::const_iterator beg=eMap._ess.begin(),end=eMap._ess.end(); beg!=end; beg++)
        if(beg->second._tris.size() == 2) {
            const vector<int>& tris=beg->second._tris;
            Vec3 na=mesh.getTN()[tris[0]];
            Vec3 nb=mesh.getTN()[tris[1]];
            if(getAngle3D<scalar>(na,nb) < deg)
                dset.joinSafe(tris[0],tris[1]);
        }
    //build neighbor triangle
    vector<vector<int> > tNeigh(mesh.getV().size());
    for(int i=0; i<(int)mesh.getI().size(); i++) {
        Vec3i ti=mesh.getI()[i];
        tNeigh[ti[0]].push_back(i);
        tNeigh[ti[1]].push_back(i);
        tNeigh[ti[2]].push_back(i);
    }
    //build new mesh
    vector<Vec3,Eigen::aligned_allocator<Vec3> > vssSeg;
    for(int i=0; i<(int)mesh.getV().size(); i++) {
        const vector<int>& tssn=tNeigh[i];
        boost::unordered_map<int,int> vGroup;
        for(int t=0; t<(int)tssn.size(); t++) {
            int root=(int)dset.find(tssn[t]);
            if(vGroup.find(root) == vGroup.end()) {
                int nrV=(int)vssSeg.size();
                vssSeg.push_back(mesh.getV()[i]);
                findReplace(mesh.getI()[tssn[t]],i,nrV);
                vGroup[root]=nrV++;
            } else findReplace(mesh.getI()[tssn[t]],i,vGroup[root]);
        }
    }
    mesh.getV()=vssSeg;
    mesh.makeUniform();
    mesh.smooth();
}
void FEMMeshFormat::genBeam(const Vec3i& nr,scalar cellSz,const std::string& os)
{
    boost::filesystem::ofstream oss(os);
    oss << "*NODE" << std::endl;

    sizeType vOff=0;
    Vec3i nrG=nr*2+Vec3i::Ones();
    boost::unordered_map<Vec3i,sizeType,Hash> vMap;
    for(sizeType i=0; i<nrG[0]; i++)
        for(sizeType j=0; j<nrG[1]; j++)
            for(sizeType k=0; k<nrG[2]; k++) {
                Vec3 pos=Vec3((scalar)i,(scalar)j,(scalar)k)*cellSz;
                vMap[Vec3i(i,j,k)]=vOff++;
                oss << vOff << ',' << pos[0] << ',' << pos[1] << ',' << pos[2] << std::endl;
            }

    vector<Vec4i,Eigen::aligned_allocator<Vec4i> > tss;
    for(sizeType i=0; i<nr[0]; i++)
        for(sizeType j=0; j<nr[1]; j++)
            for(sizeType k=0; k<nr[2]; k++) {
                Vec3i base(i*2,j*2,k*2);
                //X
                createCell(tss,vMap,base,0,1,2,0);
                createCell(tss,vMap,base,2,1,2,0);
                //Y
                createCell(tss,vMap,base,0,0,2,1);
                createCell(tss,vMap,base,2,0,2,1);
                //Z
                createCell(tss,vMap,base,0,0,1,2);
                createCell(tss,vMap,base,2,0,1,2);
            }

    oss << "*ELEMENT, type=C3D4, ELSET=PART1" << std::endl;
    for(sizeType i=0; i<(sizeType)tss.size(); i++)
        oss << (i+1) << ','
            << (tss[i][0]+1) << ',' << (tss[i][1]+1) << ','
            << (tss[i][2]+1) << ',' << (tss[i][3]+1) << std::endl;
    oss << "*ELSET,ELSET=EALL,GENERATE" << std::endl;
    oss << "1," << nr << std::endl;
}