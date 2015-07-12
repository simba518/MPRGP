#ifndef OBJ_MESH_H
#define OBJ_MESH_H

#include "Config.h"
#include "MathBasic.h"
#include "IO.h"
#include "CollisionDetection.h"

#include <set>
#include <map>
#include <vector>
#include <stack>

PRJ_BEGIN

using namespace std;

template <typename T>
class WriteObjVertex
{
public:
    static void write(char* buf,T& a,T& b,T& c) {
#ifdef _MSC_VER
        sscanf_s(buf,"%f %f %f",&a,&b,&c);
#else
        sscanf(buf,"%f %f %f",&a,&b,&c);
#endif
    }
};
template <>
class WriteObjVertex<scalarD>
{
public:
    static void write(char* buf,scalarD& a,scalarD& b,scalarD& c) {
#ifdef _MSC_VER
        sscanf_s(buf,"%lf %lf %lf",&a,&b,&c);
#else
        sscanf(buf,"%lf %lf %lf",&a,&b,&c);
#endif
    }
};

template <typename T>
class ObjMeshTpl
{
public:
    typedef typename Eigen::Matrix<T,3,1> PT3;
    typedef typename Eigen::Matrix<T,3,3> MAT3;
    struct Edge {
		//subdivision data
		//location of this edge in the 
		//original mesh: 
		//case 1: _subdId.second == -1 
		//	if this edge is inside a triangle, 
		//	in this case _subdId.first is the 
		//	offset in the original _iss list
		//case 2: _subdId.second != -1
		//	if this edge is a subedge of the 
		//	original mesh in this case _subdId 
		//	store the edge-end vertex offset in original vss
		std::pair<int,int> _subdId;
        //offset in the subdivided _vss of 
		//the vertex created for that edge
		int _interNode;

		//info data
        PT3 _nor;
        vector<int> _tris;
        Edge():_nor(0.0f,0.0f,0.0f) {}
    };
    struct EdgeMap {
        friend class ObjMeshTpl;
        struct LSS {
            bool operator()(const pair<int,int>& a,
                            const pair<int,int>& b) const {
                return (a.first < b.first) || (a.first == b.first && a.second < b.second);
            }
        };
        map<pair<int,int>,Edge,LSS> _ess;	//edge
    };
	struct PovTexture{
		virtual int nr() const=0;
		virtual int operator[](int fid) const=0;
		virtual std::string operator()(int tid) const=0;
	};
    ObjMeshTpl():_id(0),_pos(0.0f,0.0f,0.0f),_scale(1.0f),_ctrOff(0.0f,0.0f,0.0f) {
        _trans=MAT3::Identity();
        _dim=3;
    }
    ObjMeshTpl(const int& id):_id(id),_pos(0.0f,0.0f,0.0f),_scale(1.0f),_ctrOff(0.0f,0.0f,0.0f) {
        _trans=MAT3::Identity();
        _dim=3;
    }
    void setDim(int dim){
        _dim=dim;
    }
    virtual ~ObjMeshTpl() {}
	template <typename T2>
	void cast(ObjMeshTpl<T2>& other) const{
		other.getI()=getI();
		other.getV().resize(_vss.size());
		for(sizeType i=0;i<(sizeType)_vss.size();i++)
			other.getV()[i]=_vss[i].cast<T2>();
		other.smooth();
	}
    bool read(istream &is,
              bool move=true,
              bool scale=true,
              int assertFaceType=-1,
              int assertTriangle=-1,
              bool doCheckCCW=true) {
        _dim=3;
        _ctrOff=PT3(0.0f,0.0f,0.0f);

        vector<PT3,Eigen::aligned_allocator<PT3> > tmpVss;
        _vss.swap(tmpVss);

        vector<PT3,Eigen::aligned_allocator<PT3> > tmpNss;
        _nss.swap(tmpNss);

        vector<PT3,Eigen::aligned_allocator<PT3> > tmpFNss;
        _fnss.swap(tmpFNss);

        vector<Vec3i,Eigen::aligned_allocator<Vec3i> > tmpIss;
        _iss.swap(tmpIss);

        vector<int> tmpIssg;
        _issg.swap(tmpIssg);

        map<string,int> tmpGss;
        _gss.swap(tmpGss);

        map<int,string> tmpIgss;
        _igss.swap(tmpIgss);

        char buf[4096],buf2[4096];
        string gName;

        PT3 v;
        Eigen::Matrix<int,3,1> a,b;
        Eigen::Matrix<int,3,1> at,bt;
        Eigen::Matrix<int,3,1> an,bn;
        int gIndex=-1;
        int currIndex=-1;

        while(!is.eof() && is.good()) {
            int faceType=-1;
            int isTriangle=1;

            bool match=false;
            is.getline(buf,std::streamsize(4096),'\n');
            switch(buf[0]) {
            case 'g':
#ifdef _MSC_VER
                ASSERT(sscanf_s(buf,"g %s",buf2,4096) == 1);
#else
                ASSERT(sscanf(buf,"g %s",buf) == 1);
#endif
                gName=buf2;
                if(_gss.find(gName) != _gss.end()) {
                    currIndex=_gss[gName];
                } else {
                    gIndex++;
                    currIndex=gIndex;
                    _gss[gName]=currIndex;
                    _igss[currIndex]=gName;
                }
                break;
            case 'v':
                switch(buf[1]) {
                case ' ':
                    WriteObjVertex<T>::write(buf+2,v.x(),v.y(),v.z());
                    _vss.push_back(v);
                    break;
                case 'n':
                    WriteObjVertex<T>::write(buf+3,v.x(),v.y(),v.z());
                    _nss.push_back(v);
                    break;
                }
                break;
            case 'f':
                //0
                if(!match)
#ifdef _MSC_VER
                    switch(sscanf_s(buf,"f %d %d %d %d",&(a.x()),&(a.y()),&(a.z()),&(b.x())))
#else
                    switch(sscanf(buf,"f %d %d %d %d",&(a.x()),&(a.y()),&(a.z()),&(b.x())))
#endif
                    {
                    case 3:
                        _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        match=true;
                        faceType=0;
                        break;
                    case 4:
                        _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        _iss.push_back(Vec3i(a.x(),a.z(),b.x())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        match=true;
                        faceType=0;
                        isTriangle=0;
                        break;
                    }
                //1
                if(!match)
#ifdef _MSC_VER
                    switch(sscanf_s(buf,"f %d/%d %d/%d %d/%d %d/%d",&(a.x()),&(at.x()), &(a.y()),&(at.y()), &(a.z()),&(at.z()), &(b.x()),&(bt.x())))
#else
                    switch(sscanf(buf,"f %d/%d %d/%d %d/%d %d/%d",&(a.x()),&(at.x()), &(a.y()),&(at.y()), &(a.z()),&(at.z()), &(b.x()),&(bt.x())))
#endif
                    {
                    case 6:
                        _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        match=true;
                        faceType=1;
                        break;
                    case 8:
                        _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        _iss.push_back(Vec3i(a.x(),a.z(),b.x())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        match=true;
                        faceType=1;
                        isTriangle=0;
                        break;
                    }
                //2
                if(!match)
#ifdef _MSC_VER
                    switch(sscanf_s(buf,"f %d//%d %d//%d %d//%d %d//%d",&(a.x()),&(an.x()), &(a.y()),&(an.y()), &(a.z()),&(an.z()), &(b.x()),&(bn.x())))
#else
                    switch(sscanf(buf,"f %d//%d %d//%d %d//%d %d//%d",&(a.x()),&(an.x()), &(a.y()),&(an.y()), &(a.z()),&(an.z()), &(b.x()),&(bn.x())))
#endif
                    {
                    case 6:
                        _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        _fnss.push_back(_nss[an.x()-1]);
                        _fnss.push_back(_nss[an.y()-1]);
                        _fnss.push_back(_nss[an.z()-1]);
                        if(doCheckCCW)
                            checkCCW(_vss[a.x()-1],_nss[an.x()-1],
                                     _vss[a.y()-1],_nss[an.y()-1],
                                     _vss[a.z()-1],_nss[an.z()-1],_iss.back());
                        match=true;
                        faceType=2;
                        break;
                    case 8:
                        _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        _fnss.push_back(_nss[an.x()-1]);
                        _fnss.push_back(_nss[an.y()-1]);
                        _fnss.push_back(_nss[an.z()-1]);
                        if(doCheckCCW)
                            checkCCW(_vss[a.x()-1],_nss[an.x()-1],
                                     _vss[a.y()-1],_nss[an.y()-1],
                                     _vss[a.z()-1],_nss[an.z()-1],_iss.back());

                        _iss.push_back(Vec3i(a.x(),a.z(),b.x())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        _fnss.push_back(_nss[an.x()-1]);
                        _fnss.push_back(_nss[an.z()-1]);
                        _fnss.push_back(_nss[bn.x()-1]);
                        if(doCheckCCW)
                            checkCCW(_vss[a.x()-1],_nss[an.x()-1],
                                     _vss[a.z()-1],_nss[an.z()-1],
                                     _vss[b.x()-1],_nss[bn.x()-1],_iss.back());
                        match=true;
                        faceType=2;
                        isTriangle=0;
                        break;
                    }
                //3
                if(!match)
#ifdef _MSC_VER
                    switch(sscanf_s(buf,"f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d",&(a.x()),&(at.x()),&(an.x()), &(a.y()),&(at.y()),&(an.y()), &(a.z()),&(at.z()),&(an.z()), &(b.x()),&(bt.x()),&(bn.x())))
#else
                    switch(sscanf(buf,"f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d",&(a.x()),&(at.x()),&(an.x()), &(a.y()),&(at.y()),&(an.y()), &(a.z()),&(at.z()),&(an.z()), &(b.x()),&(bt.x()),&(bn.x())))
#endif
                    {
                    case 9:
                        _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        _fnss.push_back(_nss[an.x()-1]);
                        _fnss.push_back(_nss[an.y()-1]);
                        _fnss.push_back(_nss[an.z()-1]);
                        if(doCheckCCW)
                            checkCCW(_vss[a.x()-1],_nss[an.x()-1],
                                     _vss[a.y()-1],_nss[an.y()-1],
                                     _vss[a.z()-1],_nss[an.z()-1],_iss.back());
                        match=true;
                        faceType=3;
                        break;
                    case 12:
                        _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        _fnss.push_back(_nss[an.x()-1]);
                        _fnss.push_back(_nss[an.y()-1]);
                        _fnss.push_back(_nss[an.z()-1]);
                        if(doCheckCCW)
                            checkCCW(_vss[a.x()-1],_nss[an.x()-1],
                                     _vss[a.y()-1],_nss[an.y()-1],
                                     _vss[a.z()-1],_nss[an.z()-1],_iss.back());

                        _iss.push_back(Vec3i(a.x(),a.z(),b.x())-Vec3i(1,1,1));
                        _issg.push_back(currIndex);
                        _fnss.push_back(_nss[an.x()-1]);
                        _fnss.push_back(_nss[an.z()-1]);
                        _fnss.push_back(_nss[bn.x()-1]);
                        if(doCheckCCW)
                            checkCCW(_vss[a.x()-1],_nss[an.x()-1],
                                     _vss[a.z()-1],_nss[an.z()-1],
                                     _vss[b.x()-1],_nss[bn.x()-1],_iss.back());
                        match=true;
                        faceType=3;
                        isTriangle=0;
                        break;
                    }
                //4
                if(!match)
                    return false;

                if(assertFaceType != -1 && assertFaceType != faceType)
                    return false;
                if(assertTriangle != -1 && isTriangle != assertTriangle)
                    return false;
            }
        }

        if(_vss.empty() || _iss.empty()) {
            WARNING("No Vertex In Mesh")
            applyTrans();
            return true;
        }

        for(int i=0,nr=(int)_iss.size(); i<nr; i++) {
            const Vec3i& tri=_iss[i];
            if(tri.x() < (int)0 || tri.x() >= (int)_vss.size()) {
                WARNING("Triangle Index Out Of Bounds")
                return false;
            }
            if(tri.y() < (int)0 || tri.y() >= (int)_vss.size()) {
                WARNING("Triangle Index Out Of Bounds")
                return false;
            }
            if(tri.z() < (int)0 || tri.z() >= (int)_vss.size()) {
                WARNING("Triangle Index Out Of Bounds")
                return false;
            }
        }

        //BBox for vssInput
        BBox<T> box;
        for(int i=0; i<(int)_vss.size(); i++)
            box.setUnion(_vss[i]);

        //if(!(compG(box._minC,PT3(-1E3f,-1E3f,-1E3f)) && compL(box._maxC,PT3(1E3,1E3,1E3))))
        //{
        //	WARNING("Too Big A Mesh")
        //	return false;
        //}

        _trans=MAT3::Identity();
        if(scale) {
            PT3 extent=box.getExtent();
            _scale=1.0f/min(extent.x(),min(extent.y(),extent.z()));
            for(int i=0,sz=(int)_vss.size(); i<sz; i++)
                _vss[i]*=_scale;
        }
        _scale=1.0f;

        _pos=PT3(0.0f,0.0f,0.0f);
        if(move) {
            BBox<T> bb;
            for(int i=0; i<(int)_vss.size(); i++)
                bb.setUnion(_vss[i]);

            for(int i=0,sz=(int)_vss.size(); i<sz; i++)
                _vss[i]-=bb._minC;
        }
		
		if(scale || move)
			applyTrans();
		else if(_nss.size() != _vss.size())
			smooth();
        return true;
    }
    bool readBinary(istream& is) {
        if(!is.good())
            return false;
        if(!readBinaryData(_dim,is).good())
            return false;

        vector<PT3,Eigen::aligned_allocator<PT3> >& vss=getV();
        vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss=getI();
        vector<int>& issg=getIG();
        vector<PT3,Eigen::aligned_allocator<PT3> >& nss=getN();

        int szV;
        int szT;
        int szTG;
        int szG;

        if(!readBinaryData(szV,is).good())
            return false;
        if(!readBinaryData(szT,is).good())
            return false;
        if(!readBinaryData(szTG,is).good())
            return false;
        if(!readBinaryData(szG,is).good())
            return false;

        vss.resize(szV);
        nss.resize(szV);
        iss.resize(szT);
        issg.resize(szTG);

        if(!readVector(vss,is).good())
            return false;
        if(!readVector(nss,is).good())
            return false;
        if(!readVector(iss,is).good())
            return false;
        if(!readVector(issg,is).good())
            return false;

        for(int ig=0; ig<szG; ig++) {
            string name;
            int id;

            readBinaryData(name,is);
            readBinaryData(id,is);
            _gss[name]=id;
            _igss[id]=name;

            if(!is.good())
                return false;
        }

        if(!readBinaryData(_ctrOff,is).good())
            return false;

        return true;
    }
    bool write(ostream &os) const {
        if(!os.good())
            return false;
        if(_dim == 2)
            return false;

        char buf[4096];
        for(int i=0,sz=(int)_vss.size(); i<sz; i++) {
            const PT3& vd=_vss[i];
#ifdef _MSC_VER
            sprintf_s(buf,"v %f %f %f\n",vd.x(),vd.y(),vd.z());
#else
            sprintf(buf,"v %f %f %f\n",vd.x(),vd.y(),vd.z());
#endif
            os << buf;
            if(!os.good())
                return false;
        }

        if(_nss.size() == _vss.size()) {
            for(int i=0,sz=(int)_nss.size(); i<sz; i++) {
                const PT3& vn=_nss[i];
#ifdef _MSC_VER
                sprintf_s(buf,"vn %f %f %f\n",vn.x(),vn.y(),vn.z());
#else
                sprintf(buf,"vn %f %f %f\n",vn.x(),vn.y(),vn.z());
#endif
                os << buf;
                if(!os.good())
                    return false;
            }
        }

        if(_issg.empty()) {
            for(int i=0,sz=(int)_iss.size(); i<sz; i++) {
#ifdef _MSC_VER
                if(_nss.size() == _vss.size())
                    sprintf_s(buf,"f %d//%d %d//%d %d//%d\n",
                              _iss[i].x()+1,_iss[i].x()+1,
                              _iss[i].y()+1,_iss[i].y()+1,
                              _iss[i].z()+1,_iss[i].z()+1);
                else
                    sprintf_s(buf,"f %d %d %d\n",
                              _iss[i].x()+1,
                              _iss[i].y()+1,
                              _iss[i].z()+1);
#else
                if(_nss.size() == _vss.size())
                    sprintf(buf,"f %d//%d %d//%d %d//%d\n",
                            _iss[i].x()+1,_iss[i].x()+1,
                            _iss[i].y()+1,_iss[i].y()+1,
                            _iss[i].z()+1,_iss[i].z()+1);
                else
                    sprintf(buf,"f %d %d %d\n",
                            _iss[i].x()+1,
                            _iss[i].y()+1,
                            _iss[i].z()+1);
#endif
                os << buf;
            }
        } else {
            map<int,string> groupOs;
            for(int i=0,sz=(int)_iss.size(); i<sz; i++) {
                if(groupOs.find(_issg[i]) == groupOs.end()) {
                    groupOs[_issg[i]]=string();
                    if(_issg[i] != -1) {
                        string currGroupName=_igss.find(_issg[i])->second;
                        string& str=groupOs[_issg[i]];
                        str+="g ";
                        str+=currGroupName;
                        str+="\n";
                    }
                }

#ifdef _MSC_VER
                if(_nss.size() == _vss.size())
                    sprintf_s(buf,"f %d//%d %d//%d %d//%d\n",
                              _iss[i].x()+1,_iss[i].x()+1,
                              _iss[i].y()+1,_iss[i].y()+1,
                              _iss[i].z()+1,_iss[i].z()+1);
                else
                    sprintf_s(buf,"f %d %d %d\n",
                              _iss[i].x()+1,
                              _iss[i].y()+1,
                              _iss[i].z()+1);
#else
                if(_nss.size() == _vss.size())
                    sprintf(buf,"f %d//%d %d//%d %d//%d\n",
                            _iss[i].x()+1,_iss[i].x()+1,
                            _iss[i].y()+1,_iss[i].y()+1,
                            _iss[i].z()+1,_iss[i].z()+1);
                else
                    sprintf(buf,"f %d %d %d\n",
                            _iss[i].x()+1,
                            _iss[i].y()+1,
                            _iss[i].z()+1);
#endif
                groupOs[_issg[i]]+=buf;
            }

            map<int,string>::iterator
            begin=groupOs.begin(),end=groupOs.end();
            for(; begin!=end; begin++)
                os << begin->second;
        }

        return os.good();
    }
    bool write(const boost::filesystem::path& path) const {
        boost::filesystem::ofstream ofs(path);
        return write(ofs);
    }
	bool writePov(const boost::filesystem::path& path,bool normal=false,const PovTexture* tex=NULL) const
	{
		return writePov(boost::filesystem::ofstream(path),normal,tex);
	}
	bool writePov(ostream &os,bool normal=false,const PovTexture* tex=NULL) const {
        if(_dim == 2)
            return false;
#define CACHE_SIZE 1024

        os << "mesh2 {\n";

        //vertices
        {
            os << "    vertex_vectors {" << _vss.size() << "\n";
            int nr=(int)_vss.size();
            ostringstream oss;

            for(int i=0; i<(int)_vss.size(); i++) {
                char buf[CACHE_SIZE];
                const PT3& v=_vss[i];
#ifdef _MSC_VER
                if(i == nr-1)
                    sprintf_s<CACHE_SIZE>(buf,"        <%f,%f,%f>\n",v.x(),v.y(),v.z());
                else
                    sprintf_s<CACHE_SIZE>(buf,"        <%f,%f,%f>,\n",v.x(),v.y(),v.z());
#else
                if(i == nr-1)
                    sprintf(buf,"        <%f,%f,%f>\n",v.x(),v.y(),v.z());
                else
                    sprintf(buf,"        <%f,%f,%f>,\n",v.x(),v.y(),v.z());
#endif
                oss << buf;
            }
            os << oss.str() << "    }\n";
        }

        //normals
        if(normal && _nss.size() == _vss.size()) {
            os << "    normal_vectors {" << _nss.size() << "\n";
            int nr=(int)_nss.size();
            ostringstream oss;

            for(int i=0; i<nr; i++) {
                char buf[CACHE_SIZE];
                const PT3& n=_nss[i];
#ifdef _MSC_VER
                if(i == nr-1)
                    sprintf_s<CACHE_SIZE>(buf,"        <%f,%f,%f>\n",n.x(),n.y(),n.z());
                else
                    sprintf_s<CACHE_SIZE>(buf,"        <%f,%f,%f>,\n",n.x(),n.y(),n.z());
#else
                if(i == nr-1)
                    sprintf(buf,"        <%f,%f,%f>\n",n.x(),n.y(),n.z());
                else
                    sprintf(buf,"        <%f,%f,%f>,\n",n.x(),n.y(),n.z());
#endif
                oss << buf;
            }

            os << oss.str() << "    }\n";
        }

        //face indices
		if(tex){
			int nrT=(*tex).nr();
			os << "    texture_list {" << nrT << "\n";
			for(int i=0;i<nrT;i++)
			{
				os << (*tex)(i);
				if(i<nrT-1)
					os << ",";
				os << "\n";
			}
            os << "    }\n";
		}
        {
            os << "    face_indices {" << _iss.size() << "\n";
            int nr=(int)_iss.size();
            ostringstream oss;

            for(int i=0; i<nr; i++) {
                char buf[CACHE_SIZE];
                const Vec3i& f=_iss[i];
				if(tex){
#ifdef _MSC_VER
                if(i == nr-1)
                    sprintf_s<CACHE_SIZE>(buf,"        <%d,%d,%d>%d\n",f.x(),f.y(),f.z(),(*tex)[i]);
                else
                    sprintf_s<CACHE_SIZE>(buf,"        <%d,%d,%d>%d,\n",f.x(),f.y(),f.z(),(*tex)[i]);
#else
                if(i == nr-1)
                    sprintf(buf,"        <%d,%d,%d>%d\n",f.x(),f.y(),f.z(),(*tex)[i]);
                else
                    sprintf(buf,"        <%d,%d,%d>%d,\n",f.x(),f.y(),f.z(),(*tex)[i]);
#endif
				}else{
#ifdef _MSC_VER
                if(i == nr-1)
                    sprintf_s<CACHE_SIZE>(buf,"        <%d,%d,%d>\n",f.x(),f.y(),f.z());
                else
                    sprintf_s<CACHE_SIZE>(buf,"        <%d,%d,%d>,\n",f.x(),f.y(),f.z());
#else
                if(i == nr-1)
                    sprintf(buf,"        <%d,%d,%d>\n",f.x(),f.y(),f.z());
                else
                    sprintf(buf,"        <%d,%d,%d>,\n",f.x(),f.y(),f.z());
#endif
				}
                oss << buf;
            }
            os << oss.str() << "    }\n";
        }

        if(normal && _nss.size() == _vss.size()) {
            os << "    normal_indices {" << _iss.size() << "\n";
            int nr=(int)_iss.size();
            ostringstream oss;

            for(int i=0; i<nr; i++) {
                char buf[CACHE_SIZE];
                const Vec3i& f=_iss[i];
#ifdef _MSC_VER
                if(i == nr-1)
                    sprintf_s<CACHE_SIZE>(buf,"        <%d,%d,%d>\n",f.x(),f.y(),f.z());
                else
                    sprintf_s<CACHE_SIZE>(buf,"        <%d,%d,%d>,\n",f.x(),f.y(),f.z());
#else
                if(i == nr-1)
                    sprintf(buf,"        <%d,%d,%d>\n",f.x(),f.y(),f.z());
                else
                    sprintf(buf,"        <%d,%d,%d>,\n",f.x(),f.y(),f.z());
#endif
                oss << buf;
            }

            os << oss.str() << "    }\n";
        }
        os << "}\n";
        return os.good();

#undef CACHE_SIZE
    }
    bool writePBRT(ostream &os) const {
        if(!os.good())
            return false;
        if(_dim == 2)
            return false;

        os << "AttributeBegin" << endl;
        os << "Shape \"trianglemesh\"" << endl;
        if(!os.good())
            return false;

        os << "\"point P\" ["<< endl;
        for(int j=0; j<_vss.size(); j++)
            os << _vss[j].x() << " " << _vss[j].y() << " " << _vss[j].z() << endl;
        os << "]" << endl;
        if(!os.good())
            return false;

        os << "\"integer indices\" [" << endl;
        for(int j=0; j<_iss.size(); j++)
            os << _iss[j].x() << " " << _iss[j].y() << " " << _iss[j].z() << endl;
        os << "]" << endl;
        os << "AttributeEnd" << endl;
        if(!os.good())
            return false;

        return true;
    }
    bool writeVTK(const boost::filesystem::path& path,bool binary,bool normal=false,bool vertexNormal=false) const {
        VTKWriter<T> os("WaveFront Obj Mesh",path,binary);
        if(_dim == 3)
            return writeVTK3D(os,normal,vertexNormal);
        else return writeVTK2D(os,normal,vertexNormal);
    }
    bool writeVTK(VTKWriter<T>& os,bool normal=false,bool vertexNormal=false) const{
        if(_dim == 3)
            return writeVTK3D(os,normal,vertexNormal);
        else return writeVTK2D(os,normal,vertexNormal);
    }
    bool writeVTK3D(VTKWriter<T>& os,bool normal=false,bool vertexNormal=false) const{
        {
            os.setRelativeIndex();
            os.appendPoints(getV().begin(),getV().end());
            os.appendCells(getI().begin(),getI().end(),VTKWriter<T>::TRIANGLE,true);
        }
        if(_tnss.size() == _iss.size() && normal){
            os.setRelativeIndex();
            std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
            for(int i=0;i<(int)_tnss.size();i++){
                PT3 ctr=getTC((int)i);
                vss.push_back(ctr);
                vss.push_back(ctr+_tnss[i]);
            }
            os.appendPoints(vss.begin(),vss.end());
            // original is
            /*
              VTKWriter<T>::IteratorIndex<Vec3i> beg(0,2,0),end(vss.size()/2,2,0);
             */
            typename VTKWriter<T>::template IteratorIndex<Vec3i> beg(0,2,0);
            typename VTKWriter<T>::template IteratorIndex<Vec3i> end(vss.size()/2,2,0);
            os.appendCells(beg,end,VTKWriter<T>::LINE,true);
        }
		if(_nss.size() == _vss.size() && vertexNormal){
			os.setRelativeIndex();
            std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
            for(int i=0;i<(int)_vss.size();i++){
                vss.push_back(_vss[i]);
                vss.push_back(_vss[i]+_nss[i]);
            }
            os.appendPoints(vss.begin(),vss.end());
            typename VTKWriter<T>::template IteratorIndex<Vec3i> beg(0,2,0);
            typename VTKWriter<T>::template IteratorIndex<Vec3i> end(vss.size()/2,2,0);
            os.appendCells(beg,end,VTKWriter<T>::LINE,true);
		}
        return true;
    }
    bool writeVTK2D(VTKWriter<T>& os,bool normal=false,bool vertexNormal=false) const{
        {
            os.setRelativeIndex();
            os.appendPoints(getV().begin(),getV().end());
            os.appendCells(getI().begin(),getI().end(),VTKWriter<T>::LINE,true);
        }
        if(normal){
            os.setRelativeIndex();
            std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
            for(int i=0;i<(int)_iss.size();i++){
                Vec3i LI=_iss[i];
                PT3 ctr=(_vss[LI[0]]+_vss[LI[1]])/2.0f;
                vss.push_back(ctr);
                vss.push_back(ctr+_tnss[i]);
            }
            os.appendPoints(vss.begin(),vss.end());
            typename VTKWriter<T>::template IteratorIndex<Vec3i> beg(0,2,0),end(vss.size()/2,2,0);
            os.appendCells(beg,end,VTKWriter<T>::LINE,true);
        }
		if(_nss.size() == _vss.size() && vertexNormal){
			os.setRelativeIndex();
            std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
            for(int i=0;i<(int)_vss.size();i++){
                vss.push_back(_vss[i]);
                vss.push_back(_vss[i]+_nss[i]);
            }
            os.appendPoints(vss.begin(),vss.end());
            typename VTKWriter<T>::template IteratorIndex<Vec3i> beg(0,2,0);
            typename VTKWriter<T>::template IteratorIndex<Vec3i> end(vss.size()/2,2,0);
            os.appendCells(beg,end,VTKWriter<T>::LINE,true);
		}
        return true;
    }
    bool writeBinary(ostream &os) const {
        if(!os.good())
            return false;
        if(!writeBinaryData(_dim,os).good())
            return false;

        const vector<PT3,Eigen::aligned_allocator<PT3> >& vss=getV();
        const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss=getI();
        const vector<int>& issg=getIG();
        const vector<PT3,Eigen::aligned_allocator<PT3> >& nss=getN();
        const map<string,int>& gss=_gss;

        const int szV=(int)vss.size();
        const int szT=(int)iss.size();
        const int szTG=(int)issg.size();
        const int szG=(int)gss.size();

        if(!writeBinaryData(szV,os).good())
            return false;
        if(!writeBinaryData(szT,os).good())
            return false;
        if(!writeBinaryData(szTG,os).good())
            return false;
        if(!writeBinaryData(szG,os).good())
            return false;

        if(!writeVector(vss,os).good())
            return false;
        if(!writeVector(nss,os).good())
            return false;
        if(!writeVector(iss,os).good())
            return false;
        if(!writeVector(issg,os).good())
            return false;
        for(map<string,int>::const_iterator
                beg=gss.begin(),end=gss.end(); beg != end; beg++) {
            writeBinaryData(beg->first,os);
            writeBinaryData(beg->second,os);
            if(!os.good())
                return false;
        }

        if(!writeBinaryData(_ctrOff,os).good())
            return false;

        return true;
    }
    bool write(ostream &vs,ostream &fs,int &index) const {
        if(!vs.good() || !fs.good())
            return false;
        if(_dim == 2)
            return false;

        int start=index;
        WARNING("Group Info In ObjMesh Is Ignored!")

        char buf[4096];
        for(int i=0,sz=(int)_vss.size(); i<sz; i++) {
            const PT3& vd=_vss[i];
#ifdef _MSC_VER
            sprintf_s(buf,"v %f %f %f\n",vd.x(),vd.y(),vd.z());
#else
            sprintf(buf,"v %f %f %f\n",vd.x(),vd.y(),vd.z());
#endif
            index++;
            vs << buf;
            if(!vs.good())
                return false;
        }

        for(int i=0,sz=(int)_iss.size(); i<sz; i++) {
#ifdef _MSC_VER
            sprintf_s(buf,"f %d %d %d\n", start+_iss[i].x(), start+_iss[i].y(), start+_iss[i].z());
#else
            sprintf(buf,"f %d %d %d\n", start+_iss[i].x(), start+_iss[i].y(), start+_iss[i].z());
#endif
            fs << buf;
            if(!fs.good())
                return false;
        }

        return true;
    }
    void writeCsv(const boost::filesystem::path& path) const
	{
#define WRITE_POINT(ID)	os << _vss[ID][0] << "," << _vss[ID][1] << "," << _vss[ID][2] << "," << std::endl;
		boost::filesystem::ofstream os(path);
		int nrT=(int)_iss.size();
		os << nrT << "," << std::endl;
		for(int i=0;i<nrT;i++)
		{
			WRITE_POINT(_iss[i][0])
			WRITE_POINT(_iss[i][1])
			if(getDim() == 3)
				WRITE_POINT(_iss[i][2])
		}
	}
	void addMesh(const ObjMeshTpl& mesh,const std::string& g)
	{
		int vOff=(int)_vss.size();
		_vss.insert(_vss.end(),mesh._vss.begin(),mesh._vss.end());
		for(int i=0;i<(int)mesh._iss.size();i++){
			_iss.push_back(mesh._iss[i]+Vec3i::Constant(vOff));
			_issg.push_back((int)_gss.size());
		}
		ASSERT(_gss.find(g) == _gss.end())
		_gss[g]=(int)_gss.size();
	}
    PT3 getTC(int i) const {
        if(_dim == 3)
            return (_vss[_iss[i].x()]+_vss[_iss[i].y()]+_vss[_iss[i].z()])/3.0f;
        else return (_vss[_iss[i].x()]+_vss[_iss[i].y()])/2.0f;
    }
    const PT3& getV(int i) const {
        return _vss[i];
    }
    const Vec3i& getI(int i) const {
        return _iss[i];
    }
    const int& getIG(int i) const {
        return _issg[i];
    }
    const PT3& getN(int i) const {
        return _nss[i];
    }
    const PT3& getTN(int i) const {
        return _tnss[i];
    }
    map<string,int>& getGS() {
        return _gss;
    }
    map<int,string>& getIGS() {
        return _igss;
    }
    PT3& getV(int i) {
        return _vss[i];
    }
    Vec3i& getI(int i) {
        return _iss[i];
    }
    int& getIG(int i) {
        return _issg[i];
    }
    PT3& getN(int i) {
        return _nss[i];
    }
    PT3& getTN(int i) {
        return _tnss[i];
    }
    T getArea(int i) const {
        if(_dim == 3)
            return TriangleTpl<T>(_vss[_iss[i][0]],_vss[_iss[i][1]],_vss[_iss[i][2]]).area();
        else return LineSegTpl<T>(_vss[_iss[i][0]],_vss[_iss[i][1]]).length();
    }
    const vector<PT3,Eigen::aligned_allocator<PT3> >& getV() const {
        return _vss;
    }
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& getI() const {
        return _iss;
    }
    const vector<int>& getIG() const {
        return _issg;
    }
    const vector<PT3,Eigen::aligned_allocator<PT3> >& getN() const {
        return _nss;
    }
    const vector<PT3,Eigen::aligned_allocator<PT3> >& getFN() const {
        return _fnss;
    }
    const vector<PT3,Eigen::aligned_allocator<PT3> >& getTN() const {
        return _tnss;
    }
    const map<string,int>& getGS() const {
        return _gss;
    }
    const map<int,string>& getIGS() const {
        return _igss;
    }
    string getTG(int i) const {
        if(_iss[i].w() == -1)
            return string();
        else _igss.find(_iss[i].w())->second;
    }
    int getGId(const string& name) const {
        if(_gss.find(name) == _gss.end())
            return -1;
        else return _gss.find(name)->second;
    }
    vector<PT3,Eigen::aligned_allocator<PT3> >& getV() {
        return _vss;
    }
    vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& getI() {
        return _iss;
    }
    vector<int>& getIG() {
        return _issg;
    }
    vector<PT3,Eigen::aligned_allocator<PT3> >& getN() {
        return _nss;
    }
    vector<PT3,Eigen::aligned_allocator<PT3> >& getFN() {
        return _fnss;
    }
    vector<PT3,Eigen::aligned_allocator<PT3> >& getTN() {
        return _tnss;
    }
    const MAT3& getT() const {
        return _trans;
    }
    MAT3& getT() {
        return _trans;
    }
    const PT3& getPos() const {
        return _pos;
    }
    PT3& getPos() {
        return _pos;
    }
    const T& getScale() const {
        return _scale;
    }
    T& getScale() {
        return _scale;
    }
    void applyTrans() {
        BBox<T> bb;
        for(int i=0; i<(int)_vss.size(); i++)
            bb.setUnion(_vss[i]);
        PT3 ctr=(bb._minC+bb._maxC)*0.5f;
        applyTrans(ctr);
    }
    void applyTrans(const PT3& customCtr) {
        for(int i=0,sz=(int)_vss.size(); i<sz; i++)
            _vss[i]=_trans*(_vss[i]-customCtr)*_scale+customCtr+_pos;
        _ctrOff=_trans*(_ctrOff)*_scale;

        _trans=MAT3::Identity();
        _pos=PT3(0.0f,0.0f,0.0f);
        _scale=1.0f;

        smooth();
    }
    BBox<T> getBB() const {
        BBox<T> ret;
        for(int i=0; i<(int)_vss.size(); i++)
            ret.setUnion(_vss[i]);
        return ret;
    }
    const int& getId() const {
        return _id;
    }
    //physics properties
    T getVolume() const {
        T volume = 0.0;
        for(int it=0,nr=(int)_iss.size(); it<nr; it++) {
            if(_dim == 3){
                const PT3 &p1 = _vss[_iss[it].x()];
                const PT3 &p2 = _vss[_iss[it].y()];
                const PT3 &p3 = _vss[_iss[it].z()];
                volume += TriangleTpl<T>(p1, p2, p3).signedVolume();
            }else{
                const PT3 &p1 = _vss[_iss[it].x()];
                const PT3 &p2 = _vss[_iss[it].y()];
                volume += TriangleTpl<T>(PT3::Zero(), p1, p2).area();
            }
        }
        return std::abs(volume);
    }
    PT3 getCentroid() const {
        //must be convex
        PT3 numerator(0.0f,0.0f,0.0f);
        T denominator(0.0f);
        for(int i=0,nr=(int)_iss.size(); i<nr; i++) {
            T area=getArea(i);
            PT3 ctr=getTC(i);
            numerator+=ctr*area;
            denominator+=area;
        }
        return _ctrOff+numerator/denominator;
    }
    PT3 getVolumeCentroid() const {
        //must be convex
        PT3 numerator(0.0f,0.0f,0.0f);
        T denominator(0.0f);
        for(int i=0,nr=(int)_iss.size(); i<nr; i++) {
            Vec3i TI=_iss[i];
            T vol;
            if(_dim == 3){
                TetrahedronTpl<T> tet(PT3::Zero(),_vss[TI[0]],_vss[TI[1]],_vss[TI[2]]);
                vol=tet.volume();
                if(tet._swap)
                    vol*=-1.0f;
                numerator+=tet.masscenter()*vol;
            }else{
                TriangleTpl<T> tri(PT3::Zero(),_vss[TI[0]],_vss[TI[1]]);
                vol=tri.area();
                numerator+=tri.masscenter()*vol;
            }
            denominator+=vol;
        }
        return _ctrOff+numerator/denominator;
    }
    PT3& centroidOffset() {
        return _ctrOff;
    }
    const PT3& centroidOffset() const {
        return _ctrOff;
    }
    //simple utility
    T getMass(const T& dens) const {
        return getVolume()*dens;
    }
    void subdivide(const int& nrIter,vector<std::pair<int,int> >* vssInfo=NULL)
	{
		EdgeMap eMap;
        buildEdge(eMap);
		for(typename map<pair<int,int>,Edge,typename EdgeMap::LSS>::iterator
            beg=eMap._ess.begin(),end=eMap._ess.end(); beg != end; beg++)
			beg->second._subdId=beg->first;
		if(vssInfo)
			vssInfo->assign(_vss.size(),std::pair<int,int>(-1,-1));
		for(int it=0,mod=1;it<nrIter;it++,mod*=4) {
            if(_dim == 3)
				subdivideSingle3D(eMap,mod,vssInfo);
			else subdivideSingle2D();
        }
        smooth();
    }
    void subdivideSingle3D(EdgeMap& eMap,int mod,vector<std::pair<int,int> >* vssInfo=NULL) {
		EdgeMap eMapOut;

        //first for each edge insert an internode
        for(typename map<pair<int,int>,Edge,typename EdgeMap::LSS>::iterator
            beg=eMap._ess.begin(),end=eMap._ess.end(); beg != end; beg++) {
            //insert node
            beg->second._interNode=(int)_vss.size();
            _vss.push_back((_vss[beg->first.first]+_vss[beg->first.second])*0.5f);
			if(vssInfo)
				vssInfo->push_back(beg->second._subdId);
			//parent info
			{
				std::pair<int,int> PA(beg->second._interNode,beg->first.first);
				std::pair<int,int> PB(beg->second._interNode,beg->first.second);
				if(PA.first > PA.second)std::swap(PA.first,PA.second);
				if(PB.first > PB.second)std::swap(PB.first,PB.second);
				Edge EA;EA._subdId=beg->second._subdId;
				Edge EB;EB._subdId=beg->second._subdId;
				eMapOut._ess[PA]=EA;
				eMapOut._ess[PB]=EB;
			}
		}

        //second replace each triangle with four
        vector<Vec3i,Eigen::aligned_allocator<Vec3i> > issNew;
        vector<int> issgNew;
        for(int i=0,nr=(int)_iss.size(); i<nr; i++) {
            const Vec3i& iT=_iss[i];
            
            const int v0=(int)iT.x();
            const int v1=(int)iT.y();
            const int v2=(int)iT.z();
            const int v3=getE(v0,v1,eMap)._interNode;
            const int v4=getE(v1,v2,eMap)._interNode;
            const int v5=getE(v2,v0,eMap)._interNode;

            issNew.push_back(Vec3i(v0,v3,v5));
            issNew.push_back(Vec3i(v3,v1,v4));
            issNew.push_back(Vec3i(v5,v3,v4));
            issNew.push_back(Vec3i(v5,v4,v2));
            
            //make same group
            if(_issg.size() == _iss.size())
            {
                const int& igT=_issg[i];
                issgNew.push_back(igT);
                issgNew.push_back(igT);
                issgNew.push_back(igT);
                issgNew.push_back(igT);
            }

			//parent info
			{
				std::pair<int,int> PA(v3,v4);
				std::pair<int,int> PB(v3,v5);
				std::pair<int,int> PC(v4,v5);
				if(PA.first > PA.second)std::swap(PA.first,PA.second);
				if(PB.first > PB.second)std::swap(PB.first,PB.second);
				if(PC.first > PC.second)std::swap(PC.first,PC.second);
				Edge EA;EA._subdId=std::pair<int,int>(i/mod,-1);
				Edge EB;EB._subdId=std::pair<int,int>(i/mod,-1);
				Edge EC;EC._subdId=std::pair<int,int>(i/mod,-1);
				eMapOut._ess[PA]=EA;
				eMapOut._ess[PB]=EB;
				eMapOut._ess[PC]=EC;
			}
        }
		eMap=eMapOut;

        _iss=issNew;
        if(issgNew.size() == _iss.size())
            _issg=issgNew;
    }
    void subdivideSingle2D() {
        vector<Vec3i,Eigen::aligned_allocator<Vec3i> > issNew;
        for(int i=0;i<(int)_iss.size();i++){
            Vec3i LI=_iss[i];
            int back=(int)_vss.size();
            
            issNew.push_back(Vec3i(LI[0],back,0));
            issNew.push_back(Vec3i(back,LI[1],0));
            
            _vss.push_back((_vss[LI[0]]+_vss[LI[1]])/2.0f);
        }

        _iss=issNew;
    }
	void marchingCube(const std::vector<T>& MCVal,ObjMeshTpl<T>& mesh)
	{
		//marching cube on surface
        mesh=ObjMeshTpl();
		//calculate normal dot view dir
		std::vector<int> vid(_vss.size(),-1);
		for(int i=0;i<(sizeType)_vss.size();i++)
		{
			if(MCVal[i] < 0.0f)
			{
				vid[i]=(int)mesh.getV().size();
				mesh.getV().push_back(_vss[i]);
			}
		}
		//generate vertex
		EdgeMap eMap;
		buildEdge(eMap);
        for(auto
			beg=eMap._ess.begin(),end=eMap._ess.end();beg!=end;beg++)
		{
			std::pair<int,int> k=beg->first;
			Edge& e=beg->second;
			e._interNode=-1;
			if((MCVal[k.first] <  0.0f && MCVal[k.second] >= 0.0f) ||
			   (MCVal[k.first] >= 0.0f && MCVal[k.second] <  0.0f))
			{
				e._interNode=(int)mesh.getV().size();
				mesh.getV().push_back((_vss[k.second]*MCVal[k.first]-_vss[k.first]*MCVal[k.second])/
									  (MCVal[k.first]-MCVal[k.second]));
			}
		}
		//generate index
		int nrT[8]={1,2,2,1,2,1,1,0};
		int idT[8][6]=
		{
			{0,1,2, -1,-1,-1},
			{1,5,3, 1,2,5},
			{0,3,4, 0,4,2},
			{2,5,4, -1,-1,-1},

			{0,4,5, 0,1,4},
			{1,4,3, -1,-1,-1},
			{0,3,5, -1,-1,-1},
			{-1,-1,-1, -1,-1,-1},
		};
		for(int i=0;i<(int)_iss.size();i++)
		{
			const Vec3i& I=_iss[i];
			unsigned char type=0;
			if(MCVal[I[0]] >= 0.0f)type+=1;
			if(MCVal[I[1]] >= 0.0f)type+=2;
			if(MCVal[I[2]] >= 0.0f)type+=4;
            for(int t=0;t<nrT[type];t++)
			{
				Vec3i IBK;
				for(int V=0;V<3;V++)
				{
                    int VBK=idT[type][t*3+V];
					if(VBK < 3)
						IBK[V]=vid[I[VBK]];
					else{
						std::pair<int,int> e((int)I[VBK-3],(int)I[(VBK-2)%3]);
						if(e.first > e.second)std::swap(e.first,e.second);
						IBK[V]=eMap._ess[e]._interNode;
					}
				}
				ASSERT(compGE(IBK,Vec3i::Zero()));
				mesh.getI().push_back(IBK);
			}
		}
		mesh.smooth();
	}
    //smooth
    void makeUnique()
    {
        std::vector<bool> valid;
        valid.assign(_vss.size(),false);
        for(int i=0;i<int(_iss.size());i++)
        {
            valid[_iss[i][0]]=true;
            valid[_iss[i][1]]=true;
            if(getDim() == 3)
                valid[_iss[i][2]]=true;
        }

        std::map<sizeType,sizeType> remap;
        sizeType index=0;
        for(int i=0;i<int(_vss.size());i++)
            if(valid[i])
            {
                _vss[index]=_vss[i];
                remap[i]=index++;
            }
        _vss.resize(index);

        for(int i=0;i<int(_iss.size());i++)
        {
            _iss[i][0]=remap[_iss[i][0]];
            _iss[i][1]=remap[_iss[i][1]];
            if(getDim() == 3)
                _iss[i][2]=remap[_iss[i][2]];
        }
    }
	void makeUniform()
	{
		//initialize edge
		EdgeMap eMap;
		buildEdge(eMap);
		//make uniform
		vector<bool> visited(_iss.size(),false);
		for(int i=0;i<(int)visited.size();i++)
		if(!visited[i])
		{
			std::stack<int> queue;
			queue.push(i);
			visited[i]=true;
			while(!queue.empty())
			{
				int ti=queue.top();queue.pop();
				const Vec3i& I=_iss[ti];
				std::pair<int,int> e;
				for(int eid=0;eid<3;eid++)
				{
					e.first=(int)I[eid];
					e.second=(int)I[(eid+1)%3];
					if(e.first > e.second)std::swap(e.first,e.second);
					const Edge& edg=eMap._ess[e];
					for(int tid=0;tid<(int)edg._tris.size();tid++)
					{
						int tj=edg._tris[tid];
						if(tj != ti && !visited[tj])
						{
							makeUniform(ti,tj,e.first,e.second);
							queue.push(tj);
							visited[tj]=true;
						}
					}
				}
			}
		}
	}
	void makeUniform(int i,int j,int v0,int v1)
	{
		int v0i,v1i,v0j,v1j;
		for(int d=0;d<3;d++)
		{
			if(_iss[i][d] == v0)v0i=d;
			if(_iss[i][d] == v1)v1i=d;
			if(_iss[j][d] == v0)v0j=d;
			if(_iss[j][d] == v1)v1j=d;
		}
		bool isI=(v0i+1)%3 == v1i;
		bool isJ=(v0j+1)%3 == v1j;
		if(isI == isJ)std::swap(_iss[j][1],_iss[j][2]);
	}
    void smooth() {
		T feps=numeric_limits<T>::min()*1E10f;
        //force copy and restrict
        vector<PT3,Eigen::aligned_allocator<PT3> > tmpNss;
        _nss.swap(tmpNss);

        vector<PT3,Eigen::aligned_allocator<PT3> > tmpTnss;
        _tnss.swap(tmpTnss);

        _nss.resize(_vss.size());
        _tnss.resize(_iss.size());
        vector<int> nr_face(_vss.size(),0);

        for(int k=0,sz=(int)_vss.size(); k<sz; k++)
            _nss[k].setConstant(0.0f);

        //iterate face
        if(_dim == 2)
        {
            for(int k=0,sz=(int)_iss.size(); k<sz; k++) 
            {
                const Vec3i &f=_iss[k];

                const PT3 e0=_vss[f.y()]-_vss[f.x()];
                PT3 n(e0.y(),-e0.x(),0.0f);
                if(n.norm() > feps)
                    n.normalize();
                else
                    n=PT3::Zero();

                _tnss[k]=n;
                nr_face[f.x()]++;
                _nss[f.x()]+=n;
                nr_face[f.y()]++;
                _nss[f.y()]+=n;
            }
        }
        else
        {
            for(int k=0,sz=(int)_iss.size(); k<sz; k++) 
            {
                const Vec3i &f=_iss[k];

                const PT3 e0=_vss[f.y()]-_vss[f.x()];
                const PT3 e1=_vss[f.z()]-_vss[f.x()];
                const PT3 e2=_vss[f.z()]-_vss[f.y()];

                PT3 n=(e0).cross(e1);
                if(n.norm() > feps)
                    n.normalize();
                else
                    n=PT3::Zero();

                const T a0=getAngle3D<T>(_vss[f.y()]-_vss[f.x()],_vss[f.z()]-_vss[f.x()]);
                const T a1=getAngle3D<T>(_vss[f.z()]-_vss[f.y()],_vss[f.x()]-_vss[f.y()]);
                const T a2=getAngle3D<T>(_vss[f.x()]-_vss[f.z()],_vss[f.y()]-_vss[f.z()]);

                _tnss[k]=n;
                nr_face[f.x()]++;
                _nss[f.x()]+=n*a0;
                nr_face[f.y()]++;
                _nss[f.y()]+=n*a1;
                nr_face[f.z()]++;
                _nss[f.z()]+=n*a2;
            }
        }

        //compute normal
        for(int k=0,sz=(int)_nss.size(); k<sz; k++) {
            if(_nss[k].norm() > feps)
                _nss[k].normalize();
            else
                _nss[k]=PT3::Zero();
        }
    }
	//topology
    const Edge& getE(int a,int b,const EdgeMap& eMap) const {
        ASSERT(!eMap._ess.empty())
        if(a<b) return eMap._ess.find(pair<int,int>(a,b))->second;
        else return eMap._ess.find(pair<int,int>(b,a))->second;
    }
    void buildEdge(EdgeMap& eMap) const {
        //edge
        eMap._ess.clear();
        for(int k=0,sz=(int)_iss.size(); k<sz; k++) {
            const Vec3i &f=_iss[k];
            addEdge((int)f.x(),(int)f.y(),k,eMap);
            addEdge((int)f.y(),(int)f.z(),k,eMap);
            addEdge((int)f.x(),(int)f.z(),k,eMap);
        }

        //edge normal
        for(typename map<pair<int,int>,Edge,typename EdgeMap::LSS>::iterator
                begin=eMap._ess.begin(),end=eMap._ess.end(); begin != end; begin++) {
            begin->second._nor.normalize();
        }
    }
	void buildKRingV(vector<map<int,int> >& KRing,int r) const
	{
		EdgeMap eMap;
		buildEdge(eMap);
		int nrV=(int)_vss.size();
		//initialize 1-ring triangle
		vector<set<int> > oneRingV;
		oneRingV.assign(nrV,set<int>());
		KRing.assign(nrV,map<int,int>());
		for(typename map<pair<int,int>,Edge,typename EdgeMap::LSS>::const_iterator
			beg=eMap._ess.begin(),end=eMap._ess.end();beg!=end;beg++)
		{
			oneRingV[beg->first.first].insert(beg->first.second);
			oneRingV[beg->first.second].insert(beg->first.first);
			findInsertV(KRing[beg->first.first],beg->first.second,0);
			findInsertV(KRing[beg->first.second],beg->first.first,0);
		}
		//compute r ring triangle
		for(int i=1;i<r;i++)
		{
			vector<map<int,int> > KRingOld=KRing;
			for(int v=0;v<nrV;v++)
			{
				const map<int,int>& Ring=KRingOld[v];
				for(map<int,int>::const_iterator 
					beg=Ring.begin(),end=Ring.end();
					beg!=end;beg++)
				{
					findInsert(KRing[v],oneRingV[beg->first],i);
				}
			}
		}
	}
	void buildKRing(vector<map<int,int> >& KRing,int r) const
	{
		EdgeMap eMap;
		buildEdge(eMap);
		int nrV=(int)_vss.size();
		//initialize 1-ring triangle
		vector<set<int> > oneRingT;
		KRing.assign(nrV,map<int,int>());
		oneRingT.assign(nrV,set<int>());
		for(typename map<pair<int,int>,Edge,typename EdgeMap::LSS>::const_iterator
			beg=eMap._ess.begin(),end=eMap._ess.end();beg!=end;beg++)
		{
			oneRingT[beg->first.first].insert(beg->second._tris.begin(),beg->second._tris.end());
			oneRingT[beg->first.second].insert(beg->second._tris.begin(),beg->second._tris.end());
			findInsert(KRing[beg->first.first],beg->second._tris,0);
			findInsert(KRing[beg->first.second],beg->second._tris,0);
		}
		//compute r ring triangle
		for(int i=1;i<r;i++)
		{
			vector<map<int,int> > KRingOld=KRing;
			for(int v=0;v<nrV;v++)
			{
				const map<int,int>& Ring=KRingOld[v];
				for(map<int,int>::const_iterator 
					beg=Ring.begin(),end=Ring.end();
					beg!=end;beg++)
				{
					const Vec3i& t=_iss[beg->first];
					findInsert(KRing[v],oneRingT[t[0]],i);
					findInsert(KRing[v],oneRingT[t[1]],i);
					findInsert(KRing[v],oneRingT[t[2]],i);
				}
			}
		}
	}
	template <typename TC>
	void findInsert(map<int,int>& Ring,const TC& tris,int currR) const
	{
		for(typename TC::const_iterator beg=tris.begin(),end=tris.end();beg!=end;beg++)
			if(Ring.find(*beg) == Ring.end())
				Ring.insert(std::pair<int,int>(*beg,currR));
	}
	void findInsertV(map<int,int>& Ring,const int& v,int currR) const
	{
		if(Ring.find(v) == Ring.end())
			Ring.insert(std::pair<int,int>(v,currR));
	}
    //dimension
    const int& getDim() const{return _dim;}
    int& getDim(){return _dim;}
protected:
    void addEdge(int a,int b,int tri,EdgeMap& eMap) const {
        if(a>b)swap(a,b);
        Edge& edg=eMap._ess[pair<int,int>(a,b)];
        edg._nor+=_tnss[tri];
        edg._tris.push_back(tri);
    }
    void checkCCW(const PT3& v1,const PT3& n1,
                  const PT3& v2,const PT3& n2,
                  const PT3& v3,const PT3& n3,Vec3i& i) {
        PT3 n=(n1+n2+n3)/3.0f;
        if((v2-v1).cross(v3-v1).dot(n) < 0)
            swap(i.y(),i.z());
    }
protected:
    vector<PT3,Eigen::aligned_allocator<PT3> > _vss;	//vert
    vector<PT3,Eigen::aligned_allocator<PT3> > _nss;	//normal
    vector<PT3,Eigen::aligned_allocator<PT3> > _fnss;	//tree normals per face
    vector<PT3,Eigen::aligned_allocator<PT3> > _tnss;	//tri_normal
    vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _iss;//index
    vector<int> _issg;
    map<string,int> _gss;//groups
    map<int,string> _igss;//inverse groups
    MAT3 _trans;
    int _id;
    PT3 _pos;
    T _scale;
    PT3 _ctrOff;
    int _dim;
};

typedef ObjMeshTpl<scalarD> ObjMeshD;
typedef ObjMeshTpl<scalarF> ObjMeshF;
typedef ObjMeshTpl<scalar> ObjMesh;

PRJ_END

#endif
