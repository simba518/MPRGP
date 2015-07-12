#ifndef IMPLICIT_FUNC_UTILS_H
#define IMPLICIT_FUNC_UTILS_H

#include "ImplicitFunc.h"

PRJ_BEGIN

class ImplicitFuncTranspose : public ImplicitFunc<scalar>
{
public:
    ImplicitFuncTranspose(boost::shared_ptr<ImplicitFunc<scalar> > inner,const Vec3& rot,const Vec3& pos) {
        _inner=inner;
        _rot=makeRotation<scalar>(rot);
        _pos=pos;
    }
    scalar operator()(const Vec3& pos) const {
        return (*_inner)(_rot.transpose()*(pos-_pos));
    }
    boost::shared_ptr<ImplicitFunc<scalar> > _inner;
    Mat3 _rot;
    Vec3 _pos;
};
class ImplicitFuncTransposeRef : public ImplicitFunc<scalar>
{
public:
    ImplicitFuncTransposeRef(const ImplicitFunc<scalar>& inner,const Vec3& rot,const Vec3& pos) 
        :_inner(inner)
    {
        _rot=makeRotation<scalar>(rot);
        _pos=pos;
    }
    scalar operator()(const Vec3& pos) const {
        return _inner(_rot.transpose()*(pos-_pos));
    }
    const ImplicitFunc<scalar>& _inner;
    Mat3 _rot;
    Vec3 _pos;
};
template <typename T>
class ConvertToImplicitFuncCSG
{
public:
    static boost::shared_ptr<ImplicitFunc<scalar> > convert(const T& IF) {
        return IF;
    }
};
template <>
class ConvertToImplicitFuncCSG<BBox<scalar,3> >
{
public:
    static boost::shared_ptr<ImplicitFunc<scalar> > convert(const BBox<scalar,3>& IF) {
        return boost::shared_ptr<ImplicitFuncCSG>(new ImplicitFuncCSG(IF,ImplicitFuncCSG::INTERSECT));
    }
};
template <>
class ConvertToImplicitFuncCSG<std::string>
{
public:
    class ImplicitFuncHeightField : public ImplicitFuncMesh3D
    {
    public:
        ImplicitFuncHeightField(const ObjMesh& mesh)
            :ImplicitFuncMesh3D(mesh,1.0f,mesh.getBB()) {}
        scalar operator()(const Vec3& pos) const {
            scalarSLC d;
            if(_bvh->distanceToRay3D(Vec3SLC(pos.x(),pos.y(),pos.z()),Vec3SLC(pos.x(),pos.y(),pos.z()+1.0f),0,d,true))
                return -(scalar)d;
            else if(_bvh->distanceToRay3D(Vec3SLC(pos.x(),pos.y(),pos.z()),Vec3SLC(pos.x(),pos.y(),pos.z()-1.0f),0,d,true))
                return (scalar)d;
            return _bb.distTo(pos);
        }
    };
    static boost::shared_ptr<ImplicitFunc<scalar> > convert(const std::string& path) {
        ObjMesh mesh;
        boost::filesystem::ifstream is(path);
        mesh.read(is,false,false);
        return boost::shared_ptr<ImplicitFuncHeightField>(new ImplicitFuncHeightField(mesh));
    }
};
class ImplicitFuncCSGUtils
{
public:
    //version with all implicit func
    static boost::shared_ptr<ImplicitFuncCSG> getUnion(boost::shared_ptr<ImplicitFunc<scalar> > a,boost::shared_ptr<ImplicitFunc<scalar> > b) {
        boost::shared_ptr<ImplicitFuncCSG> ret(new ImplicitFuncCSG(ImplicitFuncCSG::UNION));
        ret->_a=a;
        ret->_b=b;
        return ret;
    }
    static boost::shared_ptr<ImplicitFuncCSG> getIntersect(boost::shared_ptr<ImplicitFunc<scalar> > a,boost::shared_ptr<ImplicitFunc<scalar> > b) {
        boost::shared_ptr<ImplicitFuncCSG> ret(new ImplicitFuncCSG(ImplicitFuncCSG::INTERSECT));
        ret->_a=a;
        ret->_b=b;
        return ret;
    }
    static boost::shared_ptr<ImplicitFuncCSG> getSubtract(boost::shared_ptr<ImplicitFunc<scalar> > a,boost::shared_ptr<ImplicitFunc<scalar> > b) {
        boost::shared_ptr<ImplicitFuncCSG> ret(new ImplicitFuncCSG(ImplicitFuncCSG::SUBTRACT));
        ret->_a=a;
        ret->_b=b;
        return ret;
    }
    //version with other type
    template <typename IFA,typename IFB>
    static boost::shared_ptr<ImplicitFuncCSG> getUnion(const IFA& a,const IFB& b) {
        return getUnion(ConvertToImplicitFuncCSG<IFA>::convert(a),ConvertToImplicitFuncCSG<IFB>::convert(b));
    }
    template <typename IFA,typename IFB>
    static boost::shared_ptr<ImplicitFuncCSG> getIntersect(const IFA& a,const IFB& b) {
        return getIntersect(ConvertToImplicitFuncCSG<IFA>::convert(a),ConvertToImplicitFuncCSG<IFB>::convert(b));
    }
    template <typename IFA,typename IFB>
    static boost::shared_ptr<ImplicitFuncCSG> getSubtract(const IFA& a,const IFB& b) {
        return getSubtract(ConvertToImplicitFuncCSG<IFA>::convert(a),ConvertToImplicitFuncCSG<IFB>::convert(b));
    }
    //reinitialize
    static boost::shared_ptr<ImplicitFuncReinit> getReinit(const ScalarField& s,const boost::shared_ptr<ImplicitFuncCSG>& b) {
        return boost::shared_ptr<ImplicitFuncReinit>(new ImplicitFuncReinit(s,*b));
    }
    static boost::shared_ptr<ImplicitFuncReinit> reinit(ScalarField& s,const boost::shared_ptr<ImplicitFuncCSG>& b) {
        boost::shared_ptr<ImplicitFuncReinit> tmp(new ImplicitFuncReinit(s,*b));
        GridOp<scalar,scalar>::copyFromOtherGridOfSameGeometry(s,tmp->_ls);
        return tmp;
    }
    //negate
    static boost::shared_ptr<ImplicitFuncNegate> negate(boost::shared_ptr<ImplicitFunc<scalar> > b) {
        boost::shared_ptr<ImplicitFuncNegate> ret(new ImplicitFuncNegate());
        ret->_inner=b;
        return ret;
    }
    //traspose
    template <typename IF>
    static boost::shared_ptr<ImplicitFuncTranspose> transpose(const IF& a,const Vec3& rot,const Vec3& pos) {
        boost::shared_ptr<ImplicitFuncTranspose> ret(new ImplicitFuncTranspose(ConvertToImplicitFuncCSG<IF>::convert(a),rot,pos));
        return ret;
    }
    //add box to objMesh
    static void addBox(ObjMesh& mesh,const BBox<scalar>& bb) {
        std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& vss=mesh.getV();
        std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& nss=mesh.getN();
        std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss=mesh.getI();

        //Z
        {
            iss.push_back(Vec3i(vss.size(),vss.size()+1,vss.size()+2));
            iss.push_back(Vec3i(vss.size(),vss.size()+2,vss.size()+3));

            nss.push_back(Vec3(0.0f,0.0f,-1.0f));
            nss.push_back(Vec3(0.0f,0.0f,-1.0f));
            nss.push_back(Vec3(0.0f,0.0f,-1.0f));
            nss.push_back(Vec3(0.0f,0.0f,-1.0f));

            vss.push_back(Vec3(bb._minC.x(),bb._minC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._minC.x(),bb._maxC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._maxC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._minC.y(),bb._minC.z()));
        }

        {
            iss.push_back(Vec3i(vss.size(),vss.size()+1,vss.size()+2));
            iss.push_back(Vec3i(vss.size(),vss.size()+2,vss.size()+3));

            nss.push_back(Vec3(0.0f,0.0f,1.0f));
            nss.push_back(Vec3(0.0f,0.0f,1.0f));
            nss.push_back(Vec3(0.0f,0.0f,1.0f));
            nss.push_back(Vec3(0.0f,0.0f,1.0f));

            vss.push_back(Vec3(bb._minC.x(),bb._minC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._minC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._maxC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._minC.x(),bb._maxC.y(),bb._maxC.z()));
        }

        //X
        {
            iss.push_back(Vec3i(vss.size(),vss.size()+1,vss.size()+2));
            iss.push_back(Vec3i(vss.size(),vss.size()+2,vss.size()+3));

            nss.push_back(Vec3(-1.0f,0.0f,0.0f));
            nss.push_back(Vec3(-1.0f,0.0f,0.0f));
            nss.push_back(Vec3(-1.0f,0.0f,0.0f));
            nss.push_back(Vec3(-1.0f,0.0f,0.0f));

            vss.push_back(Vec3(bb._minC.x(),bb._minC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._minC.x(),bb._minC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._minC.x(),bb._maxC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._minC.x(),bb._maxC.y(),bb._minC.z()));
        }

        {
            iss.push_back(Vec3i(vss.size(),vss.size()+1,vss.size()+2));
            iss.push_back(Vec3i(vss.size(),vss.size()+2,vss.size()+3));

            nss.push_back(Vec3(1.0f,0.0f,0.0f));
            nss.push_back(Vec3(1.0f,0.0f,0.0f));
            nss.push_back(Vec3(1.0f,0.0f,0.0f));
            nss.push_back(Vec3(1.0f,0.0f,0.0f));

            vss.push_back(Vec3(bb._maxC.x(),bb._minC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._maxC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._maxC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._minC.y(),bb._maxC.z()));
        }

        //Y
        {
            iss.push_back(Vec3i(vss.size(),vss.size()+1,vss.size()+2));
            iss.push_back(Vec3i(vss.size(),vss.size()+2,vss.size()+3));

            nss.push_back(Vec3(0.0f,-1.0f,0.0f));
            nss.push_back(Vec3(0.0f,-1.0f,0.0f));
            nss.push_back(Vec3(0.0f,-1.0f,0.0f));
            nss.push_back(Vec3(0.0f,-1.0f,0.0f));

            vss.push_back(Vec3(bb._minC.x(),bb._minC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._minC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._minC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._minC.x(),bb._minC.y(),bb._maxC.z()));
        }

        {
            iss.push_back(Vec3i(vss.size(),vss.size()+1,vss.size()+2));
            iss.push_back(Vec3i(vss.size(),vss.size()+2,vss.size()+3));

            nss.push_back(Vec3(0.0f,1.0f,0.0f));
            nss.push_back(Vec3(0.0f,1.0f,0.0f));
            nss.push_back(Vec3(0.0f,1.0f,0.0f));
            nss.push_back(Vec3(0.0f,1.0f,0.0f));

            vss.push_back(Vec3(bb._minC.x(),bb._maxC.y(),bb._minC.z()));
            vss.push_back(Vec3(bb._minC.x(),bb._maxC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._maxC.y(),bb._maxC.z()));
            vss.push_back(Vec3(bb._maxC.x(),bb._maxC.y(),bb._minC.z()));
        }
    }
    //get gradient
    static Vec3 sampleGrad(const ImplicitFunc<scalar>& f,const Vec3& pt,const Vec3& cellSz)
    {
        Vec3 ret=Vec3::Zero();
        if(cellSz.x() > EPS)
            ret.x()=(f(pt+Vec3(cellSz.x()*0.5f,0.0f,0.0f))-f(pt-Vec3(cellSz.x()*0.5f,0.0f,0.0f)))/cellSz.x();
        if(cellSz.y() > EPS)
            ret.y()=(f(pt+Vec3(0.0f,cellSz.y()*0.5f,0.0f))-f(pt-Vec3(0.0f,cellSz.y()*0.5f,0.0f)))/cellSz.y();
        if(cellSz.z() > EPS)
            ret.z()=(f(pt+Vec3(0.0f,0.0f,cellSz.z()*0.5f))-f(pt-Vec3(0.0f,0.0f,cellSz.z()*0.5f)))/cellSz.z();
        return ret;
    }
};

PRJ_END

#endif
