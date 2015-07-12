#ifndef COLLISION_DETECTION_H
#define COLLISION_DETECTION_H

#include "MathBasic.h"
#include "IO.h"

PRJ_BEGIN

template <typename T>
class SAT
{
public:
    typedef typename Eigen::Matrix<T,3,1> PT;
    typedef typename Eigen::Matrix<T,2,1> PT2;
public:
    static bool testSegment(const PT2& a,const PT2& b) {
        return a.x() > b.y() || b.x() > a.y();
    }
    template <typename T1,typename T2>
    static bool sep(const PT& d,const T1& A,const T2& B) {
        return testSegment(A.project(d),B.project(d));
    }
};
template <typename T>
class LineSegTpl
{
public:
    typedef typename Eigen::Matrix<T,3,1> PT;
public:
    LineSegTpl(const PT& x,const PT& y):_x(x),_y(y) {}
    T length() const {
        return (_y-_x).norm();
    }
    PT circumcenter() const {
        return (_x+_y)*0.5f;
    }
    PT masscenter() const {
        return (_x+_y)*0.5f;
    }
    PT normal() const{
        PT ret=_y-_x;
        return PT(ret.y(),-ret.x(),0.0f)/
               std::max(ret.norm(),ScalarUtil<T>::scalar_eps);
    }
    PT gradDir(sizeType id) const{
        PT dir=_x-_y;
        dir/=dir.squaredNorm();
        if(id == 1)
            dir*=-1.0f;
        return dir;
    }
    void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT& b) const {
        b[1]=(pt-_x).dot(_y-_x)/std::max<T>((_y-_x).squaredNorm(),ScalarUtil<T>::scalar_eps);
        b[1]=std::max<T>(0.0f,std::min<T>(1.0f,b[1]));
        b[0]=1.0f-b[1];
        cp=_x*b[0]+_y*b[1];
        sqrDistance=(cp-pt).squaredNorm();
    }
    void writeVTK(VTKWriter<T>& os) const
    {
        os.setRelativeIndex();
        std::vector<PT,Eigen::aligned_allocator<PT> > vss;
        vss.push_back(_x);
        vss.push_back(_y);
        os.appendPoints(vss.begin(),vss.end());

        std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
        iss.push_back(Vec3i(0,1,0));
        os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
    }
public:
    //data
    PT _x;
    PT _y;
};
template <typename T>
class PlaneTpl
{
public:
    typedef typename Eigen::Matrix<T,3,1> PT;
    typedef typename Eigen::Matrix<T,2,1> PT2;
    PlaneTpl() {}
    PlaneTpl(const PT& x0,const PT& n)
        :_x0(x0),_n(n) {}
    PlaneTpl(const PT& a,const PT& b,const PT& c)
        :_x0(a),_n((b-a).cross(c-a)) {}
    //intersection
    T side(const PT& p) const {
        return (p-_x0).dot(_n);
    }
    bool intersect(const BBox<T>& bb) const {
        return SAT<T>::testSegment(PT2::Constant(_x0.dot(_n)),bb.project(_n));
    }
    PT2 project(const PT& d) const {
        return (d==_n) ? PT2::Constant(_x0.dot(_n)) :
               PT2(-ScalarUtil<T>::scalar_max,ScalarUtil<T>::scalar_max);
    }
    void writeVTK(VTKWriter<T>& os) const
    {
#define SEG 128
        {
            sizeType mid;
            _n.maxCoeff(&mid);
            PT nx=(mid == 0 || mid == 1) ? PT(-_n[1],_n[0],0.0f) : PT(-_n[2],0.0f,_n[0]);nx.normalize();
            PT ny=_n.cross(nx).normalized();
            os.setRelativeIndex();
            std::vector<PT,Eigen::aligned_allocator<PT> > vss;
            std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
            for(sizeType i=0;i<SEG;i++)
            {
                scalar ANG=2.0f*M_PI*(scalar)i/(scalar)SEG;
                vss.push_back(nx*cos(ANG)+ny*sin(ANG));
                iss.push_back(Vec3i(i,(i+1)%SEG,0));
            }
            os.appendPoints(vss.begin(),vss.end());
            os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
        }
#undef SEG

        {
            os.setRelativeIndex();
            std::vector<PT,Eigen::aligned_allocator<PT> > vss;
            std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
            vss.push_back(_x0);
            vss.push_back(_n+_x0);
            iss.push_back(Vec3i(0,1,0));
            os.appendPoints(vss.begin(),vss.end());
            os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
        }
    }
public:
    //data
    PT _x0;
    PT _n;
};
template <typename T>
class TriangleTpl
{
public:
    typedef typename Eigen::Matrix<T,3,1> PT;
    typedef typename Eigen::Matrix<T,2,1> PT2;
    typedef typename Eigen::Matrix<T,2,2> MAT2;
    typedef typename Eigen::Matrix<T,3,3> MAT3;
public:
    TriangleTpl() {}
    TriangleTpl(const PT& a,const PT& b,const PT& c):_a(a),_b(b),_c(c) {}
    PT circumcenter() const {
        PT alpha=_b-_a;
        PT beta=_c-_a;

        PT2 ls(0.5f,0.0f);
        PT2 le(-alpha.dot(beta)/alpha.dot(alpha),1.0f);

        PT2 rs(0.0f,0.5f);
        PT2 re(1.0f,-alpha.dot(beta)/beta.dot(beta));

        MAT2 A;
        A.col(0)=le;
        A.col(1)=-re;
        T t=(A.inverse()*(rs-ls))[0];
        ls+=le*t;

        return _a+alpha*ls.x()+beta*ls.y();
    }
    PT masscenter() const{return (_a+_b+_c)/3.0f;}
    PT bary(const PT& pt) const {
        Eigen::Matrix<T,3,2> A;
        A.col(0)=_a-_c;
        A.col(1)=_b-_c;
        MAT2 ATA=A.transpose()*A;
        PT2 ab=ATA.inverse()*(A.transpose()*(pt-_c));
        return PT(ab.x(),ab.y(),1.0f-ab.sum());
    }
    T area() const {
        return (_b-_a).cross(_c-_a).norm()*0.5f;
    }
    T signedVolume() const {
        T v321 = _c.x() * _b.y() * _a.z();
        T v231 = _b.x() * _c.y() * _a.z();
        T v312 = _c.x() * _a.y() * _b.z();
        T v132 = _a.x() * _c.y() * _b.z();
        T v213 = _b.x() * _a.y() * _c.z();
        T v123 = _a.x() * _b.y() * _c.z();
        return (1.0f/6.0f)*(-v321+v231+v312-v132-v213+v123);
    }
    PT normal() const {
        PT dir=(_b-_a).cross(_c-_a);
        return dir/std::max(dir.norm(),ScalarUtil<T>::scalar_eps);
    }
    PT gradDir(sizeType id) const{
        switch(id){
        case 0: return height(_a,_b,_c);
        case 1: return height(_b,_a,_c);
        default:return height(_c,_a,_b);}
    }
    PT height(const PT& a,const PT& b,const PT& c) const{
        PT ret=a-b;
        PT bc=(c-b).normalized();
        ret=ret-ret.dot(bc)*bc;
        return ret/ret.squaredNorm();
    }
    //intersect
    bool isInside(const PT& pt) const {
        PT bc=bary(pt);
        PT ptPlane=bc.x()*_a+bc.y()*_b+bc.z()*_c;
        return (ptPlane-pt).norm() < ScalarUtil<T>::scalar_eps &&
               bc.x() >= 0.0f && bc.y() >= 0.0f && bc.z() >= 0.0f;
    }
    bool isPlanePointInside(const PT& pt) const {
        PT bc=bary(pt);
        return bc.x() >= 0.0f && bc.y() >= 0.0f && bc.z() >= 0.0f;
    }
    bool updateIntr(T& s,T& t,T num,T denom) const {
        if(std::abs(denom) < EPS)
            return num <= 0.0f;
        else if(denom > 0.0f)
            s=max<T>(s,num/denom);
        else
            t=min<T>(t,num/denom);
        return s<t;
    }
    bool intersect(const LineSegTpl<T>& l,T& s,T& t) const {
        //compute A,B
        MAT2 m;
        m(0,0)=(_a-_c).x();
        m(0,1)=(_b-_c).x();
        m(1,0)=(_a-_c).y();
        m(1,1)=(_b-_c).y();
        PT2 A=m.inverse()*PT2(l._x.x()-  _c.x(),l._x.y()-  _c.y());
        PT2 B=m.inverse()*PT2(l._y.x()-l._x.x(),l._y.y()-l._x.y());

        //get the expected interval
        s=0.0f,t=1.0f;
        return updateIntr(s,t,-A.x(),B.x()) &&
               updateIntr(s,t,-A.y(),B.y()) &&
               updateIntr(s,t,A.x()+A.y()-1.0f,-(B.x()+B.y()));
    }
    bool intersect(const LineSegTpl<T>& l,T& t,bool infinite=false) const {
        PT abt;
        MAT3 mat;
        mat.col(0)=_a-_c;
        mat.col(1)=_b-_c;
        mat.col(2)=l._x-l._y;
        if(std::abs(mat.determinant()) < ScalarUtil<T>::scalar_eps)
            return false;

        abt=mat.inverse()*(l._x-_c);
        t=abt.z();
        if(infinite)
            abt.z()=min<T>(abt.z(),0.5f);
        return (abt.x() >= 0.0f && abt.y() >= 0.0f && (abt.x()+abt.y()) <= 1.0f) &&	//in triangle
               (abt.z() >= 0.0f && abt.z() <= 1.0f);	//in segment
    }
    bool intersect(const BBox<T>& bb) const {
        PT d=(_b-_a).cross(_c-_a);
        if(SAT<T>::testSegment(bb.project(d),PT2::Constant(_a.dot(d))))
            return false;
        for(sizeType i=0; i<3; i++)
            if(SAT<T>::sep(PT::Unit(i),*this,bb))
                return false;
        for(sizeType i=0; i<3; i++)
            if(SAT<T>::sep(PT::Unit(i).cross(d),*this,bb))
                return false;
        return true;
    }
	bool calcLineDist(const LineSegTpl<T>& l,PT& bt,PT2& bl) const{
		bt=bary(l._x);
		PT dbt=bary(l._y)-bt;
		
		T b=0.0f,t=1.0f;
		for(char i=0;i<3;i++)
		{
			if(bt[i] < 0.0f && dbt[i] < 0.0f)
				return false;
			if(bt[i] < 0.0f && dbt[i] > 0.0f)
				b=std::max(b,(-bt[i])/dbt[i]);
			if(bt[i] > 0.0f && dbt[i] < 0.0f)
				t=std::min(t,bt[i]/(-dbt[i]));
		}
		if(b>t)return false;
		bt+=dbt*b;
		bl=PT2(1.0f-b,b);
		return true;
	}
	template <typename PT_BARY>
    void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT_BARY& b) const {
        PT diff = _a - pt;
        PT edge0 = _b - _a;
        PT edge1 = _c - _a;
        T a00 = edge0.dot(edge0);
        T a01 = edge0.dot(edge1);
        T a11 = edge1.dot(edge1);
        T b0 = diff.dot(edge0);
        T b1 = diff.dot(edge1);
        T c = diff.dot(diff);
        T det = std::abs(a00*a11 - a01*a01);
        T s = a01*b1 - a11*b0;
        T t = a01*b0 - a00*b1;

        if (s + t <= det) {
            if (s < (T)0) {
                if (t < (T)0) { // region 4
                    if (b0 < (T)0) {
                        t = (T)0;
                        if (-b0 >= a00) {
                            s = (T)1;
                            sqrDistance = a00 + ((T)2)*b0 + c;
                        } else {
                            s = -b0/a00;
                            sqrDistance = b0*s + c;
                        }
                    } else {
                        s = (T)0;
                        if (b1 >= (T)0) {
                            t = (T)0;
                            sqrDistance = c;
                        } else if (-b1 >= a11) {
                            t = (T)1;
                            sqrDistance = a11 + ((T)2)*b1 + c;
                        } else {
                            t = -b1/a11;
                            sqrDistance = b1*t + c;
                        }
                    }
                } else { // region 3
                    s = (T)0;
                    if (b1 >= (T)0) {
                        t = (T)0;
                        sqrDistance = c;
                    } else if (-b1 >= a11) {
                        t = (T)1;
                        sqrDistance = a11 + ((T)2)*b1 + c;
                    } else {
                        t = -b1/a11;
                        sqrDistance = b1*t + c;
                    }
                }
            } else if (t < (T)0) { // region 5
                t = (T)0;
                if (b0 >= (T)0) {
                    s = (T)0;
                    sqrDistance = c;
                } else if (-b0 >= a00) {
                    s = (T)1;
                    sqrDistance = a00 + ((T)2)*b0 + c;
                } else {
                    s = -b0/a00;
                    sqrDistance = b0*s + c;
                }
            } else { // region 0
                // minimum at interior point
                T invDet = ((T)1)/det;
                s *= invDet;
                t *= invDet;
                sqrDistance = s*(a00*s + a01*t + ((T)2)*b0) +
                              t*(a01*s + a11*t + ((T)2)*b1) + c;
            }
        } else {
            T tmp0, tmp1, numer, denom;

            if (s < (T)0) { // region 2
                tmp0 = a01 + b0;
                tmp1 = a11 + b1;
                if (tmp1 > tmp0) {
                    numer = tmp1 - tmp0;
                    denom = a00 - ((T)2)*a01 + a11;
                    if (numer >= denom) {
                        s = (T)1;
                        t = (T)0;
                        sqrDistance = a00 + ((T)2)*b0 + c;
                    } else {
                        s = numer/denom;
                        t = (T)1 - s;
                        sqrDistance = s*(a00*s + a01*t + ((T)2)*b0) +
                                      t*(a01*s + a11*t + ((T)2)*b1) + c;
                    }
                } else {
                    s = (T)0;
                    if (tmp1 <= (T)0) {
                        t = (T)1;
                        sqrDistance = a11 + ((T)2)*b1 + c;
                    } else if (b1 >= (T)0) {
                        t = (T)0;
                        sqrDistance = c;
                    } else {
                        t = -b1/a11;
                        sqrDistance = b1*t + c;
                    }
                }
            } else if (t < (T)0) { // region 6
                tmp0 = a01 + b1;
                tmp1 = a00 + b0;
                if (tmp1 > tmp0) {
                    numer = tmp1 - tmp0;
                    denom = a00 - ((T)2)*a01 + a11;
                    if (numer >= denom) {
                        t = (T)1;
                        s = (T)0;
                        sqrDistance = a11 + ((T)2)*b1 + c;
                    } else {
                        t = numer/denom;
                        s = (T)1 - t;
                        sqrDistance = s*(a00*s + a01*t + ((T)2)*b0) +
                                      t*(a01*s + a11*t + ((T)2)*b1) + c;
                    }
                } else {
                    t = (T)0;
                    if (tmp1 <= (T)0) {
                        s = (T)1;
                        sqrDistance = a00 + ((T)2)*b0 + c;
                    } else if (b0 >= (T)0) {
                        s = (T)0;
                        sqrDistance = c;
                    } else {
                        s = -b0/a00;
                        sqrDistance = b0*s + c;
                    }
                }
            } else { // region 1
                numer = a11 + b1 - a01 - b0;
                if (numer <= (T)0) {
                    s = (T)0;
                    t = (T)1;
                    sqrDistance = a11 + ((T)2)*b1 + c;
                } else {
                    denom = a00 - ((T)2)*a01 + a11;
                    if (numer >= denom) {
                        s = (T)1;
                        t = (T)0;
                        sqrDistance = a00 + ((T)2)*b0 + c;
                    } else {
                        s = numer/denom;
                        t = (T)1 - s;
                        sqrDistance = s*(a00*s + a01*t + ((T)2)*b0) +
                                      t*(a01*s + a11*t + ((T)2)*b1) + c;
                    }
                }
            }
        }

        // Account for numerical round-off error.
        if (sqrDistance < (T)0) {
            sqrDistance = (T)0;
        }

        cp = _a + s*edge0 + t*edge1;
        b(1) = s;
        b(2) = t;
        b(0) = (T)1 - s - t;
    }
    PT2 project(const PT& d) const {
        T a=_a.dot(d);
        T b=_b.dot(d);
        T c=_c.dot(d);
        return PT2(min<T>(a,min<T>(b,c)),max<T>(a,max<T>(b,c)));
    }
    void writeVTK(VTKWriter<T>& os) const
    {
        os.setRelativeIndex();
        std::vector<PT,Eigen::aligned_allocator<PT> > vss;
        vss.push_back(_a);
        vss.push_back(_b);
        vss.push_back(_c);
        os.appendPoints(vss.begin(),vss.end());

        std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
        iss.push_back(Vec3i(0,1,0));
        iss.push_back(Vec3i(1,2,0));
        iss.push_back(Vec3i(2,0,0));
        os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
    }
public:
    //data
    PT _a;
    PT _b;
    PT _c;
};
template <typename T>
class TetrahedronTpl
{
public:
    typedef typename Eigen::Matrix<T,2,1> PT2;
    typedef typename Eigen::Matrix<T,3,1> PT;
    typedef typename Eigen::Matrix<T,4,1> PT4;
    typedef typename Eigen::Matrix<T,6,1> PT6;
    typedef typename Eigen::Matrix<T,3,3> MAT3;
public:
    TetrahedronTpl() {}
    TetrahedronTpl(const PT& a,const PT& b,const PT& c,const PT& d)
        :_a(a),_b(b),_c(c),_d(d) {
        _swap=false;
        if(volume() < 0.0f) {
            swap(_c,_d);
            _swap=true;
        }
    }
    PT circumcenter() const {
        MAT3 A;
        A.row(0)=_b-_a;
        A.row(1)=_c-_a;
        A.row(2)=_d-_a;
        PT B=PT(A.row(0).squaredNorm(),
                A.row(1).squaredNorm(),
                A.row(2).squaredNorm())*0.5f;
        MAT3 AI=A.inverse();
        return AI*B+_a;
    }
    PT masscenter() const{return (_a+_b+_c+_d)/4.0f;}
    PT4 bary(const PT& pt) const {
        MAT3 A;
        A.col(0)=_a-_d;
        A.col(1)=_b-_d;
        A.col(2)=_c-_d;
        PT abc=A.inverse()*(pt-_d);
        return PT4(abc.x(),abc.y(),abc.z(),1.0f-abc.sum());
    }
    bool isInside(const PT& pt) const {
        PT4 bc=bary(pt);
        return bc.x() >= 0 && bc.y() >= 0 && bc.z() >= 0 && bc.w() >= 0;
    }
    T volume() const {
        return (_b-_a).cross(_c-_a).dot(_d-_a)/6.0f;
    }
    bool dualCellVolume(PT4& dv) {
        dv=PT4::Zero();
        PT cc=circumcenter();
        bool safe=isInside(cc);
        safe=safe&&dualFaceVolume(dv[0],dv[1],dv[2],_a,_b,_c,cc);
        safe=safe&&dualFaceVolume(dv[0],dv[2],dv[3],_a,_c,_d,cc);
        safe=safe&&dualFaceVolume(dv[0],dv[3],dv[1],_a,_d,_b,cc);
        safe=safe&&dualFaceVolume(dv[1],dv[2],dv[3],_b,_c,_d,cc);
        return safe;
    }
    static bool dualFaceVolume(T& va,T& vb,T& vc,const PT& a,const PT& b,const PT& c,const PT& cc) {
        TriangleTpl<T> tri(a,b,c);
        PT ccf=tri.circumcenter();
        va+=TetrahedronTpl(cc,ccf,a,(a+b)*0.5f).volume()+TetrahedronTpl(cc,ccf,a,(a+c)*0.5f).volume();
        vb+=TetrahedronTpl(cc,ccf,b,(b+c)*0.5f).volume()+TetrahedronTpl(cc,ccf,b,(b+a)*0.5f).volume();
        vc+=TetrahedronTpl(cc,ccf,c,(c+a)*0.5f).volume()+TetrahedronTpl(cc,ccf,c,(c+b)*0.5f).volume();
        return tri.isPlanePointInside(ccf);
    }
    static T dihedralAngle(const PT& a,const PT& b,const PT& c,const PT& d) {
        PT n1=(a-c).cross(b-c).normalized();
        PT n2=(a-d).cross(b-d).normalized();
        return M_PI*0.5f-asin(n1.dot(n2));
    }
    PT6 dihedralAngleTet() {
        PT6 ret;
        ret[0]=dihedralAngle(_a,_b,_c,_d);
        ret[1]=dihedralAngle(_a,_c,_b,_d);
        ret[2]=dihedralAngle(_a,_d,_b,_c);
        ret[3]=dihedralAngle(_b,_c,_a,_d);
        ret[4]=dihedralAngle(_b,_d,_a,_c);
        ret[5]=dihedralAngle(_c,_d,_a,_b);
        return ret;
    }
    //intersection
	bool calcLineDist(const LineSegTpl<T>& l,PT4& bt,PT2& bl) const{
		bt=bary(l._x);
		PT4 dbt=bary(l._y)-bt;
		
		T b=0.0f,t=1.0f;
		for(char i=0;i<4;i++)
		{
			if(bt[i] < 0.0f && dbt[i] < 0.0f)
				return false;
			if(bt[i] < 0.0f && dbt[i] > 0.0f)
				b=std::max(b,(-bt[i])/dbt[i]);
			if(bt[i] > 0.0f && dbt[i] < 0.0f)
				t=std::min(t,bt[i]/(-dbt[i]));
		}
		if(b>t)return false;
		bt+=dbt*b;
		bl=PT2(1.0f-b,b);
		return true;
	}
    void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT4& bc) const {
        scalar sqrDistanceTmp;
        PT cpTmp;
        PT b;

        // Construct the planes for the faces of the tetrahedron.  The normals
        // are outer pointing, but specified not to be unit length.  We only need
        // to know sidedness of the query point, so we will save cycles by not
        // computing unit-length normals.
        PlaneTpl<T> plane;
        Vec3i indices;

        // Determine which faces are visible to the query point.  Only these
        // need to be processed by point-to-triangle distance queries.
        sqrDistance=ScalarUtil<T>::scalar_max;
        PT minTetraClosest=PT::Zero();
        for(sizeType i=0; i<4; ++i) {
            getPlane(i,plane,indices);
            if(plane.side(pt) >= 0) {
                TriangleTpl<T> tri(getNode(indices[0]),getNode(indices[1]),getNode(indices[2]));
                tri.calcPointDist(pt,sqrDistanceTmp,cpTmp,b);
                if (sqrDistanceTmp < sqrDistance) {
                    sqrDistance=sqrDistanceTmp;
                    cp=cpTmp;
                    bc=PT4::Zero();
                    bc[indices[0]]=b[0];
                    bc[indices[1]]=b[1];
                    bc[indices[2]]=b[2];
                }
            }
        }

        if (sqrDistance == ScalarUtil<T>::scalar_max) {
            sqrDistance=(T)0.0f;
            cp=pt;
            bc=bary(pt);
        }
    }
    void getPlane(const sizeType& i,PlaneTpl<T>& p,Vec3i& ind) const {
        //a,c,b
        //a,d,c
        //a,b,d
        //b,c,d
        switch(i) {
        case 0:
            p=PlaneTpl<T>(_a,_c,_b);
            ind=Vec3i(0,2,1);
            break;
        case 1:
            p=PlaneTpl<T>(_a,_d,_c);
            ind=Vec3i(0,3,2);
            break;
        case 2:
            p=PlaneTpl<T>(_a,_b,_d);
            ind=Vec3i(0,1,3);
            break;
        case 3:
            p=PlaneTpl<T>(_b,_c,_d);
            ind=Vec3i(1,2,3);
            break;
        }
    }
    const PT& getNode(const sizeType& i) const {
        return (&_a)[i];
    }
    void writeVTK(VTKWriter<T>& os) const
    {
        os.setRelativeIndex();
        std::vector<PT,Eigen::aligned_allocator<PT> > vss;
        vss.push_back(_a);
        vss.push_back(_b);
        vss.push_back(_c);
        vss.push_back(_d);
        os.appendPoints(vss.begin(),vss.end());

        std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
        iss.push_back(Vec3i(0,1,0));
        iss.push_back(Vec3i(0,2,0));
        iss.push_back(Vec3i(0,3,0));
        iss.push_back(Vec3i(1,2,0));
        iss.push_back(Vec3i(2,3,0));
        iss.push_back(Vec3i(3,0,0));
        os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
    }
public:
    //data
    PT _a;
    PT _b;
    PT _c;
    PT _d;
    bool _swap;
};
template <typename T,int dim>
class OBBTpl;
template <typename T>
class OBBTpl<T,2>
{
public:
    typedef typename Eigen::Matrix<T,2,1> PT;
    typedef typename Eigen::Matrix<T,3,1> PT3;
    typedef typename Eigen::Matrix<T,2,2> MAT;
    typedef typename Eigen::Matrix<T,3,3> MAT3;
public:
    OBBTpl(){}
    OBBTpl(const BBox<T,2>& bb)
    {
        _ext=bb.getExtent()*0.5f;
        _rot.setIdentity();
        _trans=bb._maxC-_ext;
    }
    OBBTpl(const BBox<T,3>& bb)
    {
        new (this)OBBTpl(BBox<T,2>(PT(bb._minC[0],bb._minC[1]),PT(bb._maxC[0],bb._maxC[1])));
    }
    OBBTpl(const MAT& rot,const PT& trans,const BBox<T,2>& bb)
    {
        new (this)OBBTpl(bb);
        _rot=rot;
        _trans=_rot*_trans+trans;
    }
    OBBTpl(const MAT3& rot,const PT3& trans,const BBox<T,3>& bb)
    {
        new (this)OBBTpl(rot.block<2,2>(0,0),trans.block<2,1>(0,0),BBox<T,2>(PT(bb._minC[0],bb._minC[1]),PT(bb._maxC[0],bb._maxC[1])));
    }
    OBBTpl(const MAT& rot,const PT& trans,const PT& ext):_rot(rot),_trans(trans),_ext(ext){}
    OBBTpl(const MAT3& rot,const PT3& trans,const PT3& ext)
    {
        _rot=rot.block<2,2>(0,0);
        _trans=trans.block<2,1>(0,0);
        _ext=ext.block<2,1>(0,0);
    }
    bool intersect(const OBBTpl<T,2>& other) const
    {
        MAT rotB2A=_rot.transpose()*other._rot;
        PT transB2A=_rot.transpose()*(other._trans-_trans);
        PT dirX=rotB2A.col(0)*other._ext.x();
        PT dirY=rotB2A.col(1)*other._ext.y();
        PT dir(std::abs(dirX.x())+std::abs(dirY.x()),std::abs(dirX.y())+std::abs(dirY.y()));
        if(SAT<T>::testSegment(PT(-_ext.x(),_ext.x()),PT(transB2A.x()-dir.x(),transB2A.x()+dir.x())))
            return false;
        if(SAT<T>::testSegment(PT(-_ext.y(),_ext.y()),PT(transB2A.y()-dir.y(),transB2A.y()+dir.y())))
            return false;

        transB2A=other._rot.transpose()*(_trans-other._trans);
        dirX=rotB2A.row(0)*_ext.x();
        dirY=rotB2A.row(1)*_ext.y();
        dir=PT(std::abs(dirX.x())+std::abs(dirY.x()),std::abs(dirX.y())+std::abs(dirY.y()));
        if(SAT<T>::testSegment(PT(-other._ext.x(),other._ext.x()),PT(transB2A.x()-dir.x(),transB2A.x()+dir.x())))
            return false;
        if(SAT<T>::testSegment(PT(-other._ext.y(),other._ext.y()),PT(transB2A.y()-dir.y(),transB2A.y()+dir.y())))
            return false;
        return true;
    }
    bool intersect(const BBox<T,2>& other) const{return intersect(OBBTpl<T,2>(other));}
    void writeVTK(VTKWriter<T>& os) const
    {
        os.setRelativeIndex();
        std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
        PT pt1=_rot*PT(-_ext.x(),-_ext.y())+_trans;
        PT pt2=_rot*PT( _ext.x(),-_ext.y())+_trans;
        PT pt3=_rot*PT( _ext.x(), _ext.y())+_trans;
        PT pt4=_rot*PT(-_ext.x(), _ext.y())+_trans;
        vss.push_back(PT3(pt1.x(),pt1.y(),0.0f));
        vss.push_back(PT3(pt2.x(),pt2.y(),0.0f));
        vss.push_back(PT3(pt3.x(),pt3.y(),0.0f));
        vss.push_back(PT3(pt4.x(),pt4.y(),0.0f));
        os.appendPoints(vss.begin(),vss.end());

        std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
        iss.push_back(Vec3i(0,1,0));
        iss.push_back(Vec3i(1,2,0));
        iss.push_back(Vec3i(2,3,0));
        iss.push_back(Vec3i(3,0,0));
        os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
    }
public:
    MAT _rot;
    PT _trans;
    PT _ext;
};
template <typename T>
class OBBTpl<T,3>
{
public:
    typedef typename Eigen::Matrix<T,3,1> PT;
    typedef typename Eigen::Matrix<T,2,1> PT2;
    typedef typename Eigen::Matrix<T,3,3> MAT;
public:
    OBBTpl(){}
    OBBTpl(const BBox<T,3>& bb)
    {
        _ext=bb.getExtent()*0.5f;
        _rot.setIdentity();
        _trans=bb._maxC-_ext;
    }
    OBBTpl(const MAT& rot,const PT& trans,const BBox<T,3>& bb)
    {
        new (this)OBBTpl(bb);
        _rot=rot;
        _trans=_rot*_trans+trans;
    }
    OBBTpl(const MAT& rot,const PT& trans,const PT& ext):_rot(rot),_trans(trans),_ext(ext){}
    bool intersect(const OBBTpl<T,3>& other) const
    {
        MAT rotB2A=_rot.transpose()*other._rot;
        PT transB2A=_rot.transpose()*(other._trans-_trans);
        PT dirs[3];
        
        //set two
        PT ctrAInB=other._rot.transpose()*(_trans-other._trans);
        dirs[0]=rotB2A.row(0)*_ext.x();
        dirs[1]=rotB2A.row(1)*_ext.y();
        dirs[2]=rotB2A.row(2)*_ext.z();
        PT dir=PT(std::abs(dirs[0].x())+std::abs(dirs[1].x())+std::abs(dirs[2].x()),
                  std::abs(dirs[0].y())+std::abs(dirs[1].y())+std::abs(dirs[2].y()),
                  std::abs(dirs[0].z())+std::abs(dirs[1].z())+std::abs(dirs[2].z()));
        if(SAT<T>::testSegment(PT2(-other._ext.x(),other._ext.x()),PT2(ctrAInB.x()-dir.x(),ctrAInB.x()+dir.x())))
            return false;
        if(SAT<T>::testSegment(PT2(-other._ext.y(),other._ext.y()),PT2(ctrAInB.y()-dir.y(),ctrAInB.y()+dir.y())))
            return false;
        if(SAT<T>::testSegment(PT2(-other._ext.z(),other._ext.z()),PT2(ctrAInB.z()-dir.z(),ctrAInB.z()+dir.z())))
            return false;

        //set two
        dirs[0]=rotB2A.col(0)*other._ext.x();
        dirs[1]=rotB2A.col(1)*other._ext.y();
        dirs[2]=rotB2A.col(2)*other._ext.z();
        dir=PT(std::abs(dirs[0].x())+std::abs(dirs[1].x())+std::abs(dirs[2].x()),
               std::abs(dirs[0].y())+std::abs(dirs[1].y())+std::abs(dirs[2].y()),
               std::abs(dirs[0].z())+std::abs(dirs[1].z())+std::abs(dirs[2].z()));
        if(SAT<T>::testSegment(PT2(-_ext.x(),_ext.x()),PT2(transB2A.x()-dir.x(),transB2A.x()+dir.x())))
            return false;
        if(SAT<T>::testSegment(PT2(-_ext.y(),_ext.y()),PT2(transB2A.y()-dir.y(),transB2A.y()+dir.y())))
            return false;
        if(SAT<T>::testSegment(PT2(-_ext.z(),_ext.z()),PT2(transB2A.z()-dir.z(),transB2A.z()+dir.z())))
            return false;
        
#define CROSS0(A) 0.0f,-A[2],A[1]
#define CROSS1(A) A[2],0.0f,-A[0]
#define CROSS2(A) -A[1],A[0],0.0f

#define ABSDOT0(A,B) std::abs(A[1])*B[1]+std::abs(A[2])*B[2]
#define ABSDOT1(A,B) std::abs(A[0])*B[0]+std::abs(A[2])*B[2]
#define ABSDOT2(A,B) std::abs(A[0])*B[0]+std::abs(A[1])*B[1]

#define ABSDOTB0(A,B) std::abs(A[1].dot(B))+std::abs(A[2].dot(B))
#define ABSDOTB1(A,B) std::abs(A[0].dot(B))+std::abs(A[2].dot(B))
#define ABSDOTB2(A,B) std::abs(A[0].dot(B))+std::abs(A[1].dot(B))

#define CROSS_TEST(AXIS)    \
for(sizeType j=0;j<3;j++)    \
{    \
    PT axis(CROSS##AXIS(rotB2A.col(j)));    \
    if(axis.squaredNorm() > ScalarUtil<T>::scalar_eps)    \
    {    \
        T extAPrj=ABSDOT##AXIS(axis,_ext);    \
        T ctrBPrj=transB2A.dot(axis);    \
        T extBPrj=ABSDOTB##AXIS(dirs,axis);    \
        if(SAT<T>::testSegment(PT2(-extAPrj,extAPrj),PT2(ctrBPrj-extBPrj,ctrBPrj+extBPrj)))    \
            return false;    \
    }    \
}

#undef CROSS0
#undef CROSS1
#undef CROSS2

#undef ABSDOT0
#undef ABSDOT1
#undef ABSDOT2

#undef ABSDOTB0
#undef ABSDOTB1
#undef ABSDOTB2
        return true;
    }
    bool intersect(const BBox<T,3>& other) const{return intersect(OBBTpl<T,3>(other));}
    void writeVTK(VTKWriter<T>& os) const
    {
        os.setRelativeIndex();
        std::vector<PT,Eigen::aligned_allocator<PT> > vss;
        vss.push_back(_rot*PT(-_ext.x(),-_ext.y(),-_ext.z())+_trans);
        vss.push_back(_rot*PT( _ext.x(),-_ext.y(),-_ext.z())+_trans);
        vss.push_back(_rot*PT( _ext.x(), _ext.y(),-_ext.z())+_trans);
        vss.push_back(_rot*PT(-_ext.x(), _ext.y(),-_ext.z())+_trans);
        vss.push_back(_rot*PT(-_ext.x(),-_ext.y(), _ext.z())+_trans);
        vss.push_back(_rot*PT( _ext.x(),-_ext.y(), _ext.z())+_trans);
        vss.push_back(_rot*PT( _ext.x(), _ext.y(), _ext.z())+_trans);
        vss.push_back(_rot*PT(-_ext.x(), _ext.y(), _ext.z())+_trans);
        os.appendPoints(vss.begin(),vss.end());

        std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
        iss.push_back(Vec3i(0,1,0));
        iss.push_back(Vec3i(1,2,0));
        iss.push_back(Vec3i(2,3,0));
        iss.push_back(Vec3i(3,0,0));
        iss.push_back(Vec3i(4,5,0));
        iss.push_back(Vec3i(5,6,0));
        iss.push_back(Vec3i(6,7,0));
        iss.push_back(Vec3i(7,4,0));
        iss.push_back(Vec3i(0,4,0));
        iss.push_back(Vec3i(1,5,0));
        iss.push_back(Vec3i(2,6,0));
        iss.push_back(Vec3i(3,7,0));
        os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
    }
public:
    MAT _rot;
    PT _trans;
    PT _ext;
};
template <typename T>
class KDOP18 : public Serializable
{
public:
	typedef typename Eigen::Matrix<T,3,1> PT;
public:
	FORCE_INLINE KDOP18():Serializable(65535){empty();}
	FORCE_INLINE KDOP18(const PT& v):Serializable(65535) {
		_dist[0] = _dist[9]  = v[0];
		_dist[1] = _dist[10] = v[1];
		_dist[2] = _dist[11] = v[2];

		T d3, d4, d5, d6, d7, d8;
		getDistances(v, d3, d4, d5, d6, d7, d8);
		_dist[3] = _dist[12] = d3;
		_dist[4] = _dist[13] = d4;
		_dist[5] = _dist[14] = d5;
		_dist[6] = _dist[15] = d6;
		_dist[7] = _dist[16] = d7;
		_dist[8] = _dist[17] = d8;
 	}
	FORCE_INLINE KDOP18(const PT& a,const PT& b):Serializable(65535) {
		_dist[0]  = std::min(a[0], b[0]);
		_dist[9]  = std::max(a[0], b[0]);
		_dist[1]  = std::min(a[1], b[1]);
		_dist[10] = std::max(a[1], b[1]);
		_dist[2]  = std::min(a[2], b[2]);
		_dist[11] = std::max(a[2], b[2]);

		T ad3, ad4, ad5, ad6, ad7, ad8;
		getDistances(a, ad3, ad4, ad5, ad6, ad7, ad8);
		T bd3, bd4, bd5, bd6, bd7, bd8;
		getDistances(b, bd3, bd4, bd5, bd6, bd7, bd8);
		_dist[3]  = std::min(ad3, bd3);
		_dist[12] = std::max(ad3, bd3);
		_dist[4]  = std::min(ad4, bd4);
		_dist[13] = std::max(ad4, bd4);
		_dist[5]  = std::min(ad5, bd5);
		_dist[14] = std::max(ad5, bd5);
		_dist[6]  = std::min(ad6, bd6);
		_dist[15] = std::max(ad6, bd6);
		_dist[7]  = std::min(ad7, bd7);
		_dist[16] = std::max(ad7, bd7);
		_dist[8]  = std::min(ad8, bd8);
		_dist[17] = std::max(ad8, bd8);
	}
	FORCE_INLINE bool write(ostream& os) const{
		for(char v=0;v<18;v++)
			writeBinaryData(_dist[v],os);
		return os.good();
	}
	FORCE_INLINE bool read(istream& is) const{
		for(char v=0;v<18;v++)
			readBinaryData(_dist[v],is);
		return is.good();
	}
	FORCE_INLINE boost::shared_ptr<Serializable> copy() const {
		return boost::shared_ptr<Serializable>(new KDOP18);
	}
	FORCE_INLINE void reset(){empty();}
	FORCE_INLINE void empty() {
		for (int i=0; i<9; i++) {
			_dist[i] = ScalarUtil<T>::scalar_max;
			_dist[i+9] = -ScalarUtil<T>::scalar_max;
		}
	}
	FORCE_INLINE void enlarged(T len){
		for(int i=0;i<3;i++){
			_dist[i]-=len;
			_dist[i+9]+=len;
		}
		for(int i=3;i<9;i++){
			_dist[i]-=len*2.0f;
			_dist[i+9]+=len*2.0f;
		}
	}
	FORCE_INLINE void setPoints(const PT& a,const PT& b,const PT& c){
		new(this) KDOP18(a,b);
		setUnion(c);
	}
	FORCE_INLINE void setUnion(const PT& p){
		_dist[0]  = std::min(p[0], _dist[0]);
		_dist[9]  = std::max(p[0], _dist[9]);
		_dist[1]  = std::min(p[1], _dist[1]);
		_dist[10] = std::max(p[1], _dist[10]);
		_dist[2]  = std::min(p[2], _dist[2]);
		_dist[11] = std::max(p[2], _dist[11]);

		T d3, d4, d5, d6, d7, d8;
		getDistances(p, d3, d4, d5, d6, d7, d8);
		_dist[3]  = std::min(d3, _dist[3]);
		_dist[12] = std::max(d3, _dist[12]);
		_dist[4]  = std::min(d4, _dist[4]);
		_dist[13] = std::max(d4, _dist[13]);
		_dist[5]  = std::min(d5, _dist[5]);
		_dist[14] = std::max(d5, _dist[14]);
		_dist[6]  = std::min(d6, _dist[6]);
		_dist[15] = std::max(d6, _dist[15]);
		_dist[7]  = std::min(d7, _dist[7]);
		_dist[16] = std::max(d7, _dist[16]);
		_dist[8]  = std::min(d8, _dist[8]);
		_dist[17] = std::max(d8, _dist[17]);
	}
	FORCE_INLINE void setUnion(const KDOP18& b){
		_dist[0]  = std::min(b._dist[0], _dist[0]);
		_dist[9]  = std::max(b._dist[9], _dist[9]);
		_dist[1]  = std::min(b._dist[1], _dist[1]);
		_dist[10] = std::max(b._dist[10], _dist[10]);
		_dist[2]  = std::min(b._dist[2], _dist[2]);
		_dist[11] = std::max(b._dist[11], _dist[11]);
		_dist[3]  = std::min(b._dist[3], _dist[3]);
		_dist[12] = std::max(b._dist[12], _dist[12]);
		_dist[4]  = std::min(b._dist[4], _dist[4]);
		_dist[13] = std::max(b._dist[13], _dist[13]);
		_dist[5]  = std::min(b._dist[5], _dist[5]);
		_dist[14] = std::max(b._dist[14], _dist[14]);
		_dist[6]  = std::min(b._dist[6], _dist[6]);
		_dist[15] = std::max(b._dist[15], _dist[15]);
		_dist[7]  = std::min(b._dist[7], _dist[7]);
		_dist[16] = std::max(b._dist[16], _dist[16]);
		_dist[8]  = std::min(b._dist[8], _dist[8]);
		_dist[17] = std::max(b._dist[17], _dist[17]);
	}
	FORCE_INLINE PT minCorner() const{return PT(_dist[0],_dist[1 ],_dist[2 ]);}
	FORCE_INLINE PT maxCorner() const{return PT(_dist[9],_dist[10],_dist[11]);}
	FORCE_INLINE bool intersect(const KDOP18& b,sizeType DIM=3) const{
		ASSERT(DIM == 3)
		for (int i=0; i<9; i++) {
			if (_dist[i] > b._dist[i+9]) return false;
			if (_dist[i+9] < b._dist[i]) return false;
		}
		return true;
	}
	FORCE_INLINE bool intersect(const BBox<T>& b,sizeType DIM=3) const{
		ASSERT(DIM == 3)
		BBox<T> bb(minCorner(),maxCorner());
		return bb.intersect(b);
	}
	FORCE_INLINE bool contain(const PT& p) const
	{
		for (int i=0; i<3; i++) {
			if (p[i] < _dist[i] || p[i] > _dist[i+9])
				return false;
		}

		T d[6];
		getDistances(p, d);
		for (int i=3; i<9; i++) {
			if (d[i-3] < _dist[i] || d[i-3] > _dist[i+9])
				return false;
		}

		return true;
	}
protected:
	FORCE_INLINE static void getDistances(const PT& p,T &d3, T &d4, T &d5, T &d6, T &d7, T &d8)
	{
		d3 = p[0] + p[1];
		d4 = p[0] + p[2];
		d5 = p[1] + p[2];
		d6 = p[0] - p[1];
		d7 = p[0] - p[2];
		d8 = p[1] - p[2];
	}
	FORCE_INLINE static void getDistances(const PT& p, T d[])
	{
		d[0] = p[0] + p[1];
		d[1] = p[0] + p[2];
		d[2] = p[1] + p[2];
		d[3] = p[0] - p[1];
		d[4] = p[0] - p[2];
		d[5] = p[1] - p[2];
	}
	FORCE_INLINE static T getDistances(const PT &p, int i)
	{
		if (i == 0) return p[0]+p[1];
		if (i == 1) return p[0]+p[2];
		if (i == 2) return p[1] + p[2];
		if (i == 3) return p[0] - p[1];
		if (i == 4) return p[0] - p[2];
		if (i == 5) return p[1] - p[2];
		return 0;
	}
	T _dist[18];
};
template <typename T>
class Triangle2Triangle
{
public:
	#define ZERO 1E-12
	typedef typename Eigen::Matrix<T,3,1> PT;
	typedef typename Eigen::Matrix<T,2,1> PT2;
	typedef typename Eigen::Matrix<T,3,2> MT32;
	//-1 coplane or parallel; 0 not intersect; 1 intersect
	int intersect(const PT t0[3], const PT t1[3], PT& e0, PT& e1)
	{
		PlaneTpl<T> p0 = PlaneTpl<T>(t0[0], t0[1], t0[2]);
		PlaneTpl<T> p1 = PlaneTpl<T>(t1[0], t1[1], t1[2]);

		PT dir0 = p0._n.cross(p1._n);
		if(dir0.norm() < ZERO)
			return -1;
		dir0 /= dir0.norm();
		Eigen::Matrix<T, 3, 3> A;
		A.block<1,3>(0,0) = p0._n.transpose();
		A.block<1,3>(1,0) = p1._n.transpose();
		A.block<1,3>(2,0) = dir0.transpose();
		PT rh(p0._n.dot(p0._x0), p1._n.dot(p1._x0), 0);
		PT b0 = A.inverse()*rh;
		vector<T> tv;
		T maxp[2] = {T(-1E30f), T(-1E30f)};
		T minp[2] = {T(1E30f), T(1E30f)};
		for(int i = 0; i < 2; ++i)
		{
			PT pp[3];
			if(i == 0){pp[0] = t0[0]; pp[1] = t0[1]; pp[2] = t0[2];}
			else{pp[0] = t1[0]; pp[1] = t1[1]; pp[2] = t1[2];}
			for(int j = 0; j < 3; ++j)
			{
				PT dir1 = pp[(j+1)%3] - pp[j];
				PT b1 = pp[j];
				MT32 B;
				B.block<3,1>(0,0) = dir0;
				B.block<3,1>(0,1) = -dir1;
				PT rhb = b1 - b0;
				if(std::abs((B.transpose()*B).determinant()) < ZERO)
					continue;
				PT2 VT = (B.transpose()*B).inverse()*B.transpose()*rhb;
				//cout <<VT(0) <<' ' <<VT(1) <<endl;
				if((B*VT - rhb).norm() < ZERO && VT(1) > T(0.0f) && VT(1) < T(1.0f))
				{
					tv.push_back(VT(0));
					if(VT(0) > maxp[i])
						maxp[i] = VT(0);
					if(VT(0) < minp[i])
						minp[i] = VT(0);
				}
			}
		}
		//cout <<t0[0].transpose() <<' ' <<t0[1].transpose() <<' ' <<t0[2].transpose() <<endl;
		//cout <<t1[0].transpose() <<' ' <<t1[1].transpose() <<' ' <<t1[2].transpose() <<endl;
		if(tv.size() == 0)
			return 0;
		if(minp[0] > maxp[1] || minp[1] > maxp[0])
			return 0;
		sort(tv.begin(), tv.end());
		assert(tv.size() >= 2);
		e0 = b0 + tv[(tv.size() - 1)/2]*dir0;
		e1 = b0 + tv[tv.size()/2]*dir0;
		return 1;
	}
	#undef ZERO
};
template <typename T>
class Sphere
{
public:
	typedef typename Eigen::Matrix<T,3,1> PT;
	Sphere(){}
	Sphere(const PT& ctr,const T& rad)
		:_ctr(ctr),_rad(rad){}
	bool intersect(const PT& a,const PT& b) const{
		T A=(b-a).squaredNorm();
		if(A < ScalarUtil<T>::scalar_eps)
			return (_ctr-a).norm() < _rad;
		
		T B=2.0f*(b-a).dot(a-_ctr);
		T C=(a-_ctr).squaredNorm()-_rad*_rad;
		T delta=B*B-4.0f*A*C;
		if(delta < 0.0f)
			return false;
		
		T L=(-B-sqrt(delta))/(2*A);
		T H=(-B+sqrt(delta))/(2*A);
		return L < 1 && H > 0;
	}
	//data
	PT _ctr;
	T _rad;
};

typedef LineSegTpl<scalar> LineSeg;
typedef PlaneTpl<scalar> Plane;
typedef TriangleTpl<scalar> Triangle;
typedef TetrahedronTpl<scalar> Tetrahedron;
typedef OBBTpl<scalar,2> OBB2D;
typedef OBBTpl<scalar,3> OBB3D;

PRJ_END

#endif
