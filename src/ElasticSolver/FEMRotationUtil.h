#ifndef FEM_ROTATION_UTIL_H
#define FEM_ROTATOIN_UTIL_H

#include "MathBasic.h"

PRJ_BEGIN

template <typename T>
static FORCE_INLINE typename ScalarUtil<T>::ScalarMat3 cross(const typename ScalarUtil<T>::ScalarVec3& v)
{
    typename ScalarUtil<T>::ScalarMat3 ret;
    ret.setZero();
    ret(0,1)=-v[2];
    ret(0,2)=v[1];
    ret(1,2)=-v[0];
    ret(1,0)=v[2];
    ret(2,0)=-v[1];
    ret(2,1)=v[0];
    return ret;
}
template <typename T>
static FORCE_INLINE typename Eigen::Matrix<T,9,9> rotationCoef(const typename ScalarUtil<T>::ScalarMat3& r)
{
    typename Eigen::Matrix<T,9,9> rCoef;
    rCoef.setZero();
#define GI(I,J) J*3+I
    for(sizeType i=0; i<3; i++)
        for(sizeType j=0; j<3; j++)
            for(sizeType k=0; k<3; k++)
                rCoef(GI(i,j),GI(k,j))=r(i,k);
    return rCoef;
#undef GI
}
template <typename MAT3>
static FORCE_INLINE Mat3d expWGrad(const Vec3d& W,MAT3* DexpWDWx=NULL,MAT3* DexpWDWy=NULL,MAT3* DexpWDWz=NULL)
{
    scalarD WLx=W[0];
    scalarD WLy=W[1];
    scalarD WLz=W[2];
    Mat3d expW;

    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;
    scalarD tt19;
    scalarD tt20;
    scalarD tt21;
    scalarD tt22;
    scalarD tt23;
    scalarD tt24;
    scalarD tt25;
    scalarD tt26;
    scalarD tt27;
    scalarD tt28;
    scalarD tt29;
    scalarD tt30;
    scalarD tt31;
    scalarD tt32;
    scalarD tt33;
    scalarD tt34;
    scalarD tt35;
    scalarD tt36;
    scalarD tt37;
    scalarD tt38;
    scalarD tt39;
    scalarD tt40;
    scalarD tt41;
    scalarD tt42;
    scalarD tt43;
    scalarD tt44;
    scalarD tt45;
    scalarD tt46;
    scalarD tt47;
    scalarD tt48;
    scalarD tt49;
    scalarD tt50;
    scalarD tt51;
    scalarD tt52;
    scalarD tt53;
    scalarD tt54;
    scalarD tt55;
    scalarD tt56;
    scalarD tt57;
    scalarD tt58;
    scalarD tt59;
    scalarD tt60;
    scalarD tt61;
    scalarD tt62;
    scalarD tt63;
    scalarD tt64;
    scalarD tt65;
    scalarD tt66;
    scalarD tt67;
    scalarD tt68;
    scalarD tt69;
    scalarD tt70;
    scalarD tt71;
    scalarD tt72;
    scalarD tt73;
    scalarD tt74;
    scalarD tt75;
    scalarD tt76;
    scalarD tt77;
    scalarD tt78;
    scalarD tt79;
    scalarD tt80;
    scalarD tt81;
    scalarD tt82;
    scalarD tt83;
    scalarD tt84;

    tt1=pow(WLx,2);
    tt2=pow(WLy,2);
    tt3=pow(WLz,2);
    tt4=std::max(tt3+tt2+tt1,1E-12);
    tt5=sqrt(tt4);
    tt6=0.5*tt5;
    tt7=cos(tt6);
    tt8=pow(tt7,2);
    tt9=1/tt4;
    tt10=sin(tt6);
    tt11=pow(tt10,2);
    tt12=-tt2*tt9*tt11;
    tt13=-tt3*tt9*tt11;
    tt14=1/tt5;
    tt15=WLx*WLy*tt9*tt11;
    tt16=WLx*WLz*tt9*tt11;
    tt17=-tt1*tt9*tt11;
    tt18=WLy*WLz*tt9*tt11;
    expW(0,0)=tt13+tt12+tt1*tt9*tt11+tt8;
    expW(0,1)=2*(tt15-WLz*tt14*tt7*tt10);
    expW(0,2)=2*(tt16+WLy*tt14*tt7*tt10);
    expW(1,0)=2*(tt15+WLz*tt14*tt7*tt10);
    expW(1,1)=tt13+tt2*tt9*tt11+tt17+tt8;
    expW(1,2)=2*(tt18-WLx*tt14*tt7*tt10);
    expW(2,0)=2*(tt16-WLy*tt14*tt7*tt10);
    expW(2,1)=2*(tt18+WLx*tt14*tt7*tt10);
    expW(2,2)=tt3*tt9*tt11+tt12+tt17+tt8;
    if(!DexpWDWx)return expW;

    tt19=pow(WLx,3);
    tt20=1/pow(tt5,3);
    tt21=-1.0*WLx*tt2*tt20*tt7*tt10;
    tt22=-1.0*WLx*tt3*tt20*tt7*tt10;
    tt23=-1.0*WLx*tt14*tt7*tt10;
    tt24=1/pow(tt4,2);
    tt25=2*WLx*tt2*tt24*tt11;
    tt26=2*WLx*tt3*tt24*tt11;
    tt27=-0.5*WLx*WLz*tt9*tt8;
    tt28=1.0*tt1*WLy*tt20*tt7*tt10;
    tt29=WLx*WLz*tt20*tt7*tt10;
    tt30=-2*tt1*WLy*tt24*tt11;
    tt31=WLy*tt9*tt11;
    tt32=0.5*WLx*WLz*tt9*tt11;
    tt33=0.5*WLx*WLy*tt9*tt8;
    tt34=-WLx*WLy*tt20*tt7*tt10;
    tt35=1.0*tt1*WLz*tt20*tt7*tt10;
    tt36=-2*tt1*WLz*tt24*tt11;
    tt37=-0.5*WLx*WLy*tt9*tt11;
    tt38=WLz*tt9*tt11;
    tt39=0.5*WLx*WLz*tt9*tt8;
    tt40=-WLx*WLz*tt20*tt7*tt10;
    tt41=-0.5*WLx*WLz*tt9*tt11;
    tt42=-1.0*tt19*tt20*tt7*tt10;
    tt43=1.0*WLx*tt2*tt20*tt7*tt10;
    tt44=2*tt19*tt24*tt11;
    tt45=-2*WLx*tt2*tt24*tt11;
    tt46=-2*WLx*tt9*tt11;
    tt47=1.0*WLx*WLy*WLz*tt20*tt7*tt10;
    tt48=-tt14*tt7*tt10;
    tt49=-2*WLx*WLy*WLz*tt24*tt11;
    tt50=-0.5*WLx*WLy*tt9*tt8;
    tt51=WLx*WLy*tt20*tt7*tt10;
    tt52=0.5*WLx*WLy*tt9*tt11;
    tt53=tt14*tt7*tt10;
    tt54=1.0*WLx*tt3*tt20*tt7*tt10;
    tt55=-2*WLx*tt3*tt24*tt11;
    tt56=pow(WLy,3);
    tt57=-1.0*tt56*tt20*tt7*tt10;
    tt58=-1.0*WLy*tt3*tt20*tt7*tt10;
    tt59=-1.0*WLy*tt14*tt7*tt10;
    tt60=2*tt56*tt24*tt11;
    tt61=2*WLy*tt3*tt24*tt11;
    tt62=-2*WLy*tt9*tt11;
    tt63=-0.5*WLy*WLz*tt9*tt8;
    tt64=WLy*WLz*tt20*tt7*tt10;
    tt65=WLx*tt9*tt11;
    tt66=0.5*WLy*WLz*tt9*tt11;
    tt67=0.5*WLy*WLz*tt9*tt8;
    tt68=-WLy*WLz*tt20*tt7*tt10;
    tt69=-0.5*WLy*WLz*tt9*tt11;
    tt70=-1.0*tt1*WLy*tt20*tt7*tt10;
    tt71=2*tt1*WLy*tt24*tt11;
    tt72=1.0*tt2*WLz*tt20*tt7*tt10;
    tt73=-2*tt2*WLz*tt24*tt11;
    tt74=1.0*WLy*tt3*tt20*tt7*tt10;
    tt75=-2*WLy*tt3*tt24*tt11;
    tt76=-1.0*tt2*WLz*tt20*tt7*tt10;
    tt77=pow(WLz,3);
    tt78=-1.0*tt77*tt20*tt7*tt10;
    tt79=-1.0*WLz*tt14*tt7*tt10;
    tt80=2*tt2*WLz*tt24*tt11;
    tt81=2*tt77*tt24*tt11;
    tt82=-2*WLz*tt9*tt11;
    tt83=-1.0*tt1*WLz*tt20*tt7*tt10;
    tt84=2*tt1*WLz*tt24*tt11;
    (*DexpWDWx)(0,0)=2*WLx*tt9*tt11+tt26+tt25-2*tt19*tt24*tt11+tt23+tt22+tt21+1.0*tt19*tt20*tt7*tt10;
    (*DexpWDWx)(0,1)=2*(tt32+tt31+tt30+tt29+tt28+tt27);
    (*DexpWDWx)(0,2)=2*(tt38+tt37+tt36+tt35+tt34+tt33);
    (*DexpWDWx)(1,0)=2*(tt41+tt31+tt30+tt40+tt28+tt39);
    (*DexpWDWx)(1,1)=tt46+tt26+tt45+tt44+tt23+tt22+tt43+tt42;
    (*DexpWDWx)(1,2)=2*(0.5*tt1*tt9*tt11+tt49+tt48+tt47+tt1*tt20*tt7*tt10-0.5*tt1*tt9*tt8);
    (*DexpWDWx)(2,0)=2*(tt38+tt52+tt36+tt35+tt51+tt50);
    (*DexpWDWx)(2,1)=2*(-0.5*tt1*tt9*tt11+tt49+tt53+tt47-tt1*tt20*tt7*tt10+0.5*tt1*tt9*tt8);
    (*DexpWDWx)(2,2)=tt46+tt55+tt25+tt44+tt23+tt54+tt21+tt42;
    (*DexpWDWy)(0,0)=tt62+tt61+tt60+tt30+tt59+tt58+tt57+tt28;
    (*DexpWDWy)(0,1)=2*(tt66+tt65+tt45+tt64+tt43+tt63);
    (*DexpWDWy)(0,2)=2*(-0.5*tt2*tt9*tt11+tt49+tt53+tt47-tt2*tt20*tt7*tt10+0.5*tt2*tt9*tt8);
    (*DexpWDWy)(1,0)=2*(tt69+tt65+tt45+tt68+tt43+tt67);
    (*DexpWDWy)(1,1)=2*WLy*tt9*tt11+tt61-2*tt56*tt24*tt11+tt71+tt59+tt58+1.0*tt56*tt20*tt7*tt10+tt70;
    (*DexpWDWy)(1,2)=2*(tt38+tt52+tt73+tt72+tt51+tt50);
    (*DexpWDWy)(2,0)=2*(0.5*tt2*tt9*tt11+tt49+tt48+tt47+tt2*tt20*tt7*tt10-0.5*tt2*tt9*tt8);
    (*DexpWDWy)(2,1)=2*(tt38+tt37+tt73+tt72+tt34+tt33);
    (*DexpWDWy)(2,2)=tt62+tt75+tt60+tt71+tt59+tt74+tt57+tt70;
    (*DexpWDWz)(0,0)=tt82+tt81+tt80+tt36+tt79+tt78+tt76+tt35;
    (*DexpWDWz)(0,1)=2*(0.5*tt3*tt9*tt11+tt49+tt48+tt3*tt20*tt7*tt10+tt47-0.5*tt3*tt9*tt8);
    (*DexpWDWz)(0,2)=2*(tt69+tt65+tt55+tt54+tt68+tt67);
    (*DexpWDWz)(1,0)=2*(-0.5*tt3*tt9*tt11+tt49+tt53-tt3*tt20*tt7*tt10+tt47+0.5*tt3*tt9*tt8);
    (*DexpWDWz)(1,1)=tt82+tt81+tt73+tt84+tt79+tt78+tt72+tt83;
    (*DexpWDWz)(1,2)=2*(tt32+tt31+tt75+tt74+tt29+tt27);
    (*DexpWDWz)(2,0)=2*(tt66+tt65+tt55+tt54+tt64+tt63);
    (*DexpWDWz)(2,1)=2*(tt41+tt31+tt75+tt74+tt40+tt39);
    (*DexpWDWz)(2,2)=2*WLz*tt9*tt11-2*tt77*tt24*tt11+tt80+tt84+tt79+1.0*tt77*tt20*tt7*tt10+tt76+tt83;
    return expW;
}
static FORCE_INLINE Mat3d expWGradInte(const Vec3d& W)
{
    Mat3d RInt;
    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;

    tt1=pow(W[1],2);
    tt2=pow(W[0],2);
    tt3=pow(W[2],2);
    tt4=std::max(tt3+tt1+tt2,1E-12);
    tt5=1/tt4;
    tt6=-tt1*tt5;
    tt7=-tt3*tt5;
    tt8=sqrt(tt4);
    tt9=1-sin(tt8)/tt8;
    tt10=1-cos(tt8);
    tt11=W[0]*W[1]*tt5*tt9;
    tt12=W[0]*W[2]*tt5*tt9;
    tt13=-tt2*tt5;
    tt14=W[1]*W[2]*tt5*tt9;
    RInt(0,0)=(tt7+tt6)*tt9+1;
    RInt(0,1)=tt11-W[2]*tt5*tt10;
    RInt(0,2)=tt12+W[1]*tt5*tt10;
    RInt(1,0)=tt11+W[2]*tt5*tt10;
    RInt(1,1)=(tt7+tt13)*tt9+1;
    RInt(1,2)=tt14-W[0]*tt5*tt10;
    RInt(2,0)=tt12-W[1]*tt5*tt10;
    RInt(2,1)=tt14+W[0]*tt5*tt10;
    RInt(2,2)=(tt6+tt13)*tt9+1;
    return RInt;
}
static FORCE_INLINE Mat3d calcFGrad(const Mat3d& F,Eigen::Matrix<scalarD,9,9>* gradRSF)
{
    //input
    scalarD F11=F(0,0);
    scalarD F12=F(0,1);
    scalarD F13=F(0,2);
    scalarD F21=F(1,0);
    scalarD F22=F(1,1);
    scalarD F23=F(1,2);
    scalarD F31=F(2,0);
    scalarD F32=F(2,1);
    scalarD F33=F(2,2);
    Mat3d RSF;

    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;
    scalarD tt19;
    scalarD tt20;
    scalarD tt21;
    scalarD tt22;
    scalarD tt23;
    scalarD tt24;
    scalarD tt25;
    scalarD tt26;
    scalarD tt27;
    scalarD tt28;
    scalarD tt29;
    scalarD tt30;
    scalarD tt31;
    scalarD tt32;
    scalarD tt33;
    scalarD tt34;
    scalarD tt35;
    scalarD tt36;
    scalarD tt37;
    scalarD tt38;
    scalarD tt39;
    scalarD tt40;
    scalarD tt41;
    scalarD tt42;
    scalarD tt43;
    scalarD tt44;
    scalarD tt45;
    scalarD tt46;
    scalarD tt47;
    scalarD tt48;
    scalarD tt49;
    scalarD tt50;
    scalarD tt51;
    scalarD tt52;
    scalarD tt53;
    scalarD tt54;
    scalarD tt55;
    scalarD tt56;
    scalarD tt57;
    scalarD tt58;
    scalarD tt59;
    scalarD tt60;
    scalarD tt61;
    scalarD tt62;
    scalarD tt63;
    scalarD tt64;
    scalarD tt65;
    scalarD tt66;
    scalarD tt67;
    scalarD tt68;
    scalarD tt69;
    scalarD tt70;
    scalarD tt71;
    scalarD tt72;
    scalarD tt73;
    scalarD tt74;
    scalarD tt75;
    scalarD tt76;
    scalarD tt77;
    scalarD tt78;
    scalarD tt79;
    scalarD tt80;
    scalarD tt81;
    scalarD tt82;
    scalarD tt83;
    scalarD tt84;
    scalarD tt85;
    scalarD tt86;
    scalarD tt87;
    scalarD tt88;
    scalarD tt89;
    scalarD tt90;
    scalarD tt91;
    scalarD tt92;
    scalarD tt93;
    scalarD tt94;
    scalarD tt95;
    scalarD tt96;
    scalarD tt97;
    scalarD tt98;
    scalarD tt99;
    scalarD tt100;
    scalarD tt101;
    scalarD tt102;
    scalarD tt103;
    scalarD tt104;
    scalarD tt105;
    scalarD tt106;
    scalarD tt107;
    scalarD tt108;
    scalarD tt109;
    scalarD tt110;
    scalarD tt111;
    scalarD tt112;
    scalarD tt113;
    scalarD tt114;
    scalarD tt115;
    scalarD tt116;
    scalarD tt117;
    scalarD tt118;
    scalarD tt119;
    scalarD tt120;
    scalarD tt121;
    scalarD tt122;
    scalarD tt123;
    scalarD tt124;
    scalarD tt125;
    scalarD tt126;
    scalarD tt127;
    scalarD tt128;
    scalarD tt129;
    scalarD tt130;
    scalarD tt131;
    scalarD tt132;
    scalarD tt133;
    scalarD tt134;
    scalarD tt135;
    scalarD tt136;
    scalarD tt137;
    scalarD tt138;
    scalarD tt139;
    scalarD tt140;
    scalarD tt141;
    scalarD tt142;
    scalarD tt143;
    scalarD tt144;
    scalarD tt145;
    scalarD tt146;
    scalarD tt147;
    scalarD tt148;
    scalarD tt149;
    scalarD tt150;
    scalarD tt151;
    scalarD tt152;
    scalarD tt153;
    scalarD tt154;
    scalarD tt155;
    scalarD tt156;
    scalarD tt157;
    scalarD tt158;
    scalarD tt159;
    scalarD tt160;
    scalarD tt161;
    scalarD tt162;
    scalarD tt163;
    scalarD tt164;
    scalarD tt165;
    scalarD tt166;
    scalarD tt167;
    scalarD tt168;
    scalarD tt169;
    scalarD tt170;
    scalarD tt171;
    scalarD tt172;
    scalarD tt173;
    scalarD tt174;
    scalarD tt175;
    scalarD tt176;
    scalarD tt177;
    scalarD tt178;
    scalarD tt179;
    scalarD tt180;
    scalarD tt181;
    scalarD tt182;
    scalarD tt183;
    scalarD tt184;
    scalarD tt185;
    scalarD tt186;
    scalarD tt187;
    scalarD tt188;
    scalarD tt189;
    scalarD tt190;
    scalarD tt191;
    scalarD tt192;
    scalarD tt193;
    scalarD tt194;
    scalarD tt195;
    scalarD tt196;
    scalarD tt197;
    scalarD tt198;
    scalarD tt199;
    scalarD tt200;

    tt1=F31+F13;
    tt2=F13-F31;
    tt3=F21-F12;
    tt4=pow(tt3,2);
    tt5=pow(tt2,2);
    tt6=F32-F23;
    tt7=pow(tt6,2);
    tt8=std::max<scalarD>(tt7/4.0+tt5/4.0+tt4/4.0,1E-12);
    tt9=sqrt(tt8);
    tt10=1/tt9;
    tt11=0.5*tt9;
    tt12=cos(tt11);
    tt13=sin(tt11);
    tt14=tt2*tt10*tt12*tt13/2.0;
    tt15=1/tt8;
    tt16=pow(tt13,2);
    tt17=tt3*tt6*tt15*tt16/4.0;
    tt18=tt17+tt14;
    tt19=F21+F12;
    tt20=-tt3*tt10*tt12*tt13/2.0;
    tt21=tt2*tt6*tt15*tt16/4.0;
    tt22=tt21+tt20;
    tt23=F11+1;
    tt24=pow(tt12,2);
    tt25=-tt4*tt15*tt16/4.0;
    tt26=-tt5*tt15*tt16/4.0;
    tt27=tt7*tt15*tt16/4.0+tt26+tt25+tt24;
    tt28=F32+F23;
    tt29=F22+1;
    tt30=F33+1;
    tt31=-tt6*tt10*tt12*tt13/2.0;
    tt32=tt3*tt2*tt15*tt16/4.0;
    tt33=tt32+tt31;
    tt34=tt3*tt10*tt12*tt13/2.0;
    tt35=tt21+tt34;
    tt36=-tt7*tt15*tt16/4.0;
    tt37=tt36+tt5*tt15*tt16/4.0+tt25+tt24;
    tt38=tt6*tt10*tt12*tt13/2.0;
    tt39=tt32+tt38;
    tt40=-tt2*tt10*tt12*tt13/2.0;
    tt41=tt17+tt40;
    tt42=tt36+tt26+tt4*tt15*tt16/4.0+tt24;
    RSF(0,0)=tt23*tt27+tt19*tt22+tt1*tt18-1;
    RSF(0,1)=tt19*tt27/2.0+2*tt29*tt22+tt28*tt18;
    RSF(0,2)=2*tt18*tt30+tt1*tt27/2.0+tt28*tt22;
    RSF(1,0)=tt19*tt37/2.0+2*tt23*tt35+tt1*tt33;
    RSF(1,1)=tt29*tt37+tt19*tt35+tt28*tt33-1;
    RSF(1,2)=2*tt33*tt30+tt28*tt37/2.0+tt1*tt35;
    RSF(2,0)=tt1*tt42/2.0+2*tt23*tt41+tt19*tt39;
    RSF(2,1)=tt28*tt42/2.0+tt19*tt41+2*tt29*tt39;
    RSF(2,2)=tt42*tt30+tt1*tt41+tt28*tt39-1;
    if(!gradRSF)return RSF;

    tt43=pow(tt3,3);
    tt44=1/pow(tt9,3);
    tt45=-0.0625*tt43*tt44*tt12*tt13;
    tt46=-0.0625*tt3*tt5*tt44*tt12*tt13;
    tt47=0.0625*tt3*tt7*tt44*tt12*tt13;
    tt48=-0.25*tt3*tt10*tt12*tt13;
    tt49=1/pow(tt8,2);
    tt50=tt43*tt49*tt16/8.0;
    tt51=tt3*tt5*tt49*tt16/8.0;
    tt52=-tt3*tt7*tt49*tt16/8.0;
    tt53=-tt3*tt15*tt16/2.0;
    tt54=tt53+tt52+tt51+tt50+tt48+tt47+tt46+tt45;
    tt55=-0.0625*tt4*tt15*tt24;
    tt56=tt4*tt44*tt12*tt13/8.0;
    tt57=0.0625*tt3*tt2*tt6*tt44*tt12*tt13;
    tt58=-tt10*tt12*tt13/2.0;
    tt59=-tt3*tt2*tt6*tt49*tt16/8.0;
    tt60=0.0625*tt4*tt15*tt16;
    tt61=tt60+tt59+tt58+tt57+tt56+tt55;
    tt62=0.0625*tt3*tt2*tt15*tt24;
    tt63=-tt3*tt2*tt44*tt12*tt13/8.0;
    tt64=0.0625*tt4*tt6*tt44*tt12*tt13;
    tt65=-tt4*tt6*tt49*tt16/8.0;
    tt66=-0.0625*tt3*tt2*tt15*tt16;
    tt67=tt6*tt15*tt16/4.0;
    tt68=tt67+tt66+tt65+tt64+tt63+tt62;
    tt69=-0.0625*tt5*tt15*tt24;
    tt70=tt5*tt44*tt12*tt13/8.0;
    tt71=-0.0625*tt3*tt2*tt6*tt44*tt12*tt13;
    tt72=tt3*tt2*tt6*tt49*tt16/8.0;
    tt73=0.0625*tt5*tt15*tt16;
    tt74=tt73+tt72+tt58+tt71+tt70+tt69;
    tt75=pow(tt2,3);
    tt76=0.0625*tt75*tt44*tt12*tt13;
    tt77=0.0625*tt4*tt2*tt44*tt12*tt13;
    tt78=-0.0625*tt2*tt7*tt44*tt12*tt13;
    tt79=0.25*tt2*tt10*tt12*tt13;
    tt80=-tt75*tt49*tt16/8.0;
    tt81=-tt4*tt2*tt49*tt16/8.0;
    tt82=tt2*tt7*tt49*tt16/8.0;
    tt83=tt2*tt15*tt16/2.0;
    tt84=tt83+tt82+tt81+tt80+tt79+tt78+tt77+tt76;
    tt85=-0.0625*tt5*tt6*tt44*tt12*tt13;
    tt86=tt5*tt6*tt49*tt16/8.0;
    tt87=-tt6*tt15*tt16/4.0;
    tt88=tt87+tt66+tt86+tt85+tt63+tt62;
    tt89=0.0625*tt43*tt44*tt12*tt13;
    tt90=0.0625*tt3*tt5*tt44*tt12*tt13;
    tt91=-0.0625*tt3*tt7*tt44*tt12*tt13;
    tt92=0.25*tt3*tt10*tt12*tt13;
    tt93=-tt43*tt49*tt16/8.0;
    tt94=-tt3*tt5*tt49*tt16/8.0;
    tt95=tt3*tt7*tt49*tt16/8.0;
    tt96=tt3*tt15*tt16/2.0;
    tt97=tt96+tt95+tt94+tt93+tt92+tt91+tt90+tt89;
    tt98=0.0625*tt4*tt15*tt24;
    tt99=-tt4*tt44*tt12*tt13/8.0;
    tt100=tt10*tt12*tt13/2.0;
    tt101=-0.0625*tt4*tt15*tt16;
    tt102=tt101+tt72+tt100+tt71+tt99+tt98;
    tt103=-0.0625*tt3*tt2*tt15*tt24;
    tt104=tt3*tt2*tt44*tt12*tt13/8.0;
    tt105=-0.0625*tt4*tt6*tt44*tt12*tt13;
    tt106=tt4*tt6*tt49*tt16/8.0;
    tt107=0.0625*tt3*tt2*tt15*tt16;
    tt108=tt87+tt107+tt106+tt105+tt104+tt103;
    tt109=pow(tt6,3);
    tt110=0.0625*tt109*tt44*tt12*tt13;
    tt111=-0.25*tt6*tt10*tt12*tt13;
    tt112=-tt109*tt49*tt16/8.0;
    tt113=tt6*tt15*tt16/2.0;
    tt114=tt113+tt112+tt86+tt106+tt111+tt110+tt85+tt105;
    tt115=-0.0625*tt3*tt6*tt15*tt24;
    tt116=tt3*tt6*tt44*tt12*tt13/8.0;
    tt117=0.0625*tt2*tt7*tt44*tt12*tt13;
    tt118=-tt2*tt7*tt49*tt16/8.0;
    tt119=tt2*tt15*tt16/4.0;
    tt120=0.0625*tt3*tt6*tt15*tt16;
    tt121=tt120+tt119+tt118+tt117+tt116+tt115;
    tt122=0.0625*tt2*tt6*tt15*tt24;
    tt123=-tt2*tt6*tt44*tt12*tt13/8.0;
    tt124=tt3*tt15*tt16/4.0;
    tt125=-0.0625*tt2*tt6*tt15*tt16;
    tt126=tt125+tt124+tt52+tt47+tt123+tt122;
    tt127=0.0625*tt5*tt15*tt24;
    tt128=-tt5*tt44*tt12*tt13/8.0;
    tt129=-0.0625*tt5*tt15*tt16;
    tt130=tt129+tt59+tt100+tt57+tt128+tt127;
    tt131=-0.0625*tt75*tt44*tt12*tt13;
    tt132=-0.0625*tt4*tt2*tt44*tt12*tt13;
    tt133=-0.25*tt2*tt10*tt12*tt13;
    tt134=tt75*tt49*tt16/8.0;
    tt135=tt4*tt2*tt49*tt16/8.0;
    tt136=-tt2*tt15*tt16/2.0;
    tt137=tt136+tt118+tt135+tt134+tt133+tt117+tt132+tt131;
    tt138=0.0625*tt5*tt6*tt44*tt12*tt13;
    tt139=-tt5*tt6*tt49*tt16/8.0;
    tt140=tt67+tt107+tt139+tt138+tt104+tt103;
    tt141=-0.0625*tt109*tt44*tt12*tt13;
    tt142=0.25*tt6*tt10*tt12*tt13;
    tt143=tt109*tt49*tt16/8.0;
    tt144=-tt6*tt15*tt16/2.0;
    tt145=tt144+tt143+tt139+tt65+tt142+tt141+tt138+tt64;
    tt146=0.0625*tt3*tt6*tt15*tt24;
    tt147=-tt3*tt6*tt44*tt12*tt13/8.0;
    tt148=-tt2*tt15*tt16/4.0;
    tt149=-0.0625*tt3*tt6*tt15*tt16;
    tt150=tt149+tt148+tt82+tt78+tt147+tt146;
    tt151=-0.0625*tt2*tt6*tt15*tt24;
    tt152=tt2*tt6*tt44*tt12*tt13/8.0;
    tt153=-tt3*tt15*tt16/4.0;
    tt154=0.0625*tt2*tt6*tt15*tt16;
    tt155=tt154+tt153+tt95+tt91+tt152+tt151;
    tt156=tt53+tt95+tt94+tt50+tt48+tt91+tt90+tt45;
    tt157=tt101+tt59+tt100+tt57+tt99+tt98;
    tt158=tt120+tt119+tt81+tt116+tt77+tt115;
    tt159=tt37/2.0;
    tt160=tt136+tt118+tt81+tt134+tt79+tt117+tt77+tt131;
    tt161=tt87+tt107+tt86+tt85+tt104+tt103;
    tt162=tt125+tt153+tt51+tt123+tt46+tt122;
    tt163=tt96+tt52+tt51+tt93+tt92+tt47+tt46+tt89;
    tt164=tt60+tt72+tt58+tt71+tt56+tt55;
    tt165=tt149+tt148+tt135+tt147+tt132+tt146;
    tt166=tt144+tt143+tt139+tt106+tt111+tt141+tt138+tt105;
    tt167=tt149+tt119+tt118+tt117+tt147+tt146;
    tt168=-0.0625*tt7*tt15*tt24;
    tt169=tt7*tt44*tt12*tt13/8.0;
    tt170=0.0625*tt7*tt15*tt16;
    tt171=tt170+tt59+tt58+tt169+tt57+tt168;
    tt172=tt83+tt82+tt135+tt80+tt133+tt78+tt132+tt76;
    tt173=tt67+tt66+tt139+tt138+tt63+tt62;
    tt174=tt154+tt124+tt94+tt152+tt90+tt151;
    tt175=tt113+tt112+tt86+tt65+tt142+tt110+tt85+tt64;
    tt176=tt120+tt148+tt82+tt78+tt116+tt115;
    tt177=0.0625*tt7*tt15*tt24;
    tt178=-tt7*tt44*tt12*tt13/8.0;
    tt179=-0.0625*tt7*tt15*tt16;
    tt180=tt179+tt72+tt100+tt178+tt71+tt177;
    tt181=tt96+tt95+tt51+tt93+tt48+tt91+tt46+tt89;
    tt182=tt67+tt107+tt65+tt64+tt104+tt103;
    tt183=tt149+tt119+tt81+tt147+tt77+tt146;
    tt184=tt129+tt72+tt100+tt71+tt128+tt127;
    tt185=tt83+tt118+tt135+tt80+tt79+tt117+tt132+tt76;
    tt186=tt154+tt153+tt51+tt152+tt46+tt151;
    tt187=tt42/2.0;
    tt188=tt53+tt52+tt94+tt50+tt92+tt47+tt90+tt45;
    tt189=tt87+tt66+tt106+tt105+tt63+tt62;
    tt190=tt120+tt148+tt135+tt116+tt132+tt115;
    tt191=tt144+tt143+tt86+tt65+tt111+tt141+tt85+tt64;
    tt192=tt154+tt124+tt52+tt47+tt152+tt151;
    tt193=tt179+tt59+tt100+tt178+tt57+tt177;
    tt194=tt73+tt59+tt58+tt57+tt70+tt69;
    tt195=tt136+tt82+tt81+tt134+tt133+tt78+tt77+tt131;
    tt196=tt125+tt124+tt94+tt123+tt90+tt122;
    tt197=tt113+tt112+tt139+tt106+tt142+tt110+tt138+tt105;
    tt198=tt125+tt153+tt95+tt91+tt123+tt122;
    tt199=tt170+tt72+tt58+tt169+tt71+tt168;
    tt200=tt27/2.0;
    (*gradRSF)(0,0)=tt27;
    (*gradRSF)(0,1)=tt1*tt68+tt19*tt61+tt23*tt54+tt21+tt20;
    (*gradRSF)(0,2)=tt19*tt88+tt23*tt84+tt1*tt74+tt17+tt14;
    (*gradRSF)(0,3)=tt1*tt108+tt19*tt102+tt23*tt97+tt21+tt20;
    (*gradRSF)(0,4)=0;
    (*gradRSF)(0,5)=tt1*tt126+tt19*tt121+tt23*tt114;
    (*gradRSF)(0,6)=tt19*tt140+tt23*tt137+tt1*tt130+tt17+tt14;
    (*gradRSF)(0,7)=tt1*tt155+tt19*tt150+tt23*tt145;
    (*gradRSF)(0,8)=0;
    (*gradRSF)(1,0)=2*tt35;
    (*gradRSF)(1,1)=tt159+tt1*tt158+2*tt23*tt157+tt19*tt156/2.0;
    (*gradRSF)(1,2)=tt1*tt162+2*tt23*tt161+tt19*tt160/2.0+tt32+tt31;
    (*gradRSF)(1,3)=tt159+tt1*tt165+2*tt23*tt164+tt19*tt163/2.0;
    (*gradRSF)(1,4)=0;
    (*gradRSF)(1,5)=tt1*tt171+2*tt23*tt167+tt19*tt166/2.0;
    (*gradRSF)(1,6)=tt1*tt174+2*tt23*tt173+tt19*tt172/2.0+tt32+tt31;
    (*gradRSF)(1,7)=tt1*tt180+2*tt23*tt176+tt19*tt175/2.0;
    (*gradRSF)(1,8)=0;
    (*gradRSF)(2,0)=2*tt41;
    (*gradRSF)(2,1)=tt19*tt183+2*tt23*tt182+tt1*tt181/2.0+tt32+tt38;
    (*gradRSF)(2,2)=tt187+tt19*tt186+tt1*tt185/2.0+2*tt23*tt184;
    (*gradRSF)(2,3)=tt19*tt190+2*tt23*tt189+tt1*tt188/2.0+tt32+tt38;
    (*gradRSF)(2,4)=0;
    (*gradRSF)(2,5)=tt19*tt193+2*tt23*tt192+tt1*tt191/2.0;
    (*gradRSF)(2,6)=tt187+tt19*tt196+tt1*tt195/2.0+2*tt23*tt194;
    (*gradRSF)(2,7)=tt19*tt199+2*tt23*tt198+tt1*tt197/2.0;
    (*gradRSF)(2,8)=0;
    (*gradRSF)(3,0)=0;
    (*gradRSF)(3,1)=tt200+tt28*tt68+2*tt29*tt61+tt19*tt54/2.0;
    (*gradRSF)(3,2)=2*tt29*tt88+tt19*tt84/2.0+tt28*tt74;
    (*gradRSF)(3,3)=tt200+tt28*tt108+2*tt29*tt102+tt19*tt97/2.0;
    (*gradRSF)(3,4)=2*tt22;
    (*gradRSF)(3,5)=tt28*tt126+2*tt29*tt121+tt19*tt114/2.0+tt17+tt14;
    (*gradRSF)(3,6)=2*tt29*tt140+tt19*tt137/2.0+tt28*tt130;
    (*gradRSF)(3,7)=tt28*tt155+2*tt29*tt150+tt19*tt145/2.0+tt17+tt14;
    (*gradRSF)(3,8)=0;
    (*gradRSF)(4,0)=0;
    (*gradRSF)(4,1)=tt28*tt158+tt19*tt157+tt29*tt156+tt21+tt34;
    (*gradRSF)(4,2)=tt28*tt162+tt19*tt161+tt29*tt160;
    (*gradRSF)(4,3)=tt28*tt165+tt19*tt164+tt29*tt163+tt21+tt34;
    (*gradRSF)(4,4)=tt37;
    (*gradRSF)(4,5)=tt28*tt171+tt19*tt167+tt29*tt166+tt32+tt31;
    (*gradRSF)(4,6)=tt28*tt174+tt19*tt173+tt29*tt172;
    (*gradRSF)(4,7)=tt28*tt180+tt19*tt176+tt29*tt175+tt32+tt31;
    (*gradRSF)(4,8)=0;
    (*gradRSF)(5,0)=0;
    (*gradRSF)(5,1)=2*tt29*tt183+tt19*tt182+tt28*tt181/2.0+tt17+tt40;
    (*gradRSF)(5,2)=2*tt29*tt186+tt28*tt185/2.0+tt19*tt184;
    (*gradRSF)(5,3)=2*tt29*tt190+tt19*tt189+tt28*tt188/2.0+tt17+tt40;
    (*gradRSF)(5,4)=2*tt39;
    (*gradRSF)(5,5)=2*tt29*tt193+tt187+tt19*tt192+tt28*tt191/2.0;
    (*gradRSF)(5,6)=2*tt29*tt196+tt28*tt195/2.0+tt19*tt194;
    (*gradRSF)(5,7)=2*tt29*tt199+tt187+tt19*tt198+tt28*tt197/2.0;
    (*gradRSF)(5,8)=0;
    (*gradRSF)(6,0)=0;
    (*gradRSF)(6,1)=2*tt68*tt30+tt28*tt61+tt1*tt54/2.0;
    (*gradRSF)(6,2)=2*tt74*tt30+tt200+tt28*tt88+tt1*tt84/2.0;
    (*gradRSF)(6,3)=2*tt108*tt30+tt28*tt102+tt1*tt97/2.0;
    (*gradRSF)(6,4)=0;
    (*gradRSF)(6,5)=2*tt126*tt30+tt28*tt121+tt1*tt114/2.0+tt21+tt20;
    (*gradRSF)(6,6)=2*tt130*tt30+tt200+tt28*tt140+tt1*tt137/2.0;
    (*gradRSF)(6,7)=2*tt155*tt30+tt28*tt150+tt1*tt145/2.0+tt21+tt20;
    (*gradRSF)(6,8)=2*tt18;
    (*gradRSF)(7,0)=0;
    (*gradRSF)(7,1)=2*tt158*tt30+tt1*tt157+tt28*tt156/2.0;
    (*gradRSF)(7,2)=2*tt162*tt30+tt1*tt161+tt28*tt160/2.0+tt21+tt34;
    (*gradRSF)(7,3)=2*tt165*tt30+tt1*tt164+tt28*tt163/2.0;
    (*gradRSF)(7,4)=0;
    (*gradRSF)(7,5)=2*tt171*tt30+tt159+tt1*tt167+tt28*tt166/2.0;
    (*gradRSF)(7,6)=2*tt174*tt30+tt1*tt173+tt28*tt172/2.0+tt21+tt34;
    (*gradRSF)(7,7)=2*tt180*tt30+tt159+tt1*tt176+tt28*tt175/2.0;
    (*gradRSF)(7,8)=2*tt33;
    (*gradRSF)(8,0)=0;
    (*gradRSF)(8,1)=tt181*tt30+tt28*tt183+tt1*tt182;
    (*gradRSF)(8,2)=tt185*tt30+tt28*tt186+tt1*tt184+tt17+tt40;
    (*gradRSF)(8,3)=tt188*tt30+tt28*tt190+tt1*tt189;
    (*gradRSF)(8,4)=0;
    (*gradRSF)(8,5)=tt191*tt30+tt28*tt193+tt1*tt192+tt32+tt38;
    (*gradRSF)(8,6)=tt195*tt30+tt28*tt196+tt1*tt194+tt17+tt40;
    (*gradRSF)(8,7)=tt197*tt30+tt28*tt199+tt1*tt198+tt32+tt38;
    (*gradRSF)(8,8)=tt42;
    return RSF;
}
static FORCE_INLINE Mat3d calcFGradDU(const Mat3d& U)
{
    Mat3d dexp;
    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;

    tt1=U(1,0)-U(0,1);
    tt2=(U(1,0)+U(0,1))/2.0;
    tt3=U(0,2)-U(2,0);
    tt4=(U(2,0)+U(0,2))/2.0;
    tt5=U(2,1)-U(1,2);
    tt6=(U(2,1)+U(1,2))/2.0;
    dexp(0,0)=U(0,0);
    dexp(0,1)=tt2-tt1/2.0;
    dexp(0,2)=tt4+tt3/2.0;
    dexp(1,0)=tt2+tt1/2.0;
    dexp(1,1)=U(1,1);
    dexp(1,2)=tt6-tt5/2.0;
    dexp(2,0)=tt4-tt3/2.0;
    dexp(2,1)=tt6+tt5/2.0;
    dexp(2,2)=U(2,2);
    return dexp;
}
static FORCE_INLINE Mat3d calcFGradDDU(const Mat3d& U)
{
    Mat3d dexp;
    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;

    tt1=U(1,0)-U(0,1);
    tt2=U(1,0)+U(0,1);
    tt3=-pow(tt1,2)/2.0;
    tt4=U(0,2)-U(2,0);
    tt5=-pow(tt4,2)/2.0;
    tt6=U(2,0)+U(0,2);
    tt7=U(2,1)-U(1,2);
    tt8=tt4*tt7/4.0;
    tt9=U(2,1)+U(1,2);
    tt10=tt1*tt7/4.0;
    tt11=-pow(tt7,2)/2.0;
    tt12=tt1*tt4/4.0;
    dexp(0,0)=tt4*tt6/2.0+(tt5+tt3)/2.0-tt1*tt2/2.0;
    dexp(0,1)=tt4*tt9/2.0+tt8-tt1*U(1,1);
    dexp(0,2)=tt4*U(2,2)-tt1*tt9/2.0+tt10;
    dexp(1,0)=-tt6*tt7/2.0+tt8+U(0,0)*tt1;
    dexp(1,1)=(tt11+tt3)/2.0-tt7*tt9/2.0+tt1*tt2/2.0;
    dexp(1,2)=-tt7*U(2,2)+tt1*tt6/2.0+tt12;
    dexp(2,0)=tt2*tt7/2.0+tt10-U(0,0)*tt4;
    dexp(2,1)=U(1,1)*tt7-tt2*tt4/2.0+tt12;
    dexp(2,2)=(tt11+tt5)/2.0+tt7*tt9/2.0-tt4*tt6/2.0;
    return dexp;
}
static FORCE_INLINE Mat3d calcFGradDUDV(const Mat3d& U,const Mat3d& V)
{
    Mat3d dexp;
    //temp
    scalarD tt1;
    scalarD tt2;
    scalarD tt3;
    scalarD tt4;
    scalarD tt5;
    scalarD tt6;
    scalarD tt7;
    scalarD tt8;
    scalarD tt9;
    scalarD tt10;
    scalarD tt11;
    scalarD tt12;
    scalarD tt13;
    scalarD tt14;
    scalarD tt15;
    scalarD tt16;
    scalarD tt17;
    scalarD tt18;
    scalarD tt19;
    scalarD tt20;
    scalarD tt21;

    tt1=U(1,0)+U(0,1);
    tt2=V(1,0)-V(0,1);
    tt3=U(1,0)-U(0,1);
    tt4=V(1,0)+V(0,1);
    tt5=-tt3*tt2/2.0;
    tt6=U(0,2)-U(2,0);
    tt7=V(0,2)-V(2,0);
    tt8=-tt6*tt7/2.0;
    tt9=U(2,0)+U(0,2);
    tt10=V(2,0)+V(0,2);
    tt11=U(2,1)-U(1,2);
    tt12=tt7*tt11/8.0;
    tt13=U(2,1)+U(1,2);
    tt14=V(2,1)-V(1,2);
    tt15=tt6*tt14/8.0;
    tt16=V(2,1)+V(1,2);
    tt17=tt2*tt11/8.0;
    tt18=tt3*tt14/8.0;
    tt19=-tt11*tt14/2.0;
    tt20=tt2*tt6/8.0;
    tt21=tt3*tt7/8.0;
    dexp(0,0)=tt6*tt10/4.0+tt9*tt7/4.0+(tt8+tt5)/2.0-tt3*tt4/4.0-tt1*tt2/4.0;
    dexp(0,1)=tt6*tt16/4.0+tt15+tt7*tt13/4.0+tt12-tt3*V(1,1)/2.0-tt2*U(1,1)/2.0;
    dexp(0,2)=tt6*V(2,2)/2.0+tt7*U(2,2)/2.0-tt3*tt16/4.0+tt18-tt2*tt13/4.0+tt17;
    dexp(1,0)=-tt9*tt14/4.0+tt15-tt10*tt11/4.0+tt12+U(0,0)*tt2/2.0+V(0,0)*tt3/2.0;
    dexp(1,1)=(tt19+tt5)/2.0-tt11*tt16/4.0-tt13*tt14/4.0+tt3*tt4/4.0+tt1*tt2/4.0;
    dexp(1,2)=-tt11*V(2,2)/2.0-tt14*U(2,2)/2.0+tt3*tt10/4.0+tt21+tt2*tt9/4.0+tt20;
    dexp(2,0)=tt1*tt14/4.0+tt18+tt4*tt11/4.0+tt17-U(0,0)*tt7/2.0-V(0,0)*tt6/2.0;
    dexp(2,1)=U(1,1)*tt14/2.0+V(1,1)*tt11/2.0-tt1*tt7/4.0+tt21-tt4*tt6/4.0+tt20;
    dexp(2,2)=(tt19+tt8)/2.0+tt11*tt16/4.0+tt13*tt14/4.0-tt6*tt10/4.0-tt9*tt7/4.0;
    return dexp;
}
static FORCE_INLINE void debugCalcFGradD()
{
    Mat3d F=Mat3d::Zero();
    Eigen::Matrix<scalarD,9,9> gradRSF,gradRSF1,gradRSF2;
    calcFGrad(F,&gradRSF);

    Mat3d U=Mat3d::Random();
    Mat3d V=Mat3d::Random();
    Eigen::Map<Cold> mapU(U.data(),9);
    Eigen::Map<Cold> mapV(V.data(),9);

    //debug dexpU
    INFOV("DFGradDU: %f %f",
          (gradRSF*mapU).norm(),
          calcFGradDU(U).norm());
    //debug ddexpU
    calcFGrad(U*1E-7f,&gradRSF1);
    INFOV("DDFGradDDU: %f %f",
          ((gradRSF1-gradRSF)*mapU/1E-7f).norm(),
          calcFGradDDU(U).norm());
    //debug ddexpUV
    calcFGrad(V*1E-7f,&gradRSF2);
    INFOV("DDFGradDUDV: %f %f %f",
          ((gradRSF1-gradRSF)*mapV/1E-7f).norm(),
          ((gradRSF2-gradRSF)*mapU/1E-7f).norm(),
          calcFGradDUDV(U,V).norm());
}
static FORCE_INLINE void adjustF3D(Mat3& F,Mat3& U,Mat3& V,Vec3& S)
{
    Eigen::SelfAdjointEigenSolver<Mat3> eig(F.transpose()*F);
    //F Hat
    S=eig.eigenvalues().cwiseAbs().cwiseSqrt();
    //adjust V
    V=eig.eigenvectors();
    if(V.determinant() < 0.0f)
        V.col(0)*=-1.0f;
    //adjust U
#define F_EPS 1E-6f
    sizeType good[3];
    sizeType nrGood=0;
    for(sizeType v=0; v<3; v++)
        if(S[v] > F_EPS) {
            U.col(v)=F*V.col(v)/S[v];
            good[nrGood++]=v;
        }
    //the tet is a point
    if(nrGood == 0)
        U.setIdentity();
    //set columns of U to be orthogonal
    else for(sizeType v=0; v<3; v++)
            if(S[v] <= F_EPS) {
                if(nrGood == 1) {
                    Vec3 d;
                    scalar nd;
                    while(true) {
                        d=U.col(good[0]).cross(Vec3::Random());
                        nd=d.norm();
                        if(nd > 1E-3f) {
                            U.col(v)=d/nd;
                            break;
                        }
                    }
                } else {
                    ASSERT(nrGood == 2)
                    U.col(v)=U.col(good[0]).cross(U.col(good[1]));
                }
                good[nrGood++]=v;
            }
#undef F_EPS
    //negate negative U columns
    if(U.determinant() < 0.0f) {
        sizeType id;
        S.minCoeff(&id);
        S[id]*=-1.0f;
        U.col(id)*=-1.0f;
    }
    //rebuild F
    Mat3 Sigma=Mat3::Zero();
    Sigma.diagonal()=S;
    F=U*Sigma*V.transpose();
}
static FORCE_INLINE void adjustF2D(Mat3& F,Mat3& U,Mat3& V,Vec3& S)
{
    Eigen::SelfAdjointEigenSolver<Mat2> eig((F.transpose()*F).block<2,2>(0,0));
    //F Hat
    S.setZero();
    S.block<2,1>(0,0)=eig.eigenvalues().cwiseAbs().cwiseSqrt();
    //adjust V
    V.setZero();
    V.block<2,2>(0,0)=eig.eigenvectors();
    if(V.block<2,2>(0,0).determinant() < 0.0f)
        V.col(0)*=-1.0f;
    //adjust U
#define F_EPS 1E-6f
    U.setZero();
    sizeType good[2];
    sizeType nrGood=0;
    for(sizeType v=0; v<2; v++)
        if(S[v] > F_EPS) {
            U.col(v)=F*V.col(v)/S[v];
            good[nrGood++]=v;
        }
    //the tet is a point
    if(nrGood == 0)
        U.block<2,2>(0,0).setIdentity();
    //set columns of U to be orthogonal
    else for(sizeType v=0; v<2; v++)
            if(S[v] <= F_EPS) {
                U.col(v)=Vec3(-U(1,good[0]),U(0,good[0]),0.0f);
                good[nrGood++]=v;
            }
#undef F_EPS
    //negate negative U columns
    if(U.block<2,2>(0,0).determinant() < 0.0f) {
        sizeType id;
        S.block<2,1>(0,0).minCoeff(&id);
        S[id]*=-1.0f;
        U.col(id)*=-1.0f;
    }
    //rebuild F
    Mat3 Sigma=Mat3::Zero();
    Sigma.diagonal()=S;
    F=U*Sigma*V.transpose();
}
static FORCE_INLINE void derivSVD(Mat3& derivU,Vec3& derivD,Mat3& derivV,const Mat3& LHS,const Mat3& U,const Vec3& D,const Mat3& V)
{
#define SOLVE_SAFE(a,b,c)				\
tmp << D(c),D(b),D(b),D(c);				\
if(std::abs(tmp(0,0)-tmp(0,1)) < 1E-6f)		\
	tmpInv=(tmp*tmp+reg).inverse()*tmp;	\
else									\
	tmpInv=tmp.inverse();				\
a=tmpInv*Vec2(LHS(b,c),-LHS(c,b));

    derivD=Vec3(LHS(0,0),LHS(1,1),LHS(2,2));

    //21,31,32
    Mat2 reg=Mat2::Identity()*1E-6f,tmp,tmpInv;
    Vec2 omg10,omg20,omg21;
    SOLVE_SAFE(omg10,1,0)
    SOLVE_SAFE(omg20,2,0)
    SOLVE_SAFE(omg21,2,1)
#undef SOLVE_SAFE

    Mat3 omgU;
    omgU<<       0,-omg10(0),-omg20(0),
         omg10(0),        0,-omg21(0),
         omg20(0), omg21(0),        0;
    derivU=U*omgU;

    Mat3 omgV;
    omgV<<       0,-omg10(1),-omg20(1),
         omg10(1),        0,-omg21(1),
         omg20(1), omg21(1),        0;
    derivV=-V*omgV;
}
static FORCE_INLINE void derivRDF(Eigen::Matrix<scalarD,9,9>& DRDF,const Mat3& F,const Mat3& U,const Vec3& D,const Mat3& V)
{
    Vec3 derivD;
    Mat3 derivU,derivV;
    if(D.minCoeff() < 1E-3f)
        DRDF.setZero();
    for(sizeType r=0; r<3; r++)
        for(sizeType c=0; c<3; c++) {
            derivSVD(derivU,derivD,derivV,U.row(r).transpose()*V.row(c),U,D,V);
            Eigen::Map<Mat3d>(DRDF.col(c*3+r).data())=(derivU*V.transpose()+U*derivV.transpose()).cast<scalarD>();
        }
}

PRJ_END

#endif
