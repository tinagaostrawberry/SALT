function [all_stampVal,I_v] = stampOfIbrGFMMatrices_TD_flat( ...
    obj,all_stampVal,I_v,x_v,x_p)
% GRID-FORMING IBR stamp - flat Nx1 column form, vectorised over all IBRs.
% Terms are accumulated in a FIXED order that matches MATLAB's BLAS
% reduction order, so results are reproducible to the last bit:
%   matrix x column : left-to-right over columns
%   row x column    : pairwise tree (odd K: t1+(t2+t3) then +(t4+t5), ...)
% Do not reorder them; the dropped terms are structural zeros.
% The 2x2 rotations stamp_jw0=[0 -1/w0;1/w0 0] and stamp_j=[0 -1;1 0] are
% expanded inline.

N = obj.numIbrGFM;
RW  = obj.getRowIdx_fromIbrGFMIndices((1:N).')-1;
PS  = obj.getStampIdx_fromIbrGFMIndices((1:N).')-1 + (1:obj.ibrGFM_numStampVals);
ROW = RW + (1:obj.numIbrGFMEqns);
[qIabc,qVg,qVo,qIf,qIo,qVi,qGam,qIfR,qX,qWr,qVoQd,qPQ,qTh,qVRa,qVR,qVoR] = ...
    getIbrGFMIndicesTD;
VA = obj.getRowIdxAbc_fromBusNums(obj.ibrGFM_bus(:).');
[~,BUS] = ismembertol(obj.ibrGFM_bus,obj.bus);

% ---- states -------------------------------------------------------------
ia=x_v(RW+qIabc(1)); ib=x_v(RW+qIabc(2)); ic=x_v(RW+qIabc(3));
vgd=x_v(RW+qVg(1));  vgq=x_v(RW+qVg(2));
vod=x_v(RW+qVo(1));  voq=x_v(RW+qVo(2));
ifd=x_v(RW+qIf(1));  ifq=x_v(RW+qIf(2));
iod=x_v(RW+qIo(1));  ioq=x_v(RW+qIo(2));
vid=x_v(RW+qVi(1));  viq=x_v(RW+qVi(2));
gd =x_v(RW+qGam(1)); gq =x_v(RW+qGam(2));
frd=x_v(RW+qIfR(1)); frq=x_v(RW+qIfR(2));
xd =x_v(RW+qX(1));   xq =x_v(RW+qX(2));
wr =x_v(RW+qWr);
vrd=x_v(RW+qVoR(1)); vrq=x_v(RW+qVoR(2));
pav=x_v(RW+qPQ(1));  qav=x_v(RW+qPQ(2));
th =x_v(RW+qTh);
va =x_v(VA(1,:).');  vb =x_v(VA(2,:).');  vc =x_v(VA(3,:).');
ad =x_v(RW+qVRa(1)); aq =x_v(RW+qVRa(2));
rd =x_v(RW+qVR(1));  rq =x_v(RW+qVR(2));
vqd=x_v(RW+qVoQd);
%
vgdp=x_p(RW+qVg(1)); vgqp=x_p(RW+qVg(2));
vodp=x_p(RW+qVo(1)); voqp=x_p(RW+qVo(2));
ifdp=x_p(RW+qIf(1)); ifqp=x_p(RW+qIf(2));
iodp=x_p(RW+qIo(1)); ioqp=x_p(RW+qIo(2));
vidp=x_p(RW+qVi(1)); viqp=x_p(RW+qVi(2));
gdp =x_p(RW+qGam(1));gqp =x_p(RW+qGam(2));
frdp=x_p(RW+qIfR(1));frqp=x_p(RW+qIfR(2));
xdp =x_p(RW+qX(1));  xqp =x_p(RW+qX(2));
wrp =x_p(RW+qWr);
vrdp=x_p(RW+qVoR(1));vrqp=x_p(RW+qVoR(2));
pavp=x_p(RW+qPQ(1)); qavp=x_p(RW+qPQ(2));
thp =x_p(RW+qTh);
adp =x_p(RW+qVRa(1));aqp =x_p(RW+qVRa(2));

% ---- constants ----------------------------------------------------------
Lf=obj.ibrGFM_Lf(:);   rf=obj.ibrGFM_Rf(:);   Cf=obj.ibrGFM_Cf(:);
Rcp=obj.ibrGFM_Rcap(:);Lc=obj.ibrGFM_Lc(:);   rc=obj.ibrGFM_Rc(:);
kci=obj.ibrGFM_kC_i(:);kcp=obj.ibrGFM_kC_p(:);Gc=obj.ibrGFM_GC(:);
kvi=obj.ibrGFM_kV_i(:);kvp=obj.ibrGFM_kV_p(:);Gv=obj.ibrGFM_GV(:);
wmR=obj.ibrGFM_wmeas(:);Mp=obj.ibrGFM_Mp(:);  Mq=obj.ibrGFM_Mq(:);
pset=obj.ibrGFM_pset(:);qset=obj.ibrGFM_qset(:);voQs=obj.ibrGFM_vo_q_set(:);
wVRR=obj.ibrGFM_w_VR(:);Rv=obj.ibrGFM_Rv(:);
w0=obj.w0; dT=obj.deltaT; jw=1/w0; dTPU=dT*w0;
b2I=obj.busToIbrGFM_I(BUS(:),(1:N).')+zeros(N,1);
b2V=obj.busToIbrGFM_V(BUS(:),(1:N).')+zeros(N,1);
z=zeros(N,1); I=zeros(N,obj.numIbrGFMEqns);

% ---- Park transform: T = 2/3*[cos;sin], JT = 2/3*[-sin;cos] ------------
c1=cos(th); c2=cos(th-2*pi/3); c3=cos(th+2*pi/3);
s1=sin(th); s2=sin(th-2*pi/3); s3=sin(th+2*pi/3);
T11=2/3*c1; T12=2/3*c2; T13=2/3*c3;
T21=2/3*s1; T22=2/3*s2; T23=2/3*s3;
J11=2/3*(-s1); J12=2/3*(-s2); J13=2/3*(-s3);
J21=2/3*c1;    J22=2/3*c2;    J23=2/3*c3;

% ---- Park I -------------------------------------------------------------
jti1=((J11.*ia+J12.*ib)+J13.*ic).*b2I;  jti2=((J21.*ia+J22.*ib)+J23.*ic).*b2I;
P14=T11.*b2I; P15=T12.*b2I; P16=T13.*b2I;
P24=T21.*b2I; P25=T22.*b2I; P26=T23.*b2I;
Ti1=((T11.*ia+T12.*ib)+T13.*ic).*b2I;   Ti2=((T21.*ia+T22.*ib)+T23.*ic).*b2I;
F1=-iod+Ti1; F2=-ioq+Ti2; F0=((ia+ib)+ic).*b2I;
I(:,1)=(((((-1).*iod+z.*ioq)+jti1.*th)+P14.*ia)+P15.*ib)+P16.*ic - F1;
I(:,2)=(((((z.*iod+(-1).*ioq)+jti2.*th)+P24.*ia)+P25.*ib)+P26.*ic) - F2;
I(:,3)=(b2I.*ia+(b2I.*ib+b2I.*ic)) - F0;
all_stampVal(PS(:,1))=-1; all_stampVal(PS(:,2))=0; all_stampVal(PS(:,3))=0; all_stampVal(PS(:,4))=-1;
all_stampVal(PS(:,5))=jti1; all_stampVal(PS(:,6))=jti2;
all_stampVal(PS(:,7))=P14; all_stampVal(PS(:,8))=P24;
all_stampVal(PS(:,9))=P15; all_stampVal(PS(:,10))=P25;
all_stampVal(PS(:,11))=P16;all_stampVal(PS(:,12))=P26;
all_stampVal(PS(:,13))=b2I;all_stampVal(PS(:,14))=b2I;all_stampVal(PS(:,15))=b2I;

% ---- Park V -------------------------------------------------------------
jtv1=((J11.*va+J12.*vb)+J13.*vc).*b2V;  jtv2=((J21.*va+J22.*vb)+J23.*vc).*b2V;
V11=T11.*b2V; V12=T12.*b2V; V13=T13.*b2V;
V21=T21.*b2V; V22=T22.*b2V; V23=T23.*b2V;
Tv1=((T11.*va+T12.*vb)+T13.*vc).*b2V;   Tv2=((T21.*va+T22.*vb)+T23.*vc).*b2V;
FV1=-vgd+Tv1; FV2=-vgq+Tv2;
I(:,4)=(((-1).*vgd+z.*vgq)+jtv1.*th) + ((V11.*va+V12.*vb)+V13.*vc) - FV1;
I(:,5)=((z.*vgd+(-1).*vgq)+jtv2.*th) + ((V21.*va+V22.*vb)+V23.*vc) - FV2;
all_stampVal(PS(:,16))=-1; all_stampVal(PS(:,17))=0;
all_stampVal(PS(:,18))=0;  all_stampVal(PS(:,19))=-1;
all_stampVal(PS(:,20))=jtv1; all_stampVal(PS(:,21))=jtv2;
all_stampVal(PS(:,169))=V11; all_stampVal(PS(:,170))=V21;
all_stampVal(PS(:,171))=V12; all_stampVal(PS(:,172))=V22;
all_stampVal(PS(:,173))=V13; all_stampVal(PS(:,174))=V23;

% ---- Passive coupling ---------------------------------------------------
kC=-1/dTPU+(1./Lc).*(-rc/2);  gC=(1./Lc)/2;
A11=kC; A12=(-jw.*wr)/2; A21=(jw.*wr)/2; A22=kC;
c31=(-jw.*ioq)/2; c32=(jw.*iod)/2;
sC1=wr.*iod+wrp.*iodp; sC2=wr.*ioq+wrp.*ioqp;
FC1=-(iod-iodp)/dTPU+(1./Lc).*(-rc.*(iod+iodp)/2+(vod+vodp)/2-(vgd+vgdp)/2)+(-jw.*sC2)/2;
FC2=-(ioq-ioqp)/dTPU+(1./Lc).*(-rc.*(ioq+ioqp)/2+(voq+voqp)/2-(vgq+vgqp)/2)+( jw.*sC1)/2;
I(:,6)=((((((A11.*iod+A12.*ioq)+c31.*wr)+gC.*vod)+z.*voq)+(-gC).*vgd)+z.*vgq)-FC1;
I(:,7)=((((((A21.*iod+A22.*ioq)+c32.*wr)+z.*vod)+gC.*voq)+z.*vgd)+(-gC).*vgq)-FC2;
all_stampVal(PS(:,22))=A11; all_stampVal(PS(:,23))=A21;
all_stampVal(PS(:,24))=A12; all_stampVal(PS(:,25))=A22;
all_stampVal(PS(:,26))=c31; all_stampVal(PS(:,27))=c32;
all_stampVal(PS(:,28))=gC;  all_stampVal(PS(:,29))=0;
all_stampVal(PS(:,30))=0;   all_stampVal(PS(:,31))=gC;
all_stampVal(PS(:,32))=-gC; all_stampVal(PS(:,33))=0;
all_stampVal(PS(:,34))=0;   all_stampVal(PS(:,35))=-gC;

% ---- Passive capacitor --------------------------------------------------
kP=-1/dTPU; gI=(1./Cf)/2+Rcp/dTPU; gO=(1./Cf)*(-1/2)+Rcp*(-1/dTPU);
B11=kP; B12=(-jw.*wr)/2; B21=(jw.*wr)/2; B22=kP;
d31=(-jw.*voq)/2; d32=(jw.*vod)/2;
sP1=wr.*vod+wrp.*vodp; sP2=wr.*voq+wrp.*voqp;
FP1=-(vod-vodp)/dTPU+(1./Cf).*((ifd+ifdp)/2-(iod+iodp)/2)+(-jw.*sP2)/2+Rcp.*((ifd-ifdp)/dTPU-(iod-iodp)/dTPU);
FP2=-(voq-voqp)/dTPU+(1./Cf).*((ifq+ifqp)/2-(ioq+ioqp)/2)+( jw.*sP1)/2+Rcp.*((ifq-ifqp)/dTPU-(ioq-ioqp)/dTPU);
I(:,8)=((((((B11.*vod+B12.*voq)+d31.*wr)+gI.*ifd)+z.*ifq)+gO.*iod)+z.*ioq)-FP1;
I(:,9)=((((((B21.*vod+B22.*voq)+d32.*wr)+z.*ifd)+gI.*ifq)+z.*iod)+gO.*ioq)-FP2;
all_stampVal(PS(:,36))=B11; all_stampVal(PS(:,37))=B21;
all_stampVal(PS(:,38))=B12; all_stampVal(PS(:,39))=B22;
all_stampVal(PS(:,40))=d31; all_stampVal(PS(:,41))=d32;
all_stampVal(PS(:,42))=gI;  all_stampVal(PS(:,43))=0;
all_stampVal(PS(:,44))=0;   all_stampVal(PS(:,45))=gI;
all_stampVal(PS(:,46))=gO;  all_stampVal(PS(:,47))=0;
all_stampVal(PS(:,48))=0;   all_stampVal(PS(:,49))=gO;

% ---- Passive filter -----------------------------------------------------
kF=-1/dTPU+(1./Lf).*(-rf/2); gF=(1./Lf)/2; gFn=(1./Lf)*(-1/2);
C11=kF; C12=(-jw.*wr)/2; C21=(jw.*wr)/2; C22=kF;
e31=(-jw.*ifq)/2; e32=(jw.*ifd)/2;
sF1=wr.*ifd+wrp.*ifdp; sF2=wr.*ifq+wrp.*ifqp;
FF1=-(ifd-ifdp)/dTPU+(1./Lf).*(-rf.*(ifd+ifdp)/2+(vid+vidp)/2-(vod+vodp)/2)+(-jw.*sF2)/2;
FF2=-(ifq-ifqp)/dTPU+(1./Lf).*(-rf.*(ifq+ifqp)/2+(viq+viqp)/2-(voq+voqp)/2)+( jw.*sF1)/2;
I(:,10)=((((((C11.*ifd+C12.*ifq)+e31.*wr)+gF.*vid)+z.*viq)+gFn.*vod)+z.*voq)-FF1;
I(:,11)=((((((C21.*ifd+C22.*ifq)+e32.*wr)+z.*vid)+gF.*viq)+z.*vod)+gFn.*voq)-FF2;
all_stampVal(PS(:,50))=C11; all_stampVal(PS(:,51))=C21;
all_stampVal(PS(:,52))=C12; all_stampVal(PS(:,53))=C22;
all_stampVal(PS(:,54))=e31; all_stampVal(PS(:,55))=e32;
all_stampVal(PS(:,56))=gF;  all_stampVal(PS(:,57))=0;
all_stampVal(PS(:,58))=0;   all_stampVal(PS(:,59))=gF;
all_stampVal(PS(:,60))=gFn; all_stampVal(PS(:,61))=0;
all_stampVal(PS(:,62))=0;   all_stampVal(PS(:,63))=gFn;

% ---- Current controller -------------------------------------------------
% -stamp_jw0*Lf = [0 jw*Lf; -jw*Lf 0]
mL=jw*Lf; kG=kci/2+kcp/dT; gGc=Gc/2;
h31=(mL.*wr)/2;  h32=(-mL.*wr)/2;      % (-stamp_jw0*Lf*wr/2) entries (1,2)&(2,1)
h41=(mL.*ifq)/2; h42=(-mL.*ifd)/2;     % (-stamp_jw0*Lf*if/2)
sK1=wr.*ifd+wrp.*ifdp; sK2=wr.*ifq+wrp.*ifqp;
FK1=-(vid+vidp)/2+kci.*(gd+gdp)/2+kcp.*(gd-gdp)/dT+(mL.*sK2)/2+Gc.*(vod+vodp)/2;
FK2=-(viq+viqp)/2+kci.*(gq+gqp)/2+kcp.*(gq-gqp)/dT+(-mL.*sK1)/2+Gc.*(voq+voqp)/2;
I(:,12)=((((((((-1/2).*vid+z.*viq)+kG.*gd)+z.*gq)+z.*ifd)+h31.*ifq)+h41.*wr)+gGc.*vod)+z.*voq-FK1;
I(:,13)=((((((((z.*vid+(-1/2).*viq)+z.*gd)+kG.*gq)+h32.*ifd)+z.*ifq)+h42.*wr)+z.*vod)+gGc.*voq)-FK2;
FKi1=-(gd-gdp)/dT+(frd+frdp)/2-(ifd+ifdp)/2;
FKi2=-(gq-gqp)/dT+(frq+frqp)/2-(ifq+ifqp)/2;
mdt=-1/dT;
I(:,14)=(((((mdt.*gd+z.*gq)+(-1/2).*ifd)+z.*ifq)+(1/2).*frd)+z.*frq)-FKi1;
I(:,15)=(((((z.*gd+mdt.*gq)+z.*ifd)+(-1/2).*ifq)+z.*frd)+(1/2).*frq)-FKi2;
all_stampVal(PS(:,64))=-1/2; all_stampVal(PS(:,65))=0;
all_stampVal(PS(:,66))=0;    all_stampVal(PS(:,67))=-1/2;
all_stampVal(PS(:,68))=kG;   all_stampVal(PS(:,69))=0;
all_stampVal(PS(:,70))=0;    all_stampVal(PS(:,71))=kG;
all_stampVal(PS(:,72))=0;    all_stampVal(PS(:,73))=h32;
all_stampVal(PS(:,74))=h31;  all_stampVal(PS(:,75))=0;
all_stampVal(PS(:,76))=h41;  all_stampVal(PS(:,77))=h42;
all_stampVal(PS(:,78))=gGc;  all_stampVal(PS(:,79))=0;
all_stampVal(PS(:,80))=0;    all_stampVal(PS(:,81))=gGc;
all_stampVal(PS(:,82))=mdt;  all_stampVal(PS(:,83))=0;
all_stampVal(PS(:,84))=0;    all_stampVal(PS(:,85))=mdt;
all_stampVal(PS(:,86))=-1/2; all_stampVal(PS(:,87))=0;
all_stampVal(PS(:,88))=0;    all_stampVal(PS(:,89))=-1/2;
all_stampVal(PS(:,90))=1/2;  all_stampVal(PS(:,91))=0;
all_stampVal(PS(:,92))=0;    all_stampVal(PS(:,93))=1/2;

% ---- Voltage controller -------------------------------------------------
mC=jw*Cf; kV=kvi/2+kvp/dT; gGv=Gv/2;
p31=(mC.*wr)/2;  p32=(-mC.*wr)/2;
p41=(mC.*voq)/2; p42=(-mC.*vod)/2;
sV1=wr.*vod+wrp.*vodp; sV2=wr.*voq+wrp.*voqp;
FV1c=-(frd+frdp)/2+kvi.*(xd+xdp)/2+kvp.*(xd-xdp)/dT+(mC.*sV2)/2+Gv.*(iod+iodp)/2;
FV2c=-(frq+frqp)/2+kvi.*(xq+xqp)/2+kvp.*(xq-xqp)/dT+(-mC.*sV1)/2+Gv.*(ioq+ioqp)/2;
I(:,16)=((((((((-1/2).*frd+z.*frq)+kV.*xd)+z.*xq)+z.*vod)+p31.*voq)+p41.*wr)+gGv.*iod)+z.*ioq-FV1c;
I(:,17)=((((((((z.*frd+(-1/2).*frq)+z.*xd)+kV.*xq)+p32.*vod)+z.*voq)+p42.*wr)+z.*iod)+gGv.*ioq)-FV2c;
FVi1=-(xd-xdp)/dT+(vrd+vrdp)/2-(vod+vodp)/2;
FVi2=-(xq-xqp)/dT+(vrq+vrqp)/2-(voq+voqp)/2;
I(:,18)=(((((mdt.*xd+z.*xq)+(-1/2).*vod)+z.*voq)+(1/2).*vrd)+z.*vrq)-FVi1;
I(:,19)=(((((z.*xd+mdt.*xq)+z.*vod)+(-1/2).*voq)+z.*vrd)+(1/2).*vrq)-FVi2;
all_stampVal(PS(:,94))=-1/2; all_stampVal(PS(:,95))=0;
all_stampVal(PS(:,96))=0;    all_stampVal(PS(:,97))=-1/2;
all_stampVal(PS(:,98))=kV;   all_stampVal(PS(:,99))=0;
all_stampVal(PS(:,100))=0;   all_stampVal(PS(:,101))=kV;
all_stampVal(PS(:,102))=0;   all_stampVal(PS(:,103))=p32;
all_stampVal(PS(:,104))=p31; all_stampVal(PS(:,105))=0;
all_stampVal(PS(:,106))=p41; all_stampVal(PS(:,107))=p42;
all_stampVal(PS(:,108))=gGv; all_stampVal(PS(:,109))=0;
all_stampVal(PS(:,110))=0;   all_stampVal(PS(:,111))=gGv;
all_stampVal(PS(:,112))=mdt; all_stampVal(PS(:,113))=0;
all_stampVal(PS(:,114))=0;   all_stampVal(PS(:,115))=mdt;
all_stampVal(PS(:,116))=-1/2;all_stampVal(PS(:,117))=0;
all_stampVal(PS(:,118))=0;   all_stampVal(PS(:,119))=-1/2;
all_stampVal(PS(:,120))=1/2; all_stampVal(PS(:,121))=0;
all_stampVal(PS(:,122))=0;   all_stampVal(PS(:,123))=1/2;

% ---- Droop --------------------------------------------------------------
wm=wmR/(2*pi);
DW1=-1; DW2=-Mp; DV1=-1; DV2=-Mq;
Dp1=-1/dT+wm*(-1/2);
Dp2=(wm.*iod)/2; Dp3=(wm.*ioq)/2; Dp4=(wm.*vgd)/2; Dp5=(wm.*vgq)/2;
Dq2=(-(wm.*ioq))/2; Dq3=(wm.*iod)/2; Dq4=(wm.*vgq)/2; Dq5=(-(wm.*vgd))/2;
Dt1=-1/dT; Dt2=1/2;
qp_=vgdp.*(-ioqp)+vgqp.*iodp;
pp_=vgdp.*iodp+vgqp.*ioqp;
vgio =vgd.*iod+vgq.*ioq;
vgjio=vgd.*(-ioq)+vgq.*iod;
FW=-wr+w0+Mp.*(pset-pav);
FVd=-vqd+voQs+Mq.*(qset-qav);
FP=-(pav-pavp)/dT+wm.*((vgio+pp_)/2-(pav+pavp)/2);
FQ=-(qav-qavp)/dT+wm.*((vgjio+qp_)/2-(qav+qavp)/2);
FT=-(th-thp)/dT+(wr+wrp)/2;
I(:,20)=DW1.*wr+DW2.*pav-FW;
I(:,21)=DV1.*vqd+DV2.*qav-FVd;
I(:,22)=(Dp1.*pav+(Dp2.*vgd+Dp3.*vgq))+(Dp4.*iod+Dp5.*ioq)-FP;
I(:,23)=(Dp1.*qav+(Dq2.*vgd+Dq3.*vgq))+(Dq4.*iod+Dq5.*ioq)-FQ;
I(:,24)=Dt1.*th+Dt2.*wr-FT;
all_stampVal(PS(:,124))=DW1; all_stampVal(PS(:,125))=DW2;
all_stampVal(PS(:,126))=DV1; all_stampVal(PS(:,127))=DV2;
all_stampVal(PS(:,128))=Dp1; all_stampVal(PS(:,129))=Dp2;
all_stampVal(PS(:,130))=Dp3; all_stampVal(PS(:,131))=Dp4; all_stampVal(PS(:,132))=Dp5;
all_stampVal(PS(:,133))=Dp1; all_stampVal(PS(:,134))=Dq2;
all_stampVal(PS(:,135))=Dq3; all_stampVal(PS(:,136))=Dq4; all_stampVal(PS(:,137))=Dq5;
all_stampVal(PS(:,138))=Dt1; all_stampVal(PS(:,139))=Dt2;

% ---- Virtual resistor ---------------------------------------------------
wVR=wVRR/(2*pi);
Ra=-1/dT-wVR/2; Rb=wVR.*Rv/2;
FRa1=-(ad-adp)/dT+wVR.*(Rv.*(iod+iodp)/2-(ad+adp)/2);
FRa2=-(aq-aqp)/dT+wVR.*(Rv.*(ioq+ioqp)/2-(aq+aqp)/2);
FR1=-rd+ad-Rv.*iod;  FR2=-rq+aq-Rv.*ioq;
FRr1=-vrd+z+rd;      FRr2=-vrq+vqd+rq;
I(:,25)=(((Ra.*ad+z.*aq)+Rb.*iod)+z.*ioq)-FRa1;
I(:,26)=(((z.*ad+Ra.*aq)+z.*iod)+Rb.*ioq)-FRa2;
I(:,27)=(((((-1).*rd+z.*rq)+1.*ad)+z.*aq)+(-Rv).*iod)+z.*ioq-FR1;
I(:,28)=(((((z.*rd+(-1).*rq)+z.*ad)+1.*aq)+z.*iod)+(-Rv).*ioq)-FR2;
I(:,29)=(((-1).*vrd+z.*vrq)+1.*rd)+z.*rq-FRr1;
I(:,30)=((((z.*vrd+(-1).*vrq)+z.*rd)+1.*rq)-FRr2)+1.*vqd;
all_stampVal(PS(:,140))=Ra; all_stampVal(PS(:,141))=0;
all_stampVal(PS(:,142))=0;  all_stampVal(PS(:,143))=Ra;
all_stampVal(PS(:,144))=Rb; all_stampVal(PS(:,145))=0;
all_stampVal(PS(:,146))=0;  all_stampVal(PS(:,147))=Rb;
all_stampVal(PS(:,148))=-1; all_stampVal(PS(:,149))=0;
all_stampVal(PS(:,150))=0;  all_stampVal(PS(:,151))=-1;
all_stampVal(PS(:,152))=1;  all_stampVal(PS(:,153))=0;
all_stampVal(PS(:,154))=0;  all_stampVal(PS(:,155))=1;
all_stampVal(PS(:,156))=-Rv;all_stampVal(PS(:,157))=0;
all_stampVal(PS(:,158))=0;  all_stampVal(PS(:,159))=-Rv;
all_stampVal(PS(:,160))=-1; all_stampVal(PS(:,161))=0;
all_stampVal(PS(:,162))=0;  all_stampVal(PS(:,163))=-1;
all_stampVal(PS(:,164))=1;  all_stampVal(PS(:,165))=0;
all_stampVal(PS(:,166))=0;  all_stampVal(PS(:,167))=1;
all_stampVal(PS(:,168))=1;


I_v(ROW.') = I.';

end
