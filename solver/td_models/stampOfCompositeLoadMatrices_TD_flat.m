function [all_stampVal,I_v] = stampOfCompositeLoadMatrices_TD_flat( ...
    obj,all_stampVal,I_v,x_v,x_p)
% COMPOSITE LOAD stamp - flat Nx1 column form, vectorised over all loads.
%
% Every quantity is a plain Nx1 column and every matrix product is
% scalar-expanded over the FIXED sparsity of R, L and L_Te.
%
% Terms are accumulated in a FIXED order that matches MATLAB's BLAS
% reduction order, so results are reproducible to the last bit:
%   matrix x column : left-to-right over columns
%   row x column    : pairwise tree (odd K: t1+(t2+t3) then +(t4+t5), ...)
% Do not reorder them; the dropped terms are structural zeros.
%
% Constant sparsity:
%   R    = diag(Rs,Rs,Rr,Rr)
%   L    = [Ls 0 Lm 0; 0 Ls 0 Lm; Lm 0 Lr 0; 0 Lm 0 Lr]
%   L_Te has nonzeros at (1,4),(2,3),(3,4),(4,3) only

N = obj.numLoadComp;

% ---- Index maps, whole population in one shot ---------------------------
RW  = obj.getRowIdx_fromLoadCompIndices((1:N).')-1;                       % N x 1
PS  = obj.getStampIdx_fromLoadCompIndices((1:N).')-1 + ...
    (1:obj.loadComp_numStampVals);
ROW = RW + (1:obj.numLoadCompEqns);
[xVs,xIabc,xWsl,xThS,xIs,xIr,xWr,xPQ,xVmP,xVm,xIz] = getLoadCompIndicesTD;
VA  = obj.getRowIdxAbc_fromBusNums(obj.loadComp_bus.');                   % 3 x N
[~,BUS] = ismembertol(obj.loadComp_bus,obj.bus);

% ---- States as flat Nx1 columns -----------------------------------------
vs1=x_v(RW+xVs(1));   vs2=x_v(RW+xVs(2));
ia =x_v(RW+xIabc(1)); ib =x_v(RW+xIabc(2)); ic =x_v(RW+xIabc(3));
wsl=x_v(RW+xWsl);     thS=x_v(RW+xThS);
is1=x_v(RW+xIs(1));   is2=x_v(RW+xIs(2));
ir1=x_v(RW+xIr(1));   ir2=x_v(RW+xIr(2));
wr =x_v(RW+xWr);
PQ1=x_v(RW+xPQ(1));   PQ2=x_v(RW+xPQ(2));
VmP=x_v(RW+xVmP);     Vm =x_v(RW+xVm);
iz1=x_v(RW+xIz(1));   iz2=x_v(RW+xIz(2));
va =x_v(VA(1,:).');   vb =x_v(VA(2,:).');   vc =x_v(VA(3,:).');
wS  = x_v(obj.getRowIdx_wSys);
%
vsp1=x_p(RW+xVs(1));  vsp2=x_p(RW+xVs(2));
wslp=x_p(RW+xWsl);    thSp=x_p(RW+xThS);
isp1=x_p(RW+xIs(1));  isp2=x_p(RW+xIs(2));
irp1=x_p(RW+xIr(1));  irp2=x_p(RW+xIr(2));
wrp =x_p(RW+xWr);
VmPp=x_p(RW+xVmP);    Vmp =x_p(RW+xVm);
wSp = x_p(obj.getRowIdx_wSys);

% ---- Constants flattened to columns -------------------------------------
w0 = obj.w0;   dT = obj.deltaT;
iDt = 1/dT;    w0i = 1/w0;
R1 = reshape(obj.loadCompIM_R(1,1,:),[],1);
R3 = reshape(obj.loadCompIM_R(3,3,:),[],1);
Ls = reshape(obj.loadCompIM_L(1,1,:),[],1)*w0i;
Lm = reshape(obj.loadCompIM_L(1,3,:),[],1)*w0i;
Lr = reshape(obj.loadCompIM_L(3,3,:),[],1)*w0i;
Te14 = reshape(obj.loadCompIM_LTe(1,4,:),[],1);
Te23 = reshape(obj.loadCompIM_LTe(2,3,:),[],1);
Te34 = reshape(obj.loadCompIM_LTe(3,4,:),[],1);
Te43 = reshape(obj.loadCompIM_LTe(4,3,:),[],1);
H  = obj.loadCompIM_H(:);   T0 = obj.loadCompIM_Tm0(:);   m = obj.loadCompIM_m(:);
PQz1=obj.loadCompZip_PQz(1,:).';  PQz2=obj.loadCompZip_PQz(2,:).';
PQi1=obj.loadCompZip_PQi(1,:).';  PQi2=obj.loadCompZip_PQi(2,:).';
PQp1=obj.loadCompZip_PQp(1,:).';  PQp2=obj.loadCompZip_PQp(2,:).';
iTf = 1./obj.loadCompZip_Tau(:);
b2I = obj.busToLoadComp_I(BUS(:),(1:N).') + zeros(N,1);
b2V = obj.busToLoadComp_V(BUS(:),(1:N).') + zeros(N,1);

I = zeros(N,obj.numLoadCompEqns);

% ------------------------- Induction motor -------------------------------
% Park transform:  T = 2/3*[cos; -sin] , JT = 2/3*[-sin; -cos]
c1=cos(thS); c2=cos(thS-2*pi/3); c3=cos(thS+2*pi/3);
s1=sin(thS); s2=sin(thS-2*pi/3); s3=sin(thS+2*pi/3);
T11= 2/3*c1;  T12= 2/3*c2;  T13= 2/3*c3;
T21=2/3*(-s1);T22=2/3*(-s2);T23=2/3*(-s3);
J11=2/3*(-s1);J12=2/3*(-s2);J13=2/3*(-s3);
J21=2/3*(-c1);J22=2/3*(-c2);J23=2/3*(-c3);

% --- Park V ---   tmpV = (JT*vabc)*b2V   (gemv: left-to-right over columns)
tmpV1 = ((J11.*va + J12.*vb) + J13.*vc).*b2V;
tmpV2 = ((J21.*va + J22.*vb) + J23.*vc).*b2V;
PV11=T11.*b2V; PV12=T12.*b2V; PV13=T13.*b2V;
PV21=T21.*b2V; PV22=T22.*b2V; PV23=T23.*b2V;
PVv1 = (PV11.*va + PV12.*vb) + PV13.*vc;
PVv2 = (PV21.*va + PV22.*vb) + PV23.*vc;
FParkV1 = -vs1 + PVv1;      FParkV2 = -vs2 + PVv2;
I(:,1) = -vs1 + tmpV1.*thS + PVv1 - FParkV1;
I(:,2) = -vs2 + tmpV2.*thS + PVv2 - FParkV2;
all_stampVal(PS(:,1)) = -1;   all_stampVal(PS(:,2)) = 0;
all_stampVal(PS(:,3)) =  0;   all_stampVal(PS(:,4)) = -1;
all_stampVal(PS(:,5)) = tmpV1; all_stampVal(PS(:,6)) = tmpV2;

% --- Park I ---
tmpI1 = ((J11.*ia + J12.*ib) + J13.*ic).*b2I;
tmpI2 = ((J21.*ia + J22.*ib) + J23.*ic).*b2I;
PI11=T11.*b2I; PI12=T12.*b2I; PI13=T13.*b2I;
PI21=T21.*b2I; PI22=T22.*b2I; PI23=T23.*b2I;
PIi1 = (PI11.*ia + PI12.*ib) + PI13.*ic;
PIi2 = (PI21.*ia + PI22.*ib) + PI23.*ic;
PI0  = b2I;                                   % ones(1,3)*b2I, each entry
FParkI1 = -(is1 + iz1) + ((T11.*ia + T12.*ib) + T13.*ic).*b2I;
FParkI2 = -(is2 + iz2) + ((T21.*ia + T22.*ib) + T23.*ic).*b2I;
FParkI0 = ((ia + ib) + ic).*b2I;
I(:,3) = -is1 - iz1 + tmpI1.*thS + PIi1 - FParkI1;
I(:,4) = -is2 - iz2 + tmpI2.*thS + PIi2 - FParkI2;
% 1x3 row-times-column dot (odd K): t1+(t2+t3)
I(:,5) = (PI0.*ia + (PI0.*ib + PI0.*ic)) - FParkI0;
all_stampVal(PS(:, 7)) = -1;   all_stampVal(PS(:, 8)) = 0;
all_stampVal(PS(:, 9)) =  0;   all_stampVal(PS(:,10)) = -1;
all_stampVal(PS(:,11)) = -1;   all_stampVal(PS(:,12)) = 0;
all_stampVal(PS(:,13)) =  0;   all_stampVal(PS(:,14)) = -1;
all_stampVal(PS(:,15)) = tmpI1; all_stampVal(PS(:,16)) = tmpI2;
all_stampVal(PS(:,17)) = PI11; all_stampVal(PS(:,18)) = PI21;
all_stampVal(PS(:,19)) = PI12; all_stampVal(PS(:,20)) = PI22;
all_stampVal(PS(:,21)) = PI13; all_stampVal(PS(:,22)) = PI23;
all_stampVal(PS(:,23)) = PI0;  all_stampVal(PS(:,24)) = PI0;
all_stampVal(PS(:,25)) = PI0;

% --- Slip frequency ---
Fslip = -wsl + wS - wr;
I(:,6) = -wsl - wr + wS - Fslip;
all_stampVal(PS(:,26)) = -1;   all_stampVal(PS(:,27)) = -1;

% --- System angle ---
Fsys = -(thS - thSp)*iDt + (wS + wSp)*0.5;
I(:,7) = (-iDt)*thS + 0.5*wS - Fsys;
all_stampVal(PS(:,28)) = -iDt;

% --- Stator / rotor ---   phisr = L_PU*isr  (gemv, zeros dropped exactly)
phi1 = Ls.*is1 + Lm.*ir1;   phi2 = Ls.*is2 + Lm.*ir2;
phi3 = Lm.*is1 + Lr.*ir1;   phi4 = Lm.*is2 + Lr.*ir2;
phip1= Ls.*isp1 + Lm.*irp1; phip2= Ls.*isp2 + Lm.*irp2;
phip3= Lm.*isp1 + Lr.*irp1; phip4= Lm.*isp2 + Lr.*irp2;
% dFsr_mid = (R + wMat*L_PU)*0.5 + L_PU*iDt   (wMat has one nonzero per row)
m11 = R1*0.5 + Ls*iDt;          m12 = (-wS.*Ls)*0.5;
m13 = Lm*iDt;                   m14 = (-wS.*Lm)*0.5;
m21 = (wS.*Ls)*0.5;             m22 = R1*0.5 + Ls*iDt;
m23 = (wS.*Lm)*0.5;             m24 = Lm*iDt;
m31 = Lm*iDt;                   m32 = (-wsl.*Lm)*0.5;
m33 = R3*0.5 + Lr*iDt;          m34 = (-wsl.*Lr)*0.5;
m41 = (wsl.*Lm)*0.5;            m42 = Lm*iDt;
m43 = (wsl.*Lr)*0.5;            m44 = R3*0.5 + Lr*iDt;
c37 = -phi4*0.5;                c47 =  phi3*0.5;
dFs_dwS1 = -phi2*0.5;           dFs_dwS2 =  phi1*0.5;
% Fsr = -(vsr+vsr_p)/2 + R*(isr+isr_p)/2 + L_PU*(isr-isr_p)/dT
%       + (wMat*phisr + wMat_p*(L_PU*isr_p))/2
dis1=is1-isp1; dis2=is2-isp2; dir1=ir1-irp1; dir2=ir2-irp2;
si1=is1+isp1;  si2=is2+isp2;  ri1=ir1+irp1;  ri2=ir2+irp2;
Fsr1 = -(vs1+vsp1)*0.5 + R1.*si1*0.5 + (Ls.*dis1 + Lm.*dir1)*iDt + ...
    ((-wS.*phi2) + (-wSp.*phip2))*0.5;
Fsr2 = -(vs2+vsp2)*0.5 + R1.*si2*0.5 + (Ls.*dis2 + Lm.*dir2)*iDt + ...
    (( wS.*phi1) + ( wSp.*phip1))*0.5;
Fsr3 =                   R3.*ri1*0.5 + (Lm.*dis1 + Lr.*dir1)*iDt + ...
    ((-wsl.*phi4) + (-wslp.*phip4))*0.5;
Fsr4 =                   R3.*ri2*0.5 + (Lm.*dis2 + Lr.*dir2)*iDt + ...
    (( wsl.*phi3) + ( wslp.*phip3))*0.5;
% I(8:11) = dFsr(:,1:2)*vs + mid*isr + dFsr(:,7)*wSlip - Fsr  (+ dFs_dwSys*wSys)
I(:,8)  = (-0.5)*vs1 + (((m11.*is1 + m12.*is2) + m13.*ir1) + m14.*ir2) ...
    + 0.*wsl - Fsr1 + dFs_dwS1*wS;
I(:,9)  = (-0.5)*vs2 + (((m21.*is1 + m22.*is2) + m23.*ir1) + m24.*ir2) ...
    + 0.*wsl - Fsr2 + dFs_dwS2*wS;
I(:,10) =              (((m31.*is1 + m32.*is2) + m33.*ir1) + m34.*ir2) ...
    + c37.*wsl - Fsr3;
I(:,11) =              (((m41.*is1 + m42.*is2) + m43.*ir1) + m44.*ir2) ...
    + c47.*wsl - Fsr4;
z = zeros(N,1);
all_stampVal(PS(:,29))=-0.5; all_stampVal(PS(:,30))=0;    all_stampVal(PS(:,31))=0; all_stampVal(PS(:,32))=0;
all_stampVal(PS(:,33))=0;    all_stampVal(PS(:,34))=-0.5; all_stampVal(PS(:,35))=0; all_stampVal(PS(:,36))=0;
all_stampVal(PS(:,37))=m11;  all_stampVal(PS(:,38))=m21;  all_stampVal(PS(:,39))=m31; all_stampVal(PS(:,40))=m41;
all_stampVal(PS(:,41))=m12;  all_stampVal(PS(:,42))=m22;  all_stampVal(PS(:,43))=m32; all_stampVal(PS(:,44))=m42;
all_stampVal(PS(:,45))=m13;  all_stampVal(PS(:,46))=m23;  all_stampVal(PS(:,47))=m33; all_stampVal(PS(:,48))=m43;
all_stampVal(PS(:,49))=m14;  all_stampVal(PS(:,50))=m24;  all_stampVal(PS(:,51))=m34; all_stampVal(PS(:,52))=m44;
all_stampVal(PS(:,53))=z;    all_stampVal(PS(:,54))=z;    all_stampVal(PS(:,55))=c37; all_stampVal(PS(:,56))=c47;

% --- Swing ---
dwr  = -(2*H*w0i)*iDt - T0.*(m.*wr.^(m-1))./(wS.^m)*0.5;
% dFswing_disr = isr.'*(L_Te+L_Te.')*0.5, isr = [is1;is2;ir1;ir2].
% L_Te nonzeros: (1,4)=Te14 (2,3)=Te23 (3,4)=Te34 (4,3)=Te43, so
% M = L_Te+L_Te.' has (1,4)=Te14 (2,3)=Te23 (3,4)=(4,3)=Te34+Te43
%                     (4,1)=Te14 (3,2)=Te23
% columns of M:  c1<-row4 ; c2<-row3 ; c3<-rows 2,4 ; c4<-rows 1,3
% row-times-matrix uses the pairwise dot tree, even K=4: (t1+t2)+(t3+t4)
Tsum = Te34 + Te43;
d1 = (ir2.*Te14)*0.5;
d2 = (ir1.*Te23)*0.5;
d3 = (is2.*Te23 + ir2.*Tsum)*0.5;
d4 = (is1.*Te14 + ir1.*Tsum)*0.5;
dwS  = T0.*(m.*wr.^m)./(wS.^(m+1))*0.5;
% isr.'*L_Te*isr : columns of L_Te are c3<-rows 2,4 and c4<-rows 1,3
rowL3 = is2.*Te23 + ir2.*Te43;   rowL4 = is1.*Te14 + ir1.*Te34;
iLi   = (z.*is1 + z.*is2) + (rowL3.*ir1 + rowL4.*ir2);
rowP3 = isp2.*Te23 + irp2.*Te43; rowP4 = isp1.*Te14 + irp1.*Te34;
ipLip = (z.*isp1 + z.*isp2) + (rowP3.*irp1 + rowP4.*irp2);
Fsw = -(2*H*w0i).*(wr-wrp)*iDt + (iLi + ipLip)*0.5 - ...
    T0.*((wr./wS).^m + (wrp./wSp).^m)*0.5;
I(:,12) = dwr.*wr + ((d1.*is1 + d2.*is2) + (d3.*ir1 + d4.*ir2)) + dwS*wS - Fsw;
all_stampVal(PS(:,57))=dwr;
all_stampVal(PS(:,58))=d1; all_stampVal(PS(:,59))=d2;
all_stampVal(PS(:,60))=d3; all_stampVal(PS(:,61))=d4;

% ---------------------------- ZIP model ----------------------------------
dVm1 = 2*PQz1.*Vm + PQi1;   dVm2 = 2*PQz2.*Vm + PQi2;
FZ1 = -PQ1 + PQz1.*Vm.^2 + PQi1.*Vm + PQp1;
FZ2 = -PQ2 + PQz2.*Vm.^2 + PQi2.*Vm + PQp2;
I(:,13) = -PQ1 + dVm1.*Vm - FZ1;
I(:,14) = -PQ2 + dVm2.*Vm - FZ2;
all_stampVal(PS(:,62))=-1; all_stampVal(PS(:,63))=0;
all_stampVal(PS(:,64))= 0; all_stampVal(PS(:,65))=-1;
all_stampVal(PS(:,66))=dVm1; all_stampVal(PS(:,67))=dVm2;

% --- Vm pre-filter ---
vsn = sqrt(vs1.^2 + vs2.^2);
dP1 = vs1./vsn;   dP2 = vs2./vsn;
FVmP = -VmP + vsn;
I(:,15) = -VmP + (dP1.*vs1 + dP2.*vs2) - FVmP;
all_stampVal(PS(:,68))=-1; all_stampVal(PS(:,69))=dP1; all_stampVal(PS(:,70))=dP2;

% --- Vm filter ---
dFVm  = -iDt - iTf*0.5;   dFVmP = iTf*0.5;
FVm = -(Vm - Vmp)*iDt + iTf.*((VmP + VmPp)*0.5 - (Vm + Vmp)*0.5);
I(:,16) = dFVm.*Vm + dFVmP.*VmP - FVm;
all_stampVal(PS(:,71))=dFVm; all_stampVal(PS(:,72))=dFVmP;

% --- Izip ---
ize1 = vs1.*PQ1 + vs2.*PQ2;
ize2 = vs2.*PQ1 - vs1.*PQ2;
V2i = 1./(Vm.^2);   V3i = -2./(Vm.^3);
j73= vs1.*V2i; j74= vs2.*V2i; j75= vs2.*V2i; j76=-vs1.*V2i;
j77= PQ1.*V2i; j78=-PQ2.*V2i; j79= PQ2.*V2i; j80= PQ1.*V2i;
j81= ize1.*V3i; j82= ize2.*V3i;
FI1 = -iz1 + ize1.*V2i;   FI2 = -iz2 + ize2.*V2i;
all_stampVal(PS(:,73))=j73; all_stampVal(PS(:,74))=j74; all_stampVal(PS(:,75))=j75;
all_stampVal(PS(:,76))=j76; all_stampVal(PS(:,77))=j77; all_stampVal(PS(:,78))=j78;
all_stampVal(PS(:,79))=j79; all_stampVal(PS(:,80))=j80; all_stampVal(PS(:,81))=j81;
all_stampVal(PS(:,82))=j82;
all_stampVal(PS(:,83))=-1; all_stampVal(PS(:,84))=0;
all_stampVal(PS(:,85))= 0; all_stampVal(PS(:,86))=-1;
I(:,17) = j73.*PQ1 + j75.*PQ2 + j77.*vs1 + j79.*vs2 + j81.*Vm - iz1 - FI1;
I(:,18) = j74.*PQ1 + j76.*PQ2 + j78.*vs1 + j80.*vs2 + j82.*Vm - iz2 - FI2;

% ---- "Outside" index ----------------------------------------------------
all_stampVal(PS(:,87))=PV11; all_stampVal(PS(:,88))=PV21;
all_stampVal(PS(:,89))=PV12; all_stampVal(PS(:,90))=PV22;
all_stampVal(PS(:,91))=PV13; all_stampVal(PS(:,92))=PV23;
all_stampVal(PS(:,93))=1;    all_stampVal(PS(:,94))=0.5;
all_stampVal(PS(:,95))=dFs_dwS1; all_stampVal(PS(:,96))=dFs_dwS2;
all_stampVal(PS(:,97))=dwS;

I_v(ROW.') = I.';
end
