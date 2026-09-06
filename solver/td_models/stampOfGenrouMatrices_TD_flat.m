function [all_stampVal,I_v] = stampOfGenrouMatrices_TD_flat( ...
    obj,all_stampVal,I_v,x_v,x_p)
% GENROU stamp - flat Nx1 column form, vectorised over all machines.
% Terms are accumulated in a FIXED order that matches MATLAB's BLAS
% reduction order, so results are reproducible to the last bit:
%   matrix x column : left-to-right over columns
%   row x column    : left-to-right
% Structurally-zero terms are dropped, which is exact: adding 0*x to a
% left-to-right accumulation changes nothing.
%
% Constant sparsity:
%   R   = diag(R1..R7)
%   L   : rows 1,4,5 -> cols {1,4,5};  rows 2,6,7 -> cols {2,6,7};  row 3 -> {3}
%   Lsp : row 1 -> cols {2,6,7};  row 2 -> cols {1,4,5};  rows 3..7 all zero
%
% Local stamp numbering (column-major, exactly what the scalar produced):
%    1: 56  stator [-R+2L/dT_PU+Lsp*wr/ws , Lsp*i/ws]   (7x8)
%   57: 72  -eye(4)
%   73: 84  Park V [invT , JinvT*e(1:3)]                (3x4)
%   85: 93  swing   95: 94 speed
%   96:116  Park I [invT , JinvT*i(1:3) , -eye(3)*b2I]  (3x7)
%  117:125  Park V outside  -eye(3)*b2V                 (3x3)

N   = obj.numGenrou;
RW  = obj.getRowIdx_fromGenrouIndices((1:N).')-1;      % N x 1, device row base
SB  = obj.getStampIdx_fromGenrouIndices((1:N).')-1;    % N x 1, stamp base
VA  = obj.getRowIdxAbc_fromBusNums(obj.genrou_bus.');  % 3 x N
[~,BUS] = ismembertol(obj.genrou_bus,obj.bus);

% ---- States as flat Nx1 columns -----------------------------------------
i1=x_v(RW+1); i2=x_v(RW+2); i3=x_v(RW+3); i4=x_v(RW+4);
i5=x_v(RW+5); i6=x_v(RW+6); i7=x_v(RW+7);
e1=x_v(RW+8); e2=x_v(RW+9); e3=x_v(RW+10); e4=x_v(RW+11);
wr=x_v(RW+12); th=x_v(RW+13);
ia=x_v(RW+14); ib=x_v(RW+15); ic=x_v(RW+16); Tm=x_v(RW+17);
va=x_v(VA(1,:).'); vb=x_v(VA(2,:).'); vc=x_v(VA(3,:).');
%
p1=x_p(RW+1); p2=x_p(RW+2); p3=x_p(RW+3); p4=x_p(RW+4);
p5=x_p(RW+5); p6=x_p(RW+6); p7=x_p(RW+7);
q1=x_p(RW+8); q2=x_p(RW+9); q3=x_p(RW+10); q4=x_p(RW+11);
wrp=x_p(RW+12); thp=x_p(RW+13); Tmp=x_p(RW+17);

% ---- Constants flattened to columns -------------------------------------
Rc = reshape(obj.genrou_RMatrix,  49,N).';
Lc = reshape(obj.genrou_LMatrix,  49,N).';
Sc = reshape(obj.genrou_LspMatrix,49,N).';
% column-major: entry (r,c) is column (c-1)*7+r
R1=Rc(:,1); R2=Rc(:,9); R3=Rc(:,17); R4=Rc(:,25);
R5=Rc(:,33); R6=Rc(:,41); R7=Rc(:,49);
L11=Lc(:, 1); L14=Lc(:,22); L15=Lc(:,29);
L22=Lc(:, 9); L26=Lc(:,37); L27=Lc(:,44);
L33=Lc(:,17);
L41=Lc(:, 4); L44=Lc(:,25); L45=Lc(:,32);
L51=Lc(:, 5); L54=Lc(:,26); L55=Lc(:,33);
L62=Lc(:,13); L66=Lc(:,41); L67=Lc(:,48);
L72=Lc(:,14); L76=Lc(:,42); L77=Lc(:,49);
S12=Sc(:, 8); S16=Sc(:,36); S17=Sc(:,43);
S21=Sc(:, 2); S24=Sc(:,23); S25=Sc(:,30);
%
H  = obj.genrou_H(:);   KD = obj.genrou_KD(:);
ws = obj.w0;         dT = obj.deltaT;
b2V = obj.busToGenrou_V(BUS(:),(1:N).') + zeros(N,1);
b2I = obj.busToGenrou_I(BUS(:),(1:N).') + zeros(N,1);
% deltaT_PU = deltaT*[ws;ws;ws;1;1;1;1]
d3 = dT*ws;   d4 = dT;

% ---- Stator: E = -R + 2L/dT_PU ,  D = -R - 2L/dT_PU ---------------------
E11=-R1+2*L11/d3; E14=    2*L14/d3; E15=    2*L15/d3;
E22=-R2+2*L22/d3; E26=    2*L26/d3; E27=    2*L27/d3;
E33=-R3+2*L33/d3;
E41=    2*L41/d4; E44=-R4+2*L44/d4; E45=    2*L45/d4;
E51=    2*L51/d4; E54=    2*L54/d4; E55=-R5+2*L55/d4;
E62=    2*L62/d4; E66=-R6+2*L66/d4; E67=    2*L67/d4;
E72=    2*L72/d4; E76=    2*L76/d4; E77=-R7+2*L77/d4;
%
D11=-R1-2*L11/d3; D14=   -2*L14/d3; D15=   -2*L15/d3;
D22=-R2-2*L22/d3; D26=   -2*L26/d3; D27=   -2*L27/d3;
D33=-R3-2*L33/d3;
D41=   -2*L41/d4; D44=-R4-2*L44/d4; D45=   -2*L45/d4;
D51=   -2*L51/d4; D54=   -2*L54/d4; D55=-R5-2*L55/d4;
D62=   -2*L62/d4; D66=-R6-2*L66/d4; D67=   -2*L67/d4;
D72=   -2*L72/d4; D76=   -2*L76/d4; D77=-R7-2*L77/d4;

% Lsp*wr/ws  is  (Lsp(r,c)*wr)/ws
A12=(S12.*wr)/ws; A16=(S16.*wr)/ws; A17=(S17.*wr)/ws;
A21=(S21.*wr)/ws; A24=(S24.*wr)/ws; A25=(S25.*wr)/ws;
% Lsp*i and Lsp*i_prev  (gemv, left-to-right over columns)
u1 = (S12.*i2 + S16.*i6) + S17.*i7;
u2 = (S21.*i1 + S24.*i4) + S25.*i5;
v1 = (S12.*p2 + S16.*p6) + S17.*p7;
v2 = (S21.*p1 + S24.*p4) + S25.*p5;
c1 = u1/ws;   c2 = u2/ws;            % column 8 of the Jacobian

% ---- Stator stamps (7x8, column-major) ----------------------------------
all_stampVal(SB+ 1)=E11;  all_stampVal(SB+ 8)=A12;  all_stampVal(SB+22)=E14;
all_stampVal(SB+29)=E15;  all_stampVal(SB+36)=A16;  all_stampVal(SB+43)=A17;
all_stampVal(SB+ 2)=A21;  all_stampVal(SB+ 9)=E22;  all_stampVal(SB+23)=A24;
all_stampVal(SB+30)=A25;  all_stampVal(SB+37)=E26;  all_stampVal(SB+44)=E27;
all_stampVal(SB+17)=E33;
all_stampVal(SB+ 4)=E41;  all_stampVal(SB+25)=E44;  all_stampVal(SB+32)=E45;
all_stampVal(SB+ 5)=E51;  all_stampVal(SB+26)=E54;  all_stampVal(SB+33)=E55;
all_stampVal(SB+13)=E62;  all_stampVal(SB+41)=E66;  all_stampVal(SB+48)=E67;
all_stampVal(SB+14)=E72;  all_stampVal(SB+42)=E76;  all_stampVal(SB+49)=E77;
all_stampVal(SB+50)=c1;   all_stampVal(SB+51)=c2;
% -eye(4)
all_stampVal(SB+57)=-1; all_stampVal(SB+62)=-1;
all_stampVal(SB+67)=-1; all_stampVal(SB+72)=-1;

% ---- Stator RHS:  J*[i;wr] + [-e(1:4);0;0;0] + -F -----------------------
Jx1 = (((((E11.*i1 + A12.*i2) + E14.*i4) + E15.*i5) + A16.*i6) + A17.*i7) + c1.*wr;
Jx2 = (((((A21.*i1 + E22.*i2) + A24.*i4) + A25.*i5) + E26.*i6) + E27.*i7) + c2.*wr;
Jx3 = E33.*i3;
Jx4 = (E41.*i1 + E44.*i4) + E45.*i5;
Jx5 = (E51.*i1 + E54.*i4) + E55.*i5;
Jx6 = (E62.*i2 + E66.*i6) + E67.*i7;
Jx7 = (E72.*i2 + E76.*i6) + E77.*i7;
% F = ((((-e + E*i) + (Lsp*i*wr)/ws) + -e_prev) + D*i_prev) + (Lsp*i_p*wr_p)/ws
Ei1 = (E11.*i1 + E14.*i4) + E15.*i5;
Ei2 = (E22.*i2 + E26.*i6) + E27.*i7;
Ei3 =  E33.*i3;
Ei4 = (E41.*i1 + E44.*i4) + E45.*i5;
Ei5 = (E51.*i1 + E54.*i4) + E55.*i5;
Ei6 = (E62.*i2 + E66.*i6) + E67.*i7;
Ei7 = (E72.*i2 + E76.*i6) + E77.*i7;
Dp1 = (D11.*p1 + D14.*p4) + D15.*p5;
Dp2 = (D22.*p2 + D26.*p6) + D27.*p7;
Dp3 =  D33.*p3;
Dp4 = (D41.*p1 + D44.*p4) + D45.*p5;
Dp5 = (D51.*p1 + D54.*p4) + D55.*p5;
Dp6 = (D62.*p2 + D66.*p6) + D67.*p7;
Dp7 = (D72.*p2 + D76.*p6) + D77.*p7;
z = zeros(N,1);
F1 = ((((-e1 + Ei1) + (u1.*wr)/ws) + -q1) + Dp1) + (v1.*wrp)/ws;
F2 = ((((-e2 + Ei2) + (u2.*wr)/ws) + -q2) + Dp2) + (v2.*wrp)/ws;
F3 = ((((-e3 + Ei3) +  z         ) + -q3) + Dp3) + z;
F4 = ((((-e4 + Ei4) +  z         ) + -q4) + Dp4) + z;
F5 = ((( z   + Ei5)              ) + Dp5);
F6 = ((( z   + Ei6)              ) + Dp6);
F7 = ((( z   + Ei7)              ) + Dp7);
I_v(RW+1) = (Jx1 + -e1) + -F1;
I_v(RW+2) = (Jx2 + -e2) + -F2;
I_v(RW+3) = (Jx3 + -e3) + -F3;
I_v(RW+4) = (Jx4 + -e4) + -F4;
I_v(RW+5) =  Jx5        + -F5;
I_v(RW+6) =  Jx6        + -F6;
I_v(RW+7) =  Jx7        + -F7;

% ---- Park transformation, voltage:  vabc = invT*edq0 --------------------
ca=cos(th); cb=cos(th-2*pi/3); cc=cos(th+2*pi/3);
sa=sin(th); sb=sin(th-2*pi/3); sc=sin(th+2*pi/3);
% invT = [c -s 1] ; JinvT = [-s -c 0]
pva = ((-sa).*e1 + (-ca).*e2) + z;
pvb = ((-sb).*e1 + (-cb).*e2) + z;
pvc = ((-sc).*e1 + (-cc).*e2) + z;
all_stampVal(SB+73)= ca;  all_stampVal(SB+74)= cb;  all_stampVal(SB+75)= cc;
all_stampVal(SB+76)=-sa;  all_stampVal(SB+77)=-sb;  all_stampVal(SB+78)=-sc;
all_stampVal(SB+79)=  1;  all_stampVal(SB+80)=  1;  all_stampVal(SB+81)=  1;
all_stampVal(SB+82)=pva;  all_stampVal(SB+83)=pvb;  all_stampVal(SB+84)=pvc;
iTea = (ca.*e1 + (-sa).*e2) + e3;      % invT*e(1:3)
iTeb = (cb.*e1 + (-sb).*e2) + e3;
iTec = (cc.*e1 + (-sc).*e2) + e3;
I_v(RW+ 8) = (((ca.*e1 + (-sa).*e2) + e3) + pva.*th) + (-b2V).*va + -((-va.*b2V) + iTea);
I_v(RW+ 9) = (((cb.*e1 + (-sb).*e2) + e3) + pvb.*th) + (-b2V).*vb + -((-vb.*b2V) + iTeb);
I_v(RW+10) = (((cc.*e1 + (-sc).*e2) + e3) + pvc.*th) + (-b2V).*vc + -((-vc.*b2V) + iTec);

% Row 11 (Efd) is left at zero - the exciter writes it.

% ---- Swing / speed ------------------------------------------------------
% M = Lsp + Lsp.' ; column 3 of i.'*M is structurally zero
M12 = S12 + S21;    % M(1,2) and M(2,1) are both S12+S21
% NOTE: row-vector x matrix and row x column use MATLAB's pairwise
% reduction tree, NOT left-to-right:  acc = t1+(t2+t3); +(t4+t5); +(t6+t7).
% Structural zeros must be kept in mind when grouping - they contribute
% nothing numerically but they do determine which terms pair up.
sw1 = i2.*M12 + (i6.*S16 + i7.*S17);
sw2 = i1.*M12 + (i4.*S24 + i5.*S25);
sw4 = i2.*S24;   sw5 = i2.*S25;   sw6 = i1.*S16;   sw7 = i1.*S17;
Ksw = (4*H/dT+KD)/ws;
all_stampVal(SB+85)=sw1; all_stampVal(SB+86)=sw2; all_stampVal(SB+88)=sw4;
all_stampVal(SB+89)=sw5; all_stampVal(SB+90)=sw6; all_stampVal(SB+91)=sw7;
all_stampVal(SB+92)=Ksw; all_stampVal(SB+93)=-1;
all_stampVal(SB+94)=-dT; all_stampVal(SB+95)= 2;
% (i.')*Lsp*i  is  ((i.')*Lsp)*i
r1 = i2.*S21;  r2 = i1.*S12;  r4 = i2.*S24;
r5 = i2.*S25;  r6 = i1.*S16;  r7 = i1.*S17;
qq  = ((r1.*i1 + r2.*i2) + (r4.*i4 + r5.*i5)) + (r6.*i6 + r7.*i7);
t1 = p2.*S21;  t2 = p1.*S12;  t4 = p2.*S24;
t5 = p2.*S25;  t6 = p1.*S16;  t7 = p1.*S17;
qqp = ((t1.*p1 + t2.*p2) + (t4.*p4 + t5.*p5)) + (t6.*p6 + t7.*p7);
% 9-term row x column, same tree: t1+(t2+t3); +(t4+t5); +(t6+t7); +(t8+t9)
JxSw = (((sw1.*i1 + sw2.*i2) + (sw4.*i4 + sw5.*i5)) + (sw6.*i6 + sw7.*i7)) ...
    + (Ksw.*wr + (-1).*Tm);
Fsw  = (((((Ksw.*wr + qq) + -Tm) + (-4*H/(ws*dT)).*wrp) ...
    + (KD/ws).*(wrp-2*ws)) + qqp) + -Tmp;
I_v(RW+12) = JxSw + -Fsw;
Fsp = (2*th - dT*wr) + (-2*thp - dT*wrp);
I_v(RW+13) = ((-dT)*wr + 2*th) + -Fsp;

% ---- Park transformation, current:  iabc = invT*idq0 -------------------
pia = ((-sa).*i1 + (-ca).*i2) + z;
pib = ((-sb).*i1 + (-cb).*i2) + z;
pic = ((-sc).*i1 + (-cc).*i2) + z;
all_stampVal(SB+ 96)= ca;  all_stampVal(SB+ 97)= cb;  all_stampVal(SB+ 98)= cc;
all_stampVal(SB+ 99)=-sa;  all_stampVal(SB+100)=-sb;  all_stampVal(SB+101)=-sc;
all_stampVal(SB+102)=  1;  all_stampVal(SB+103)=  1;  all_stampVal(SB+104)=  1;
all_stampVal(SB+105)=pia;  all_stampVal(SB+106)=pib;  all_stampVal(SB+107)=pic;
all_stampVal(SB+108)=-b2I; all_stampVal(SB+112)=-b2I; all_stampVal(SB+116)=-b2I;
iTia = (ca.*i1 + (-sa).*i2) + i3;
iTib = (cb.*i1 + (-sb).*i2) + i3;
iTic = (cc.*i1 + (-sc).*i2) + i3;
I_v(RW+14) = ((((ca.*i1 + (-sa).*i2) + i3) + pia.*th) + (-b2I).*ia) + -((-ia.*b2I) + iTia);
I_v(RW+15) = ((((cb.*i1 + (-sb).*i2) + i3) + pib.*th) + (-b2I).*ib) + -((-ib.*b2I) + iTib);
I_v(RW+16) = ((((cc.*i1 + (-sc).*i2) + i3) + pic.*th) + (-b2I).*ic) + -((-ic.*b2I) + iTic);

% Row 17 (Tm) is left at zero - the governor writes it.

% ---- "Outside" index: -eye(3)*busToGenV ---------------------------------
all_stampVal(SB+117)=-b2V; all_stampVal(SB+121)=-b2V; all_stampVal(SB+125)=-b2V;
end
