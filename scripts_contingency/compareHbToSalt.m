function cmp = compareHbToSalt(hb,salt)
%COMPAREHBTOSALT  Does HB's solution agree with SALT's on this case?
%
%   cmp = compareHbToSalt(hb,salt) compares two CONVERGED objects built from
%   the same data and returns
%
%     cmp.worstRel    the largest relative difference found
%     cmp.worstName   which quantity it was on
%     cmp.rel         a struct of the per-quantity relative differences
%     cmp.freqHz      [HB, SALT] fundamental frequency, Hz
%
%   WHAT IS COMPARED, AND WHY NOT MORE. The quantities below are the ones
%   the two solvers both represent AND the ones the contingency verdict is
%   computed from: the fundamental frequency, every bus voltage phasor, and
%   every generating-unit, load and branch current. If these agree, the two
%   models flag the same violations, which is the claim this function
%   exists to support.
%
%   The internal device states are NOT compared here, because the two
%   layouts are not the same: HB carries all seven machine windings where
%   SALT keeps three, every HB quantity is a (real,imag) pair at every
%   harmonic, and matching them up takes an explicit slot-by-slot map. That
%   map, and the structural assertions that go with it (no harmonic content
%   where the Park transforms say there can be none), live in
%   validation/validateHbVsSalt, which checks every state of all five device
%   models. This function is the cheap in-sweep version of that test.
%
%   RELATIVE TO WHAT. Each difference is scaled by the norm of SALT's own
%   value, floored at 1, so a quantity that is legitimately near zero is not
%   held to an impossible relative tolerance.
%
%   See also CALCHBCURRENTS, CALCSALTCURRENTS, VALIDATEHBVSSALT.

assert(hb.numBus==salt.numBus && isequal(hb.bus,salt.bus), ...
    'compareHbToSalt: the two objects are not the same network')

cmp.rel = struct();

% ---- Fundamental frequency ----
w_hb   = hb.x_v(hb.getRowIdx_w0);
w_salt = salt.x_v(salt.getRowIdx_w);
cmp.freqHz = [w_hb,w_salt]/(2*pi);
cmp.rel.f0 = relDiff(w_hb,w_salt);

% ---- Bus voltages ----
% SALT keeps one phasor per node. HB keeps a coefficient per harmonic, and
% the voltage belongs at the fundamental, where the phasor is 2/N times the
% stored coefficient (the 2 is the conjugate split, the 1/N the unnormalized
% DFT convention).
vS = toPhasors(salt.x_v(salt.getRowIdxAbc_fromBusNums(salt.bus)));
vH = toPhasors(hb.x_v(hb.getRowIdxAbc_fromBusNums(hb.bus,1)))*(2/hb.N);
cmp.rel.busV = relDiff(vH,vS);
% ...and the DC bin must be empty: a balanced fundamental set driven through
% the Park transforms produces no DC network voltage, and if it did, the
% comparison above would be hiding content rather than matching it.
vH0 = toPhasors(hb.x_v(hb.getRowIdxAbc_fromBusNums(hb.bus,0)))/hb.N;
cmp.rel.busV_dc = max(abs(vH0));

% ---- Currents ----
[IgH,IlH,IsH] = calcHbCurrents(hb);
[IgS,IlS,IsS] = calcSaltCurrents(salt);
% The bus lists must line up before the values mean anything
assert(isequal(IgH(:,1),IgS(:,1)),'generating-unit bus lists differ')
assert(isequal(IlH(:,1),IlS(:,1)),'load bus lists differ')
assert(isequal(IsH(:,1:2),IsS(:,1:2)),'branch lists differ')
cmp.rel.I_gen    = relDiff(IgH(:,2),IgS(:,2));
cmp.rel.I_load   = relDiff(IlH(:,2),IlS(:,2));
cmp.rel.I_series = relDiff(IsH(:,3:5),IsS(:,3:5));

% ---- Worst of the lot ----
names = fieldnames(cmp.rel);
vals  = cellfun(@(f)cmp.rel.(f),names);
[cmp.worstRel,iw] = max(vals);
cmp.worstName = names{iw};

end


function d = relDiff(got,want)
d = norm(got(:)-want(:))/max(norm(want(:)),1);
end


function c = toPhasors(v)
%TOPHASORS  A column of interleaved (re,im) entries to complex.
c = v(1:2:end) + 1j*v(2:2:end);
end
