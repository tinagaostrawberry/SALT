function [yRow,yCol,yVal,iRow,iVal] = stampOfPfLinearMatrices_SS( ...
    ... % Branches (PI line)
    br_idxVn_r,br_idxVn_i,br_idxVm_r,br_idxVm_i, br_r,br_x,br_b, ...
    ... % Transformers (PI line behind a complex tap)
    xf_idxVn_r,xf_idxVn_i,xf_idxVm_r,xf_idxVm_i, xf_r,xf_x,xf_b, ...
    xf_tap,xf_shift, ...
    ... % Shunts
    sh_idxVr,sh_idxVi, sh_Gs,sh_Bs, ...
    ... % Slack
    slack_idxVr,slack_idxVi, slack_node_Ir,slack_node_Ii, ...
    slack_Vset,slack_AngSet, slack_imScalConst,node_FreqSlack, ...
    ... % System parameters
    useFreqDeviat)
%STAMPOFPFLINEARMATRICES_SS  Constant (voltage-independent) network stamps.
%
%   Branches, transformers, shunts and the slack row are linear, so their
%   triplets are built once and reused for every Newton iteration. The
%   caller is expected to have dropped out-of-service devices already.
%
%   Real and imaginary parts are carried as separate rows and columns, so a
%   scalar admittance G+jB becomes the 2x2 block [G -B; B G]. Everything
%   below is that block written out.
%
%   Duplicate (row,col) pairs are fine: SPARSE sums them, which is exactly
%   what the accumulating per-device stamps used to do into a dense matrix.

% ---------------------------------------------------------------- Branch
br_den = br_r.^2 + br_x.^2;
br_G   = br_r ./ br_den;
br_B   = -br_x ./ br_den;
br_Bsh = br_b / 2;

brRow = [ ...
    br_idxVn_r; br_idxVn_r; br_idxVn_i; br_idxVn_i; ...   % self, from
    br_idxVn_r; br_idxVn_r; br_idxVn_i; br_idxVn_i; ...   % mutual n->m
    br_idxVm_r; br_idxVm_r; br_idxVm_i; br_idxVm_i; ...   % self, to
    br_idxVm_r; br_idxVm_r; br_idxVm_i; br_idxVm_i];      % mutual m->n
brCol = [ ...
    br_idxVn_r; br_idxVn_i; br_idxVn_r; br_idxVn_i; ...
    br_idxVm_r; br_idxVm_i; br_idxVm_r; br_idxVm_i; ...
    br_idxVm_r; br_idxVm_i; br_idxVm_r; br_idxVm_i; ...
    br_idxVn_r; br_idxVn_i; br_idxVn_r; br_idxVn_i];
brVal = [ ...
     br_G; -(br_B + br_Bsh);  (br_B + br_Bsh);  br_G; ...
    -br_G;   br_B;            -br_B;           -br_G; ...
     br_G; -(br_B + br_Bsh);  (br_B + br_Bsh);  br_G; ...
    -br_G;   br_B;            -br_B;           -br_G];

% ----------------------------------------------------------- Transformer
% The tap sits on the "from" side, so the from-side self block is scaled by
% 1/tau^2 and the two mutual blocks are divided by conj(t) and t
% respectively - which is what makes the off-diagonal blocks asymmetric
% whenever the phase shift is nonzero.
xf_den = xf_r.^2 + xf_x.^2;
xf_G   = xf_r ./ xf_den;
xf_B   = -xf_x ./ xf_den;
xf_Bsh = xf_b / 2;

xf_shiftRad = xf_shift * pi / 180;
xf_cr   = xf_tap .* cos(xf_shiftRad);
xf_ci   = xf_tap .* sin(xf_shiftRad);
xf_tau2 = xf_tap.^2;

xf_G_nn = xf_G ./ xf_tau2;
xf_B_nn = (xf_B + xf_Bsh) ./ xf_tau2;
xf_G_mm = xf_G;
xf_B_mm = xf_B + xf_Bsh;
xf_G_nm = (-xf_G.*xf_cr + xf_B.*xf_ci) ./ xf_tau2;
xf_B_nm = (-xf_B.*xf_cr - xf_G.*xf_ci) ./ xf_tau2;
xf_G_mn = (-xf_G.*xf_cr - xf_B.*xf_ci) ./ xf_tau2;
xf_B_mn = (-xf_B.*xf_cr + xf_G.*xf_ci) ./ xf_tau2;

xfRow = [ ...
    xf_idxVn_r; xf_idxVn_r; xf_idxVn_i; xf_idxVn_i; ...
    xf_idxVn_r; xf_idxVn_r; xf_idxVn_i; xf_idxVn_i; ...
    xf_idxVm_r; xf_idxVm_r; xf_idxVm_i; xf_idxVm_i; ...
    xf_idxVm_r; xf_idxVm_r; xf_idxVm_i; xf_idxVm_i];
xfCol = [ ...
    xf_idxVn_r; xf_idxVn_i; xf_idxVn_r; xf_idxVn_i; ...
    xf_idxVm_r; xf_idxVm_i; xf_idxVm_r; xf_idxVm_i; ...
    xf_idxVm_r; xf_idxVm_i; xf_idxVm_r; xf_idxVm_i; ...
    xf_idxVn_r; xf_idxVn_i; xf_idxVn_r; xf_idxVn_i];
xfVal = [ ...
    xf_G_nn; -xf_B_nn; xf_B_nn; xf_G_nn; ...
    xf_G_nm; -xf_B_nm; xf_B_nm; xf_G_nm; ...
    xf_G_mm; -xf_B_mm; xf_B_mm; xf_G_mm; ...
    xf_G_mn; -xf_B_mn; xf_B_mn; xf_G_mn];

% ----------------------------------------------------------------- Shunt
shRow = [sh_idxVr; sh_idxVr; sh_idxVi; sh_idxVi];
shCol = [sh_idxVr; sh_idxVi; sh_idxVr; sh_idxVi];
shVal = [ sh_Gs;   -sh_Bs;    sh_Bs;    sh_Gs];

% ----------------------------------------------------------------- Slack
% NOTE: in the original Slack.stamp_linear these were ASSIGNMENTS (=), not
% accumulations (+=), unlike every other device. That is safe to reproduce
% as ordinary triplets because no other device ever writes these cells: the
% slack-current rows/columns and the frequency-slack row belong to the slack
% alone. Emitting them as triplets therefore leaves Y unchanged, while
% keeping the "this row is defined, not accumulated" intent explicit here.
if useFreqDeviat
    % Pin the slack angle by fixing the imag/real ratio of its voltage,
    % which is what makes the extra frequency unknown well posed:
    %   0 = imScalConst*Vr - Vi
    slRow = [node_FreqSlack; node_FreqSlack];
    slCol = [slack_idxVr;    slack_idxVi];
    slVal = [slack_imScalConst; -1];

    % I(node_FreqSlack) = 0 - carried explicitly so the row is documented;
    % a zero entry adds nothing to the accumulated RHS.
    slRow_I = node_FreqSlack;
    slVal_I = 0;
else
    % The slack voltage is imposed outright and the current it injects
    % becomes the unknown, so KCL at the bus gains a -1 on that current and
    % two new rows pin Vr and Vi to their set points.
    Vr_set = slack_Vset * cosd(slack_AngSet);
    Vi_set = slack_Vset * sind(slack_AngSet);

    slRow = [slack_idxVr;   slack_idxVi;   slack_node_Ir; slack_node_Ii];
    slCol = [slack_node_Ir; slack_node_Ii; slack_idxVr;   slack_idxVi];
    slVal = [-1;            -1;             1;             1];

    slRow_I = [slack_node_Ir; slack_node_Ii];
    slVal_I = [Vr_set;        Vi_set];
end

% ------------------------------------------------------------- Assemble
yRow = [brRow; xfRow; shRow; slRow];
yCol = [brCol; xfCol; shCol; slCol];
yVal = [brVal; xfVal; shVal; slVal];

iRow = slRow_I;
iVal = slVal_I;
end
