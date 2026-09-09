function [obj,v_init] = pfValInitializations(obj)
%PFVALINITIALIZATIONS  Build the Newton starting point and the ZIP shares.
%
%   [obj,v_init] = pfValInitializations(pf)
%
%   Three things happen here, in the order the original initialize.m used,
%   because the later steps read the values the earlier ones wrote:
%
%     1. Bus voltages are seeded from the case's own Vm/Va (a "flat start"
%        only in the sense that no network solve has happened yet).
%     2. PV generators seed their reactive-power unknown, and ZIP loads
%        turn their (Z,I,P) weights into coefficients normalised at the
%        seeded voltage - which is why the coefficients cannot be derived
%        in the constructor.
%     3. The slack is seeded last: with frequency fixed there is nothing to
%        do but set w, and otherwise its injected-current unknowns are
%        seeded by summing the current every incident shunt, branch and
%        transformer draws at the seeded voltages. Starting the slack
%        current at zero instead would put the first Newton step far from
%        the solution on any nontrivially loaded case.

v_init = zeros(obj.numNodes,1);

%% 1. Bus voltages (Va is in degrees, as MATPOWER stores it)
v_init(obj.bus_node_Vr) = obj.bus_Vm .* cosd(obj.bus_Va);
v_init(obj.bus_node_Vi) = obj.bus_Vm .* sind(obj.bus_Va);

%% 2a. Generator reactive power
% Q is stored with the load sign convention (see PF.setGenParams), so the
% seed is the negated case value.
if obj.numGen > 0
    idxBusGen = obj.getBusIdx(obj.gen_bus);
    v_init(obj.bus_node_Q(idxBusGen)) = obj.gen_Qinit;
end

%% 2b. ZIP coefficients
% The weights are shares of the load at its nominal voltage, so dividing by
% Vm^2 / Vm turns them into the coefficients of a polynomial in |V| that
% reproduces exactly P (and Q) when |V| == Vm.
obj.load_Pz_IM  = zeros(obj.numLoad,1);
obj.load_Pi_IM  = zeros(obj.numLoad,1);
obj.load_Pp_IM  = zeros(obj.numLoad,1);
obj.load_Qz_IM  = zeros(obj.numLoad,1);
obj.load_Qi_IM  = zeros(obj.numLoad,1);
obj.load_Qp_IM  = zeros(obj.numLoad,1);
obj.load_Pz_ZIP = zeros(obj.numLoad,1);
obj.load_Pi_ZIP = zeros(obj.numLoad,1);
obj.load_Pp_ZIP = zeros(obj.numLoad,1);
obj.load_Qz_ZIP = zeros(obj.numLoad,1);
obj.load_Qi_ZIP = zeros(obj.numLoad,1);
obj.load_Qp_ZIP = zeros(obj.numLoad,1);

idxZip = find(obj.load_UseZip);
if ~isempty(idxZip)
    Vm = obj.bus_Vm(obj.getBusIdx(obj.load_bus(idxZip)));
    P  = obj.load_P(idxZip);
    Q  = obj.load_Q(idxZip);

    obj.load_Pz_IM(idxZip)  = P.*obj.load_pw_IM(idxZip,1)./(Vm.^2);
    obj.load_Pi_IM(idxZip)  = P.*obj.load_pw_IM(idxZip,2)./Vm;
    obj.load_Pp_IM(idxZip)  = P.*obj.load_pw_IM(idxZip,3);
    obj.load_Qz_IM(idxZip)  = Q.*obj.load_qw_IM(idxZip,1)./(Vm.^2);
    obj.load_Qi_IM(idxZip)  = Q.*obj.load_qw_IM(idxZip,2)./Vm;
    obj.load_Qp_IM(idxZip)  = Q.*obj.load_qw_IM(idxZip,3);

    obj.load_Pz_ZIP(idxZip) = P.*obj.load_pw_ZIP(idxZip,1)./(Vm.^2);
    obj.load_Pi_ZIP(idxZip) = P.*obj.load_pw_ZIP(idxZip,2)./Vm;
    obj.load_Pp_ZIP(idxZip) = P.*obj.load_pw_ZIP(idxZip,3);
    obj.load_Qz_ZIP(idxZip) = Q.*obj.load_qw_ZIP(idxZip,1)./(Vm.^2);
    obj.load_Qi_ZIP(idxZip) = Q.*obj.load_qw_ZIP(idxZip,2)./Vm;
    obj.load_Qp_ZIP(idxZip) = Q.*obj.load_qw_ZIP(idxZip,3);
end

%% 3. Slack, last
if obj.useFreqDeviat
    v_init(obj.node_FreqSlack) = obj.w0;
    return
end

idxSlackBus = obj.getBusIdx(obj.slack_bus);
Vbus = v_init(obj.bus_node_Vr) + 1i*v_init(obj.bus_node_Vi);
Vs   = Vbus(idxSlackBus);

I_slack = 0 + 0i;

% --- Shunts at the slack bus -------------------------------------------
mSh = obj.shunt_status & (obj.shunt_bus == obj.slack_bus);
if any(mSh)
    Y_sh = obj.shunt_Gs(mSh) + 1i*obj.shunt_Bs(mSh);
    I_slack = I_slack + sum(Y_sh * Vs);
end

% --- Branches incident on the slack bus --------------------------------
% A branch that starts and ends at the slack bus is counted once, on its
% "from" side, matching the original if/else.
mFrom = obj.branch_status & (obj.branch_busFrom == obj.slack_bus);
mTo   = obj.branch_status & (obj.branch_busTo   == obj.slack_bus) & ~mFrom;
mInc  = mFrom | mTo;
if any(mInc)
    busOther = obj.branch_busTo;
    busOther(mTo) = obj.branch_busFrom(mTo);
    Vm_other = Vbus(obj.getBusIdx(busOther(mInc)));

    y_series = 1 ./ (obj.branch_r(mInc) + 1i*obj.branch_x(mInc));
    b_shunt  = 1i*(obj.branch_b(mInc)/2);

    % KCL: I = (y_series + b_shunt)*V_slack - y_series*V_neighbour
    I_slack = I_slack + sum((y_series + b_shunt)*Vs - y_series.*Vm_other);
end

% --- Transformers incident on the slack bus ----------------------------
% Asymmetric: which side the slack sits on decides whether the tap scaling
% applies to the self term.
xfFrom = obj.xfmr_status & (obj.xfmr_busFrom == obj.slack_bus);
xfTo   = obj.xfmr_status & (obj.xfmr_busTo   == obj.slack_bus) & ~xfFrom;

if any(xfFrom)
    y_series = 1 ./ (obj.xfmr_r(xfFrom) + 1i*obj.xfmr_x(xfFrom));
    b_shunt  = 1i*(obj.xfmr_b(xfFrom)/2);
    t = obj.xfmr_tap(xfFrom) .* exp(1i*obj.xfmr_shift(xfFrom)*pi/180);
    Vm_other = Vbus(obj.getBusIdx(obj.xfmr_busTo(xfFrom)));

    Y_ff = (y_series + b_shunt) ./ (abs(t).^2);
    Y_ft = -y_series ./ conj(t);
    I_slack = I_slack + sum(Y_ff*Vs + Y_ft.*Vm_other);
end
if any(xfTo)
    y_series = 1 ./ (obj.xfmr_r(xfTo) + 1i*obj.xfmr_x(xfTo));
    b_shunt  = 1i*(obj.xfmr_b(xfTo)/2);
    t = obj.xfmr_tap(xfTo) .* exp(1i*obj.xfmr_shift(xfTo)*pi/180);
    Vm_other = Vbus(obj.getBusIdx(obj.xfmr_busFrom(xfTo)));

    Y_tt = y_series + b_shunt;
    Y_tf = -y_series ./ t;
    I_slack = I_slack + sum(Y_tt*Vs + Y_tf.*Vm_other);
end

v_init(obj.slack_node_Ir) = real(I_slack);
v_init(obj.slack_node_Ii) = imag(I_slack);
end
