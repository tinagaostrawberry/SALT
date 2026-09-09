function [mpc,new_slack_bus_id] = setGenContingency(mpc, busNum_genLoss)
%SETGENCONTINGENCY Trip every generator at the given bus(es), demote those
% buses to PQ, and reassign the slack if it was among them.
%
%   [mpc,new_slack_bus_id] = setGenContingency(mpc, busNum_genLoss)
%
% busNum_genLoss may be a SCALAR or a VECTOR of bus numbers. Each bus is
% handled independently, so a mixed set (slack + PV, say) is fine.
%
% new_slack_bus_id is [] when the slack was NOT tripped, and the bus number
% of the replacement slack when it was. Callers rely on that convention to
% decide whether to re-point their own slack (e.g. emtSS.slack_busNum), so
% it must never come back empty after the slack has actually moved.
%
% WHY PQ AND NOT "ISOLATED": tripping a generator removes an INJECTION, not
% the node. The bus keeps its branches, its load and its shunt, so it
% becomes an ordinary PQ bus (possibly zero-injection, which is fine).
% Marking it isolated (BUS_TYPE 4) would delete the bus AND every branch
% attached to it - that would sever a chain-junction bus and island whole
% sections of a concatenated system. Note MATPOWER's bustypes() already
% treats any bus with no online generator as PQ regardless of BUS_TYPE
% ("pq = find(bus(:,BUS_TYPE)==PQ | ~bus_gen_status)"); setting the type
% explicitly just keeps mpc consistent for the solvers that read BUS_TYPE
% directly instead of going through bustypes().

define_constants;

new_slack_bus_id = [];

busNum_genLoss = busNum_genLoss(:);
if isempty(busNum_genLoss)
    return;
end

%% 1. Locate the generators
% Use a LOGICAL MASK over gen rows, not ismember(busNum_genLoss,...) - that
% form returns only the FIRST matching gen row per bus, so a bus with two
% units would have had just one of them tripped while still being demoted
% to PQ.
genMask = ismember(mpc.gen(:, GEN_BUS), busNum_genLoss);

% Warn about (and drop) bus numbers that carry no generator at all. The old
% "isempty(gen_rows)" test could never fire, because ismember returns 0 for
% a non-match rather than [], which then indexed mpc.gen(0,:) and errored.
hasGen = ismember(busNum_genLoss, mpc.gen(:, GEN_BUS));
if ~all(hasGen)
    warning('setGenContingency:noGenAtBus', ...
        'No generator at bus %s. Ignoring.', ...
        mat2str(busNum_genLoss(~hasGen).'));
    busNum_genLoss = busNum_genLoss(hasGen);
    if isempty(busNum_genLoss)
        return;
    end
end

%% 2. Trip the generators
mpc.gen(genMask, GEN_STATUS) = 0;

%% 3. Reconfigure bus types - ONE BUS AT A TIME
% The old code tested the whole vector at once ("if current_bus_type == 3"
% is all(...==3) for a vector), so a mixed slack+PV set matched NEITHER
% branch and silently did nothing: the slack kept BUS_TYPE 3 with no online
% generator and new_slack_bus_id came back []. MATPOWER would then quietly
% pick its own replacement slack internally, leaving PF and SALT solving on
% different slack buses with no error raised.
slackWasTripped = false;
for k = 1:numel(busNum_genLoss)
    b = busNum_genLoss(k);
    bus_row = find(mpc.bus(:, BUS_I) == b, 1);
    if isempty(bus_row)
        warning('setGenContingency:noSuchBus', ...
            'Bus %g is not in mpc.bus. Ignoring.', b);
        continue;
    end

    % Only demote once EVERY generator at this bus is off. Since the input
    % is a BUS number, step 2 trips all units at that bus, so this holds by
    % construction - it is kept as a defensive check (and to make the
    % invariant explicit) rather than because it can currently fire.
    if any(mpc.gen(mpc.gen(:, GEN_BUS) == b, GEN_STATUS) > 0)
        continue;
    end

    switch mpc.bus(bus_row, BUS_TYPE)
        case REF
            mpc.bus(bus_row, BUS_TYPE) = PQ;
            slackWasTripped = true;
        case PV
            mpc.bus(bus_row, BUS_TYPE) = PQ;
    end
end

%% 4. Reassign the slack, if we lost it
% Done AFTER the loop so the choice is made against the final post-trip
% state - picking mid-loop could hand the slack to a machine that a later
% iteration then trips.
if slackWasTripped
    online_gens_idx = find(mpc.gen(:, GEN_STATUS) > 0);
    if isempty(online_gens_idx)
        error('setGenContingency:blackout', ...
            'SYSTEM BLACKOUT: No active generators remain.');
    end

    % Highest-capacity online machine.
    [~,max_p_idx]    = max(mpc.gen(online_gens_idx, PMAX));
    best_gen_row     = online_gens_idx(max_p_idx);
    new_slack_bus_id = mpc.gen(best_gen_row, GEN_BUS);

    new_slack_row = find(mpc.bus(:, BUS_I) == new_slack_bus_id, 1);
    mpc.bus(new_slack_row, BUS_TYPE) = REF;

    % The new slack needs a sane voltage setpoint. The old code threw a
    % bare error('Debugging') here; warn and fall back to 1.0 pu instead so
    % a contingency sweep is not aborted by one odd starting voltage.
    vmNew = mpc.bus(new_slack_row, VM);
    if vmNew < 0.9 || vmNew > 1.1
        warning('setGenContingency:slackVm', ...
            ['New slack bus %g had Vm = %.4f pu, outside [0.9,1.1]. ' ...
             'Resetting to 1.0 pu.'], new_slack_bus_id, vmNew);
        mpc.bus(new_slack_row, VM) = 1.0;
    end
end

%% 5. Clear stale internal ordering
% MATPOWER caches internal indexing in mpc.order. Because BUS_TYPE changed,
% the system dimensions changed too; leaving mpc.order in place makes runpf
% use the old map and crash while extracting results.
if isfield(mpc, 'order')
    mpc = rmfield(mpc, 'order');
end

end
