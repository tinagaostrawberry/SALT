function [bus_node_Vr,bus_node_Vi,bus_node_Q, ...
    slack_node_Ir,slack_node_Ii,node_FreqSlack,numNodes] = ...
    calcPfNodeIndices(bus_type,useFreqDeviat)
%CALCPFNODEINDICES  Lay out the power-flow unknown vector.
%
%   Every bus contributes a real and an imaginary voltage unknown. A PV bus
%   adds a reactive-power unknown, because its voltage magnitude is pinned
%   instead of its Q injection. The slack bus is the odd one out and depends
%   on whether frequency is a variable:
%
%     useFreqDeviat == false  the slack voltage is fixed outright, so the
%                             two slack equations are traded for two
%                             injected-current unknowns appended after all
%                             the bus unknowns.
%     useFreqDeviat == true   the slack machine is dispatched like any other
%                             generator, so the slack bus gets a Q unknown,
%                             and one extra global unknown (the system
%                             frequency) is pinned by a slack-angle equation.
%
%   The walk is deliberately sequential rather than vectorised: the original
%   Buses/Slack classes handed out node numbers one at a time from a
%   persistent counter, and the resulting interleaving (Vr,Vi,[Q],[w] per
%   bus) is what every stamp index below is built against.
%
%   bus_node_Q is 0 for buses that carry no reactive-power unknown, and
%   node_FreqSlack / slack_node_I* are [] when they do not exist.

numBus = numel(bus_type);

bus_node_Vr = zeros(numBus,1);
bus_node_Vi = zeros(numBus,1);
bus_node_Q  = zeros(numBus,1);
node_FreqSlack = [];

nodeCount = 0;
for i = 1:numBus
    nodeCount = nodeCount + 1;  bus_node_Vr(i) = nodeCount;
    nodeCount = nodeCount + 1;  bus_node_Vi(i) = nodeCount;

    switch bus_type(i)
        case 2                                  % PV
            nodeCount = nodeCount + 1;  bus_node_Q(i) = nodeCount;
        case 3                                  % slack / reference
            if useFreqDeviat
                nodeCount = nodeCount + 1;  bus_node_Q(i) = nodeCount;
                nodeCount = nodeCount + 1;  node_FreqSlack = nodeCount;
            end
        otherwise
            % PQ (1) and isolated (4) buses carry voltage unknowns only.
    end
end

% The slack current unknowns come after every bus unknown, because the
% original assigned them in a second pass over the slack list.
if useFreqDeviat
    slack_node_Ir = [];
    slack_node_Ii = [];
else
    nodeCount = nodeCount + 1;  slack_node_Ir = nodeCount;
    nodeCount = nodeCount + 1;  slack_node_Ii = nodeCount;
end

numNodes = nodeCount;
end
