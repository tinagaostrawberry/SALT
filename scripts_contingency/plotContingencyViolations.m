%% Per-contingency violation maps: SALT | Enh-PF | PF
%
% Called from inside the outage loop of run_SALT_contingency.m, after the
% three checkContingencyViolations calls. Draws the network three times for
% the current contingency and highlights that model's violations:
%   red edge      - rated-line violation
%   red node      - bus voltage violation
%   red ring      - generator Q-limit violation
% Node shape marks the bus role (generator circle / load square / other
% dot). All panels share one cached layout so they can be compared
% directly.
%
% A NEW FIGURE IS CREATED FOR EVERY CONTINGENCY, so the sweep leaves one
% window per outage to page through rather than a single window showing
% only the last one.

%% ---- cached geometry (built once, reused every contingency) ----
% Only the LAYOUT is cached. The force-directed solve is the expensive part
% and it does not depend on the contingency, so recomputing it per outage
% would also move every node between figures and make them incomparable.
if ~exist('vizG','var') || isempty(vizG)
    vizG = graph(mpc_0.branch(:,1),mpc_0.branch(:,2));
    vizFtmp = figure('Visible','off');
    vizHtmp = plot(axes(vizFtmp),vizG,'Layout','force');
    vizXY = [vizHtmp.XData(:), vizHtmp.YData(:)];
    close(vizFtmp)
end

%% ---- bus roles for this contingency ----
vizGenBus  = I_gen_contingency_SALT(:,1);
vizLoadBus = setdiff(I_load_contingency_SALT(:,1),vizGenBus);
vizMarker  = repmat({'.'},numnodes(vizG),1);
vizMarker(vizGenBus)  = {'o'};
vizMarker(vizLoadBus) = {'s'};
vizSize = 5*ones(numnodes(vizG),1);
vizSize(vizLoadBus) = 8;
vizSize(vizGenBus)  = 9;

%% ---- what to draw in each panel ----
vizToFT = @(li) sort(reshape(bus(li),[],2),2);
vizPanel = { ...
    'SALT',   busV_salt(:),  busQ_salt(:),  vizToFT(lineIdx_salt),  converge_SS; ...
    'Enh-PF', busV_enhPf(:), busQ_enhPf(:), vizToFT(lineIdx_enhPf), enhPfResult.success; ...
    'PF',     busV_pf(:),    busQ_pf(:),    vizToFT(lineIdx_pf),    pfResult.success};

%% ---- figure (one per contingency) ----
vizFig = figure('Color','w', ...
    'Name',sprintf('Violations by model - contingency %d',idxOutage));
set(vizFig,'Position',[40 80 1900 760]);
vizTL = tiledlayout(vizFig,1,3,'TileSpacing','compact','Padding','compact');

vizHdr = sprintf('Contingency %d of %d  (generator outage', ...
    idxOutage,numContingencies);
if exist('outagedGenBus','var') && ~isnan(outagedGenBus(idxOutage))
    vizHdr = sprintf('%s at bus %d)',vizHdr,outagedGenBus(idxOutage));
else
    vizHdr = [vizHdr ')'];
end
title(vizTL,vizHdr,'FontWeight','bold')
subtitle(vizTL,['red edge = rated-line violation,  red node = voltage violation,  ' ...
    'red ring = generator Q-limit violation'])

for vizP = 1:3
    vizAx = nexttile(vizTL);
    vizH = plot(vizAx,vizG,'XData',vizXY(:,1),'YData',vizXY(:,2));
    vizH.NodeColor    = [0.35 0.35 0.35];
    vizH.EdgeColor    = [0.78 0.78 0.78];
    vizH.EdgeAlpha    = 1;
    vizH.LineWidth    = 1.5;
    vizH.Marker       = vizMarker;
    vizSizeLocal      = vizSize;   % enlarged below for violated buses
    vizH.MarkerSize   = vizSizeLocal;
    vizH.NodeLabel    = 1:numnodes(vizG);
    vizH.NodeFontSize = 11;

    vizName = vizPanel{vizP,1};
    vizV = vizPanel{vizP,2}; vizQ = vizPanel{vizP,3};
    vizL = vizPanel{vizP,4}; vizOK = vizPanel{vizP,5};

    if ~vizOK
        title(vizAx,sprintf('%s  --  no solution',vizName),'Color',[0.6 0 0])
    else
        % rated-line violations -> red edges
        if ~isempty(vizL)
            vizE = findedge(vizG,vizL(:,1),vizL(:,2));
            vizE = vizE(vizE>0);
            if ~isempty(vizE)
                highlight(vizH,'Edges',vizE,'EdgeColor',[0.85 0 0],'LineWidth',5)
            end
        end
        % voltage violations -> red nodes (node IDs are positional)
        if ~isempty(vizV)
            vizSizeLocal(vizV) = max(vizSizeLocal(vizV),11);
            vizH.MarkerSize = vizSizeLocal;
            highlight(vizH,vizV,'NodeColor',[0.85 0 0])
        end
        % Q-limit violations -> red ring around the generator
        if ~isempty(vizQ)
            hold(vizAx,'on')
            scatter(vizAx,vizXY(vizQ,1),vizXY(vizQ,2),320,'o', ...
                'MarkerEdgeColor',[0.85 0 0],'LineWidth',2.5);
            hold(vizAx,'off')
        end
        title(vizAx,sprintf('%s  --  V:%d  Q:%d  Line:%d', ...
            vizName,numel(vizV),numel(vizQ),size(vizL,1)))
    end
    axis(vizAx,'equal'); axis(vizAx,'off');
end
set(findall(vizFig,'-property','FontSize'),'FontSize',14);
drawnow
