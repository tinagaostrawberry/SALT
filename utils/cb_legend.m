function cb_legend(~,evt)
if strcmp(evt.Peer.Visible,'on')
    evt.Peer.Visible = 'off';
else
    evt.Peer.Visible = 'on';
end
end