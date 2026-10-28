function pacOnlineConnectionCheck(configFile)
%PACONLINECONNECTIONCHECK Model-process UDP preflight; no model/hardware run.
% Worker must already be running in the OTHER MATLAB process.
d=load(configFile,'rt');c=d.rt;
diagrams=find_system('type','block_diagram');
for k=1:numel(diagrams)
    if strcmp(get_param(diagrams{k},'BlockDiagramType'),'model')
        assert(strcmp(get_param(diagrams{k},'SimulationStatus'),'stopped'), ...
            'Stop simulation before connection preflight.');
    end
end
clear pacOnlineUdpBridge
cleanup=onCleanup(@resetBridge); %#ok<NASGU>
p=zeros(121,1);p(1)=c.session;p(2)=1;p(4)=-1;p(5)=-1;
clockId=tic;
while toc(clockId)<3
    rx=pacOnlineUdpBridge(p,c.blockPort,c.workerPort,c.session);
    assert(rx(4)~=-1,'pacTune:UdpPortBusy', ...
        'Model port %d is unavailable. Clear pacOnlineUdpBridge in the OTHER MATLAB process before starting its worker.',c.blockPort);
    if rx(1)==c.session&&rx(2)>0&&rx(3)==1&&rx(4)==1&& ...
            all(abs(rx(5:7)-c.baselineGamma)<1e-12)
        fprintf('PASS: model port %d <-> worker port %d, matching session %.17g and baseline. No pressure/model run.\n', ...
            c.blockPort,c.workerPort,c.session);return;
    end
    pause(.05);
end
error('pacTune:WorkerConnection', ...
    'No matching worker reply. Start the matching worker in the OTHER MATLAB process, then retry before enabling pressure.');
end
function resetBridge()
clear pacOnlineUdpBridge
end
