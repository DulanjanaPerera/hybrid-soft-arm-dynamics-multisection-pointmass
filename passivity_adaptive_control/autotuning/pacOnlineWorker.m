function pacOnlineWorker(configFile)
%PACONLINEWORKER Launch in a SEPARATE MATLAB process, before starting model.
% Only localhost UDP and logs; never connects to sensors/valves or starts SLX.
d=load(configFile,'rt','c','snapshot');rt=d.rt;c=d.c;
folder=fileparts(configFile);logFile=fullfile(folder,'worker_log.mat');
assert(~isfile(logFile),'Worker logs exist. Prepare a new session, do not overwrite.');
assert(isempty(dir(fullfile(folder,'batch_*.mat'))),'Batch checkpoints exist; use a new session.');
lock=fullfile(folder,'worker.lock');
f=java.io.File(lock);assert(f.createNewFile(),'A worker already owns this session.');
release=onCleanup(@()delete(lock)); %#ok<NASGU>
% This MATLAB instance may previously have run a model. Its persistent bridge
% otherwise retains the model port and silently blocks the other instance.
diagrams=find_system('type','block_diagram');
for k=1:numel(diagrams)
    if strcmp(get_param(diagrams{k},'BlockDiagramType'),'model')
        assert(strcmp(get_param(diagrams{k},'SimulationStatus'),'stopped'), ...
            'pacTune:WorkerProcess','Run worker in a separate MATLAB process with no active model.');
    end
end
clear pacOnlineUdpBridge
channel=java.nio.channels.DatagramChannel.open();
closeChannel=onCleanup(@()channel.close()); %#ok<NASGU>
channel.bind(java.net.InetSocketAddress('127.0.0.1',int32(rt.workerPort)));
channel.configureBlocking(false);buffer=java.nio.ByteBuffer.allocate(16384);
s=pacOnlineWorkerStep([],zeros(16,1),c,rt);
connectTimeout_s=1800; % Preflight/first compilation may take several minutes.
if isfield(c,'workerConnectTimeout_s'),connectTimeout_s=c.workerConnectTimeout_s;end
packets=zeros(60000,121);count=0;sequence=0;clockId=tic;lastRx=-1;lastSend=-1;
lastTime=-1;oldCompleted=0;failure='';terminal=false;
runtimeReason=0;lastPacket=zeros(121,1);
fprintf('Worker session %d listening on localhost:%d. Start the model manually.\n',rt.session,rt.workerPort);
while toc(clockId)<3600
    now=toc(clockId);
    for k=1:64
        buffer.clear();address=channel.receive(buffer);
        if isempty(address),break;end
        if address.getPort()~=rt.blockPort,continue;end
        n=double(buffer.position());raw=buffer.array();p=sscanf(char(raw(1:n))','%f,');
        if numel(p)~=121||p(1)~=rt.session||~isfinite(p(5)),continue;end
        if p(4)==-1 % Preflight probe: reply without starting telemetry timeout.
            sequence=sequence+1;reply=[rt.session;sequence;s.requestId;s.action;s.gamma];
            bytes=java.nio.ByteBuffer.wrap(int8(sprintf('%.17g,',reply)));
            channel.send(bytes,java.net.InetSocketAddress('127.0.0.1',int32(rt.blockPort)));
            continue;
        end
        if p(5)<0,continue;end
        lastRx=now;
        lastPacket=p;runtimeReason=p(6);
        if runtimeReason>0&&runtimeReason<20&&s.action==1
            s.action=4;s.reason=['runtime_' num2str(runtimeReason)];
        end
        if p(5)>lastTime
            if count>=size(packets,1),s.action=4;s.reason='log_capacity';
            else,count=count+1;packets(count,:)=p';end
            lastTime=p(5);
        end
        if s.action==1
            try,s=pacOnlineWorkerStep(s,p,c,rt);
            catch err,s.action=4;s.reason='optimizer_exception';failure=getReport(err);end
        end
        if numel(s.attempts)~=oldCompleted
            oldCompleted=numel(s.attempts);
            checkpoint=struct('state',s,'rt',rt,'c',c,'packetCount',count); %#ok<NASGU>
            checkpoint.rawPackets=packets(1:count,:);
            save(fullfile(folder,sprintf('batch_%03d.mat',oldCompleted)),'checkpoint');
        end
        terminal=p(4)==10;
    end
    if lastRx>=0&&now-lastSend>=.1
        sequence=sequence+1;
        reply=[rt.session;sequence;s.requestId;s.action;s.gamma];
        bytes=java.nio.ByteBuffer.wrap(int8(sprintf('%.17g,',reply)));
        channel.send(bytes,java.net.InetSocketAddress('127.0.0.1',int32(rt.blockPort)));
        lastSend=now;
    end
    if terminal,break;end
    if lastRx>=0&&now-lastRx>10,failure='telemetry_timeout';break;end
    if lastRx<0&&now>connectTimeout_s,failure='model_never_connected';break;end
    pause(.01);
end
packets=packets(1:count,:);save(logFile,'s','rt','c','packets','failure','terminal', ...
    'runtimeReason','lastPacket','-v7.3');
fprintf('Worker stopped: %s; %d accepted observations, %d completed-motion results, %d packets.\n', ...
    s.reason,height(s.observations),numel(s.attempts),count);
end
