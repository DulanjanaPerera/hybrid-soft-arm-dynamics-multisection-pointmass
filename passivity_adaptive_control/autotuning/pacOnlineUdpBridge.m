function rx=pacOnlineUdpBridge(tx,blockPort,workerPort,session)
%PACONLINEUDPBRIDGE Normal-simulation, loopback-only nonblocking transport.
% Called extrinsically. Never run optimization from the simulation thread.
persistent channel buffer lastRx owner cleanup warnedSession
if isempty(owner)||owner~=session
    try
    if ~isempty(channel),channel.close();end
    channel=java.nio.channels.DatagramChannel.open();
    channel.bind(java.net.InetSocketAddress('127.0.0.1',int32(blockPort)));
    channel.configureBlocking(false);
    buffer=java.nio.ByteBuffer.allocate(16384);
    localChannel=channel;cleanup=onCleanup(@()localChannel.close()); %#ok<NASGU>
    owner=session;lastRx=zeros(7,1);
    catch err
        if ~isempty(channel)
            try,channel.close();catch,end
        end
        channel=[];
        owner=[];rx=[session;0;0;-1;0;0;0];
        if isempty(warnedSession)||warnedSession~=session
            warning('pacTune:UdpBind', ...
                ['Cannot initialize model UDP port %d: %s. Stop the worker, ' ...
                 'clear pacOnlineUdpBridge in its MATLAB process, then restart ' ...
                 'the worker and run pacOnlineConnectionCheck in model MATLAB.'],blockPort,err.message);
            warnedSession=session;
        end
        return;
    end
end
rx=lastRx;
try
    b=java.nio.ByteBuffer.wrap(int8(sprintf('%.17g,',tx)));
    channel.send(b,java.net.InetSocketAddress('127.0.0.1',int32(workerPort)));
    for k=1:8 % Bound work per simulation sample; never block on receive.
        buffer.clear();from=channel.receive(buffer);
        if isempty(from),break;end
        if from.getPort()~=workerPort,continue;end
        n=double(buffer.position());raw=buffer.array();
        p=sscanf(char(raw(1:n))','%f,');
        if numel(p)==7&&all(isfinite(p))&&p(1)==session&&p(2)>lastRx(2)
            lastRx=p(:);
        end
    end
    rx=lastRx;
catch
    % Retain the old sequence. Runtime watchdog notices lack of progress.
end
end
