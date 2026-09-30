classdef AuroraReader_s2_Position_only < matlab.System
    % Outputs 1:8: base/tip XYZ (mm), indicators; then validity, frames,
    % host completion time (s), deadline overrun, and poll-start time (s).
    % One BX per step.
    properties (Nontunable)
        Port = 'COM12'
        BaudRate = 921600
        SamplePeriod = 0.05
        BaseSensorId = '0A'
        TipSensorId = '0B'
    end
    properties(Access = private)
        aurora
        clockStart
        nextPoll = 0
        previousFrame = [NaN NaN]
    end
    methods
        function obj = AuroraReader_s2_Position_only(varargin)
            cfg = ndiAcquisitionSettings();
            obj.Port=cfg.Port; obj.BaudRate=cfg.BaudRate;
            obj.SamplePeriod=cfg.SamplePeriod;
            obj.BaseSensorId=cfg.BaseSensorId; obj.TipSensorId=cfg.TipSensorId;
            setProperties(obj,nargin,varargin{:});
        end
    end
    methods(Access = protected)
        function validatePropertiesImpl(obj)
            validateattributes(obj.SamplePeriod,{'double'},{'scalar','finite','positive'});
        end
        function setupImpl(obj)
            obj.aurora = AuroraDriver_2(obj.Port);
            try
                obj.aurora.openSerialPort();
                obj.aurora.init();
                obj.aurora.setBaudRate(obj.BaudRate);
                obj.aurora.detectAndAssignPortHandles();
                obj.aurora.sensorIndices({obj.BaseSensorId,obj.TipSensorId});
                obj.aurora.initPortHandleAll();
                obj.aurora.enablePortHandleDynamicAll();
                obj.aurora.startTracking();
                pause(0.1);
                obj.clockStart = [];
                obj.nextPoll = 0;
                obj.previousFrame = [NaN NaN];
            catch e
                releaseImpl(obj);
                rethrow(e);
            end
        end
        function [x1,y1,z1,e1,x2,y2,z2,e2,valid,frame,hostTime,overrun,pollStart] = stepImpl(obj)
            % Pace physical polls, not just simulated time. Never burst to catch up.
            if isempty(obj.clockStart), obj.clockStart=tic; end
            remaining = obj.nextPoll-toc(obj.clockStart);
            if remaining > 0, pause(remaining); end
            pollStart = toc(obj.clockStart);
            [x,y,z,ep,valid,frame] = obj.aurora.measureTipPositionAll( ...
                {obj.BaseSensorId,obj.TipSensorId});
            hostTime = toc(obj.clockStart);
            [p,valid] = ndiValidSample([x;y;z],valid,frame,obj.previousFrame);
            obj.previousFrame = frame;
            overrun = hostTime > obj.nextPoll+obj.SamplePeriod;
            obj.nextPoll = max(obj.nextPoll+obj.SamplePeriod,pollStart+obj.SamplePeriod);
            x1=p(1,1); y1=p(2,1); z1=p(3,1); e1=ep(1);
            x2=p(1,2); y2=p(2,2); z2=p(3,2); e2=ep(2);
        end
        function releaseImpl(obj)
            if ~isempty(obj.aurora)
                delete(obj.aurora);
                obj.aurora = [];
            end
        end
        function sts = getSampleTimeImpl(obj)
            sts = createSampleTime(obj,'Type','Discrete','SampleTime',obj.SamplePeriod);
        end
        function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13] = getOutputSizeImpl(~)
            o1=[1 1]; o2=[1 1]; o3=[1 1]; o4=[1 1]; o5=[1 1]; o6=[1 1]; o7=[1 1]; o8=[1 1]; o9=[1 2]; o10=[1 2]; o11=[1 1]; o12=[1 1]; o13=[1 1];
        end
        function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13] = getOutputDataTypeImpl(~)
            o1='double'; o2='double'; o3='double'; o4='double'; o5='double'; o6='double'; o7='double'; o8='double'; o9='logical'; o10='double'; o11='double'; o12='logical'; o13='double';
        end
        function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13] = isOutputComplexImpl(~)
            o1=false; o2=false; o3=false; o4=false; o5=false; o6=false; o7=false; o8=false; o9=false; o10=false; o11=false; o12=false; o13=false;
        end
        function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13] = isOutputFixedSizeImpl(~)
            o1=true; o2=true; o3=true; o4=true; o5=true; o6=true; o7=true; o8=true; o9=true; o10=true; o11=true; o12=true; o13=true;
        end
    end
end
