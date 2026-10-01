classdef AuroraReader_s2_Position_only < matlab.System
    % Outputs 1:13 retain the original XYZ and acquisition diagnostics.
    % Outputs 14:19 add calibrated orientation theta/phi, inferred lengths,
    % orientation validity, theta observability, and non-bending rotation.
    % One BX per step supplies both XYZ and quaternions.
    properties (Nontunable)
        Port = 'COM12'
        BaudRate = 921600
        SamplePeriod = 0.05
        BaseSensorId = '0A'
        TipSensorId = '0B'
        CalibrationFile = ''
        PreflightSeconds = 5
        ActuatorRadius = 0.013
    end
    properties(Access = private)
        aurora
        clockStart
        nextPoll = 0
        previousFrame = [NaN NaN]
        calibration
        thetaPrior = 0
        lastThetaOrientation = NaN
        lastPhiOrientation = NaN
        lastLengthChange_m = nan(3,1)
        lastNonBendingRotation_rad = NaN
    end
    methods
        function obj = AuroraReader_s2_Position_only(varargin)
            cfg = ndiSupportFnc.acquisitionSettings();
            obj.Port=cfg.Port; obj.BaudRate=cfg.BaudRate;
            obj.SamplePeriod=cfg.SamplePeriod;
            obj.BaseSensorId=cfg.BaseSensorId; obj.TipSensorId=cfg.TipSensorId;
            setProperties(obj,nargin,varargin{:});
        end
    end
    methods(Access = protected)
        function validatePropertiesImpl(obj)
            validateattributes(obj.SamplePeriod,{'double'},{'scalar','finite','positive'});
            validateattributes(obj.ActuatorRadius,{'double'}, ...
                {'scalar','finite','positive'});
            validateattributes(obj.PreflightSeconds,{'double'}, ...
                {'scalar','finite','positive'});
        end
        function setupImpl(obj)
            % Validate an explicit calibration before acquiring the port.
            obj.calibration = ndiSupportFnc.loadCalibration( ...
                obj.CalibrationFile, {obj.BaseSensorId,obj.TipSensorId});
            obj.aurora = AuroraDriver_2(obj.Port);
            try
                obj.aurora.openSerialPort();
                obj.aurora.init(obj.BaudRate);
                obj.aurora.detectAndAssignPortHandles();
                obj.aurora.sensorIndices({obj.BaseSensorId,obj.TipSensorId});
                obj.aurora.initPortHandleAll();
                obj.aurora.enablePortHandleDynamicAll();
                obj.aurora.startTracking();
                pause(0.1);
                if ~isempty(obj.calibration)
                    fprintf('NDI: hold the arm in its straight reference pose for preflight.\n');
                    cfg = ndiSupportFnc.acquisitionSettings();
                    cfg.Port = obj.Port;
                    cfg.BaudRate = obj.BaudRate;
                    cfg.SamplePeriod = obj.SamplePeriod;
                    cfg.BaseSensorId = obj.BaseSensorId;
                    cfg.TipSensorId = obj.TipSensorId;
                    report = ndiSupportFnc.preflight( ...
                        obj.calibration,obj.PreflightSeconds,obj.aurora,cfg);
                    assert(report.passed,'NDI:PreflightFailed', ...
                        'NDI preflight failed: %s',report.reason);
                end
                obj.clockStart = [];
                obj.nextPoll = 0;
                obj.previousFrame = [NaN NaN];
                obj.thetaPrior = 0;
                obj.lastThetaOrientation = NaN;
                obj.lastPhiOrientation = NaN;
                obj.lastLengthChange_m = nan(3,1);
                obj.lastNonBendingRotation_rad = NaN;
            catch e
                releaseImpl(obj);
                rethrow(e);
            end
        end
        function [x1,y1,z1,e1,x2,y2,z2,e2,valid,frame,hostTime,overrun,pollStart, ...
                thetaOrientation,phiOrientation,lengthChange_m,orientationValid, ...
                thetaObservable,nonBendingRotation_rad] = stepImpl(obj)
            % Pace physical polls, not just simulated time. Never burst to catch up.
            if isempty(obj.clockStart), obj.clockStart=tic; end
            remaining = obj.nextPoll-toc(obj.clockStart);
            if remaining > 0, pause(remaining); end
            pollStart = toc(obj.clockStart);
            [x,y,z,ep,valid,frame,quaternion] = obj.aurora.measureTipPositionAll( ...
                {obj.BaseSensorId,obj.TipSensorId});
            hostTime = toc(obj.clockStart);
            [p,valid] = ndiSupportFnc.validSample( ...
                [x;y;z],valid,frame,obj.previousFrame);
            obj.previousFrame = frame;
            overrun = hostTime > obj.nextPoll+obj.SamplePeriod;
            obj.nextPoll = max(obj.nextPoll+obj.SamplePeriod,pollStart+obj.SamplePeriod);
            x1=p(1,1); y1=p(2,1); z1=p(3,1); e1=ep(1);
            x2=p(1,2); y2=p(2,2); z2=p(3,2); e2=ep(2);
            % Plot/control outputs hold their last estimate through a lost
            % frame. The validity outputs still identify it as stale data.
            thetaOrientation = obj.lastThetaOrientation;
            phiOrientation = obj.lastPhiOrientation;
            lengthChange_m = obj.lastLengthChange_m;
            orientationValid = false; thetaObservable = false;
            nonBendingRotation_rad = obj.lastNonBendingRotation_rad;
            if ~all(valid) || isempty(obj.calibration), return; end
            est = ndiSupportFnc.orientationFromPose( ...
                quaternion,obj.calibration,obj.thetaPrior);
            if ~est.valid, return; end
            orientationValid = true;
            thetaObservable = est.thetaObservable;
            if thetaObservable
                thetaOrientation = est.theta;
                obj.lastThetaOrientation = est.theta;
            end
            phiOrientation = est.phi;
            obj.lastPhiOrientation = est.phi;
            lengthChange_m = ndiSupportFnc.anglesToLengths( ...
                est.theta,est.phi,obj.ActuatorRadius);
            obj.lastLengthChange_m = lengthChange_m;
            nonBendingRotation_rad = est.nonBendingRotation_rad;
            obj.lastNonBendingRotation_rad = nonBendingRotation_rad;
            obj.thetaPrior = est.theta;
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
        function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13, ...
                o14,o15,o16,o17,o18,o19] = getOutputSizeImpl(~)
            o1=[1 1]; o2=[1 1]; o3=[1 1]; o4=[1 1];
            o5=[1 1]; o6=[1 1]; o7=[1 1]; o8=[1 1];
            o9=[1 2]; o10=[1 2]; o11=[1 1]; o12=[1 1]; o13=[1 1];
            o14=[1 1]; o15=[1 1]; o16=[3 1];
            o17=[1 1]; o18=[1 1]; o19=[1 1];
        end
        function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13, ...
                o14,o15,o16,o17,o18,o19] = getOutputDataTypeImpl(~)
            o1='double'; o2='double'; o3='double'; o4='double';
            o5='double'; o6='double'; o7='double'; o8='double';
            o9='logical'; o10='double'; o11='double'; o12='logical'; o13='double';
            o14='double'; o15='double'; o16='double';
            o17='logical'; o18='logical'; o19='double';
        end
        function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13, ...
                o14,o15,o16,o17,o18,o19] = isOutputComplexImpl(~)
            o1=false; o2=false; o3=false; o4=false;
            o5=false; o6=false; o7=false; o8=false;
            o9=false; o10=false; o11=false; o12=false; o13=false;
            o14=false; o15=false; o16=false;
            o17=false; o18=false; o19=false;
        end
        function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13, ...
                o14,o15,o16,o17,o18,o19] = isOutputFixedSizeImpl(~)
            o1=true; o2=true; o3=true; o4=true;
            o5=true; o6=true; o7=true; o8=true;
            o9=true; o10=true; o11=true; o12=true; o13=true;
            o14=true; o15=true; o16=true;
            o17=true; o18=true; o19=true;
        end
    end
end
