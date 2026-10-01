classdef ndiSupportFnc
    % Orientation and constant-curvature geometry shared by the NDI reader.
    % NDI quaternions are scalar-first and map sensor-local axes to NDI axes.
    methods (Static)
        function L = sensorGeometry()
            % [Base offset; module length; tip offset], metres.
            L = [0.055;0.17411;0.05];
        end

        function cfg = acquisitionSettings()
            % Shared MATLAB and Simulink NDI acquisition defaults.
            cfg.Port = 'COM12';
            cfg.BaudRate = 921600;
            cfg.SamplePeriod = 0.05;
            cfg.BaseSensorId = '0A';
            cfg.TipSensorId = '0B';
        end

        function [theta,phi,residual,valid] = configuration(p,phiPrior,thetaPrior)
            % XYZ comparison only; never pass a missing reading into IK.
            theta = NaN; phi = NaN; residual = NaN;
            valid = all(isfinite(p));
            if ~valid, return; end
            if ~isfinite(phiPrior), phiPrior = 0.0001; end
            if ~isfinite(thetaPrior), thetaPrior = 0; end
            [theta,phi,residual] = f20260219_2_task2config_withLenExt( ...
                p,ndiSupportFnc.sensorGeometry(),1,phiPrior,thetaPrior);
            valid = isreal(theta) && isreal(phi) && isreal(residual) ...
                && all(isfinite([theta,phi,residual]));
        end

        function [position,valid] = validSample(position,valid,frame,previousFrame)
            % Reject missing, nonfinite, repeated, or mismatched frames.
            valid = logical(valid) & all(isfinite(position),1) ...
                & isfinite(frame) & (frame ~= previousFrame);
            if frame(1) ~= frame(2), valid(:) = false; end
            position(:,~valid) = NaN;
        end

        function [cal, file, data] = runCalibration(recordingDir)
            % Run manually with the arm straight, untwisted, and still.
            % One MAT-file retains both the raw evidence and the fitted cal.
            if nargin < 1 || isempty(recordingDir)
                recordingDir = fullfile(fileparts(mfilename('fullpath')), ...
                    'ndi_orientation_recordings');
            end
            reply = input(['Hold the arm STRAIGHT and UNTWISTED with both sensors ' ...
                'still. Small fixed X/Y mounting offsets are acceptable. ' ...
                'Type STRAIGHT to capture five seconds: '], 's');
            assert(strcmp(strtrim(reply), 'STRAIGHT'), ...
                'NDI:CalibrationCancelled', 'Calibration cancelled.');
            if ~isfolder(recordingDir), mkdir(recordingDir); end
            stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
            file = fullfile(recordingDir, ['calibration_' stamp '.mat']);
            data = ndiSupportFnc.capturePose(5);
            save(file, 'data'); % Preserve partial capture if fitting fails.
            cal = ndiSupportFnc.fitCalibration(data, true);
            save(file, 'cal', 'data');
            fprintf('Saved calibration and raw capture: %s\n', file);
            fprintf('Accepted %d/%d paired poses; XYZ reference error %.2f mm.\n', ...
                cal.acceptedCount, size(data.quaternion,3), ...
                1000*cal.straightPositionError_m);
            fprintf(['Maximum motion: orientation [%.2f %.2f] deg, ' ...
                'position [%.2f %.2f] mm; overruns %d.\n'], ...
                cal.maxOrientationDeviation_deg, ...
                cal.maxPositionDeviation_mm,sum(data.overrun));
            fprintf('Verify this calibration with independent bends before control use.\n');
        end

        function data = capturePose(seconds, driver, cfg)
            % One BX poll per sample; owned ports are always closed on exit.
            if nargin < 1 || isempty(seconds), seconds = 5; end
            validateattributes(seconds, {'numeric'}, ...
                {'scalar','finite','positive'});
            if nargin < 3 || isempty(cfg)
                cfg = ndiSupportFnc.acquisitionSettings();
            end
            count = max(1, round(seconds/cfg.SamplePeriod));
            owned = nargin < 2 || isempty(driver);
            if owned
                driver = AuroraDriver_2(cfg.Port);
                cleanup = onCleanup(@() delete(driver));
                driver.openSerialPort();
                driver.init(cfg.BaudRate);
                driver.detectAndAssignPortHandles();
                driver.sensorIndices({cfg.BaseSensorId,cfg.TipSensorId});
                driver.initPortHandleAll();
                driver.enablePortHandleDynamicAll();
                driver.startTracking();
                pause(0.1);
            end
            data.sensorIds = {cfg.BaseSensorId,cfg.TipSensorId};
            data.settings = cfg;
            data.geometry_m = ndiSupportFnc.sensorGeometry();
            data.createdUtc = char(datetime('now','TimeZone','UTC'));
            data.quaternionConvention = 'q0 qx qy qz; sensor-local to NDI';
            data.position_mm = nan(3,2,count);
            data.quaternion = nan(4,2,count);
            data.valid = false(count,2);
            data.orientationValid = false(count,2);
            data.frame = nan(count,2);
            data.indicator = nan(count,2);
            data.pollStart_s = nan(count,1);
            data.hostTime_s = nan(count,1);
            data.overrun = false(count,1);
            data.completedSamples = 0;
            data.error = '';
            clockStart = tic; nextPoll = 0; previousFrame = [NaN NaN];
            try
                for k = 1:count
                    remaining = nextPoll-toc(clockStart);
                    if remaining > 0, pause(remaining); end
                    started = toc(clockStart);
                    [x,y,z,e,v,f,q] = driver.measureTipPositionAll(data.sensorIds);
                    finished = toc(clockStart);
                    [p,v] = ndiSupportFnc.validSample( ...
                        [x;y;z],v,f,previousFrame);
                    previousFrame = f;
                    q(:,~v) = NaN;
                    for j = 1:2
                        [~,ok] = ndiSupportFnc.quaternionRotation(q(:,j));
                        data.orientationValid(k,j) = v(j) && ok;
                        if ~ok, q(:,j) = NaN; end
                    end
                    data.position_mm(:,:,k) = p;
                    data.quaternion(:,:,k) = q;
                    data.valid(k,:) = v;
                    data.frame(k,:) = f;
                    data.indicator(k,:) = e;
                    data.pollStart_s(k) = started;
                    data.hostTime_s(k) = finished;
                    data.overrun(k) = finished > nextPoll+cfg.SamplePeriod;
                    nextPoll = max(nextPoll+cfg.SamplePeriod, ...
                        started+cfg.SamplePeriod);
                    data.completedSamples = k;
                end
            catch e
                data.error = getReport(e,'extended','hyperlinks','off');
                warning('NDI:PoseCaptureStopped', ...
                    'Capture stopped after %d/%d samples: %s', ...
                    data.completedSamples,count,e.message);
            end
            data.position_m = data.position_mm*1e-3;
        end

        function cal = fitCalibration(data, straightPoseConfirmed)
            % A known straight pose establishes mounting rotations. A single
            % unknown bent pose cannot separate mounting error from bend.
            assert(nargin >= 2 && isequal(straightPoseConfirmed,true), ...
                'NDI:StraightPoseRequired', ...
                'Confirm a known straight, untwisted pose explicitly.');
            cfg = ndiSupportFnc.acquisitionSettings();
            assert(isequal(data.sensorIds,{cfg.BaseSensorId,cfg.TipSensorId}), ...
                'NDI:SensorOrder','Capture must use the confirmed base/tip IDs.');
            assert(isequal(data.geometry_m,ndiSupportFnc.sensorGeometry()), ...
                'NDI:CalibrationGeometry', ...
                'Capture geometry differs from the current arm geometry.');
            assert(isempty(data.error) && ...
                data.completedSamples == size(data.quaternion,3), ...
                'NDI:IncompleteCalibration', ...
                'Capture failed or ended early; repeat calibration.');
            % Arm X=NDI X, arm Y=NDI Z, arm Z=-NDI Y.
            R_NDI_from_arm = [1 0 0;0 0 -1;0 1 0];
            [accepted, rotations] = ndiSupportFnc.acceptedPoses(data);
            assert(sum(accepted) >= 20 && mean(accepted) >= 0.9, ...
                'NDI:CalibrationQuality', ...
                'Need at least 20 paired poses and 90 percent acceptance.');
            meanR = nan(3,3,2); maxAngle = zeros(1,2);
            maxPosition = zeros(1,2);
            for j = 1:2
                M = mean(rotations(:,:,j,accepted),4);
                [U,~,V] = svd(M);
                meanR(:,:,j) = U*diag([1 1 det(U*V')])*V';
                for k = find(accepted)'
                    d = meanR(:,:,j)'*rotations(:,:,j,k);
                    maxAngle(j) = max(maxAngle(j), ...
                        acos(max(-1,min(1,(trace(d)-1)/2))));
                end
                p = reshape(data.position_mm(:,j,accepted),3,[]);
                maxPosition(j) = max(vecnorm(p-mean(p,2)));
            end
            assert(all(maxAngle < 2*pi/180) && all(maxPosition < 2), ...
                'NDI:CalibrationMoved', ...
                'Pose moved by over 2 deg or 2 mm during capture.');
            cal.schemaVersion = 1;
            cal.createdUtc = char(datetime('now','TimeZone','UTC'));
            cal.sensorIds = data.sensorIds;
            cal.geometry_m = ndiSupportFnc.sensorGeometry();
            cal.R_NDI_from_arm_reference = R_NDI_from_arm;
            cal.R_baseSensor_from_arm = meanR(:,:,1)'*R_NDI_from_arm;
            cal.R_tipSensor_from_tip = meanR(:,:,2)'*R_NDI_from_arm;
            cal.meanSensorRotations = meanR;
            cal.acceptedMask = accepted;
            cal.acceptedCount = sum(accepted);
            cal.acceptedFraction = mean(accepted);
            cal.maxOrientationDeviation_deg = maxAngle*180/pi;
            cal.maxPositionDeviation_mm = maxPosition;
            cal.straightAngleTol_rad = 0.5*pi/180;
            cal.quaternionConvention = data.quaternionConvention;
            cal.straightPoseConfirmed = true;
            delta = mean(data.position_mm(:,2,accepted) ...
                -data.position_mm(:,1,accepted),3)*1e-3;
            cal.straightPosition_arm_m = R_NDI_from_arm'*delta;
            cal.straightPositionError_m = norm( ...
                cal.straightPosition_arm_m-[0;0;sum(cal.geometry_m)]);
            % Fixed mounting offsets of roughly +/-5 mm in each lateral axis
            % are retained as a measured reference, never zeroed away.
            if cal.straightPositionError_m > 0.015
                warning('NDI:StraightGeometryMismatch', ...
                    'Reference XYZ differs from expected by %.2f mm; inspect alignment.', ...
                    1000*cal.straightPositionError_m);
            end
        end

        function [report, data] = preflight(calibrationFile, seconds, driver, cfg)
            % Check a held straight pose before a run; never recalibrate it.
            if nargin < 2 || isempty(seconds), seconds = 5; end
            if nargin < 4 || isempty(cfg)
                cfg = ndiSupportFnc.acquisitionSettings();
            end
            if isstruct(calibrationFile)
                cal = calibrationFile; % Already validated by the reader.
            else
                cal = ndiSupportFnc.loadCalibration(calibrationFile, ...
                    {cfg.BaseSensorId,cfg.TipSensorId});
            end
            assert(~isempty(cal),'NDI:PreflightCalibration', ...
                'Preflight needs a verified calibration file.');
            if nargin >= 3
                data = ndiSupportFnc.capturePose(seconds,driver,cfg);
            else
                data = ndiSupportFnc.capturePose(seconds,[],cfg);
            end
            report = ndiSupportFnc.evaluatePreflight(data,cal);
            fprintf(['NDI preflight %s/%s: valid [%d %d], paired %d/%d; ' ...
                'frames [%g %g] to [%g %g]; %d overruns; ' ...
                'median interval %.4f s; bend %.2f deg; ' ...
                'extra rotation %.2f deg; XYZ change %.2f mm; %s.\n'], ...
                data.sensorIds{1},data.sensorIds{2}, ...
                report.validCounts,report.validPairs,report.samples, ...
                report.firstFrame,report.lastFrame,report.overruns, ...
                report.medianInterval_s,report.medianBend_deg, ...
                report.medianExtraRotation_deg, ...
                report.referencePositionChange_mm,report.reason);
        end

        function report = evaluatePreflight(data,cal)
            % Diagnostic gates for quality and reference consistency, not
            % manufacturer accuracy specifications.
            [accepted,~] = ndiSupportFnc.acceptedPoses(data);
            n = size(data.quaternion,3);
            report.samples = n;
            report.validPairs = sum(accepted);
            report.validCounts = sum(data.valid,1);
            report.overruns = sum(data.overrun);
            report.firstFrame = [NaN NaN];
            report.lastFrame = [NaN NaN];
            if data.completedSamples > 0
                report.firstFrame = data.frame(1,:);
                report.lastFrame = data.frame(data.completedSamples,:);
            end
            report.medianInterval_s = NaN;
            report.medianBend_deg = NaN;
            report.medianExtraRotation_deg = NaN;
            report.referencePositionChange_mm = NaN;
            report.p95BendDeviation_deg = NaN;
            report.p95PositionDeviation_mm = NaN;
            report.reason = '';
            report.passed = false;
            if data.completedSamples < n || ~isempty(data.error)
                report.reason = 'Capture ended early or raised an error'; return;
            end
            times = data.pollStart_s(isfinite(data.pollStart_s));
            if numel(times) > 1
                report.medianInterval_s = median(diff(times));
            end
            if n < 20 || report.validPairs < 20 || mean(accepted) < 0.9
                report.reason = 'Too few valid paired poses'; return;
            end
            if ~isfinite(report.medianInterval_s) || ...
                    report.medianInterval_s > 1.5*data.settings.SamplePeriod
                report.reason = 'Polling slower than the configured period'; return;
            end
            assert(isfield(cal,'straightPosition_arm_m') && ...
                isequal(size(cal.straightPosition_arm_m),[3 1]), ...
                'NDI:PreflightCalibration', ...
                'Calibration has no straight-pose XYZ reference.');
            bend = nan(n,1); extraRotation = nan(n,1);
            relativePosition = nan(3,n);
            for k = find(accepted)'
                est = ndiSupportFnc.orientationFromPose( ...
                    data.quaternion(:,:,k),cal,0);
                [B,ok] = ndiSupportFnc.quaternionRotation( ...
                    data.quaternion(:,1,k));
                if ~est.valid || ~ok, continue; end
                bend(k) = est.phi*180/pi;
                extraRotation(k) = est.nonBendingRotation_rad*180/pi;
                A = B*cal.R_baseSensor_from_arm;
                relativePosition(:,k) = A'*( ...
                    data.position_m(:,2,k)-data.position_m(:,1,k));
            end
            usable = accepted & isfinite(bend);
            if sum(usable) < 20
                report.reason = 'Too few calibrated orientation poses'; return;
            end
            report.medianBend_deg = median(bend(usable));
            report.medianExtraRotation_deg = median(extraRotation(usable));
            position = median(relativePosition(:,usable),2);
            report.referencePositionChange_mm = 1000*norm( ...
                position-cal.straightPosition_arm_m);
            bendDeviation = sort(abs(bend(usable)-report.medianBend_deg));
            positionDeviation = sort(1000*vecnorm( ...
                relativePosition(:,usable)-position));
            p95 = ceil(0.95*sum(usable));
            report.p95BendDeviation_deg = bendDeviation(p95);
            report.p95PositionDeviation_mm = positionDeviation(p95);
            if report.medianBend_deg > 5
                report.reason = 'Tip orientation differs from straight reference';
            elseif report.medianExtraRotation_deg > 5
                report.reason = 'Relative sensor twist differs from reference';
            elseif report.referencePositionChange_mm > 15
                report.reason = 'Relative XYZ differs from calibration reference';
            elseif report.p95BendDeviation_deg > 2 || ...
                    report.p95PositionDeviation_mm > 3
                report.reason = 'Reference pose was not held still';
            else
                report.passed = true;
                report.reason = 'Reference pose and acquisition quality passed';
            end
        end

        function [accepted,rotations] = acceptedPoses(data)
            % Preserve invalid/missing/repeated-frame readings as invalid.
            n = size(data.quaternion,3);
            accepted = false(n,1);
            rotations = nan(3,3,2,n);
            previousFrame = [NaN NaN];
            for k = 1:n
                ok = false(1,2);
                for j = 1:2
                    [rotations(:,:,j,k),ok(j)] = ...
                        ndiSupportFnc.quaternionRotation( ...
                        data.quaternion(:,j,k));
                end
                frame = data.frame(k,:);
                accepted(k) = all(data.valid(k,:)) && ...
                    all(data.orientationValid(k,:)) && all(ok) && ...
                    all(isfinite(data.position_mm(:,:,k)),'all') && ...
                    all(isfinite(frame)) && frame(1) == frame(2) && ...
                    ~any(frame == previousFrame);
                previousFrame = frame;
            end
        end

        function cal = loadCalibration(file, sensorIds)
            % An empty path leaves orientation disabled while XYZ stays usable.
            cal = [];
            if isempty(file) || (isstring(file) && isscalar(file) ...
                    && strlength(file) == 0), return; end
            assert(ischar(file) || (isstring(file) && isscalar(file)), ...
                'NDI:CalibrationFile', 'CalibrationFile must be a file path.');
            file = char(file);
            % A path copied from the MATLAB command window can include
            % display quotes. They are not part of a Windows filename.
            if numel(file) >= 2 && ...
                    ((file(1) == '''' && file(end) == '''') || ...
                     (file(1) == '"' && file(end) == '"'))
                file = file(2:end-1);
            end
            assert(isfile(file), 'NDI:CalibrationFile', ...
                'Orientation calibration file not found: %s', file);
            saved = load(file, 'cal');
            assert(isfield(saved, 'cal') && isstruct(saved.cal), ...
                'NDI:CalibrationFile', 'The file must contain a cal struct.');
            cal = saved.cal;
            fields = {'sensorIds','geometry_m','straightPoseConfirmed', ...
                'R_baseSensor_from_arm','R_tipSensor_from_tip', ...
                'straightAngleTol_rad'};
            assert(all(isfield(cal, fields)), 'NDI:CalibrationFile', ...
                'Orientation calibration is missing required fields.');
            assert(isequal(cal.sensorIds, sensorIds), 'NDI:SensorOrder', ...
                'Calibration sensor IDs must match the reader base/tip IDs.');
            assert(isequal(cal.geometry_m, ndiSupportFnc.sensorGeometry()), ...
                'NDI:CalibrationGeometry', ...
                'Calibration offsets or arm length do not match this project.');
            assert(isequal(cal.straightPoseConfirmed, true), ...
                'NDI:CalibrationPose', ...
                'Calibration must come from a confirmed straight pose.');
            assert(isscalar(cal.straightAngleTol_rad) ...
                && isfinite(cal.straightAngleTol_rad) ...
                && cal.straightAngleTol_rad > 0 ...
                && cal.straightAngleTol_rad < pi/2, 'NDI:CalibrationFile', ...
                'Invalid straight-angle observability threshold.');
            rotations = {cal.R_baseSensor_from_arm, cal.R_tipSensor_from_tip};
            for k = 1:2
                R = rotations{k};
                assert(isequal(size(R), [3 3]) && all(isfinite(R(:))) ...
                    && norm(R.'*R-eye(3), 'fro') < 1e-8 && det(R) > 0, ...
                    'NDI:CalibrationFile', ...
                    'Calibration contains an invalid sensor mounting rotation.');
            end
        end

        function est = orientationFromPose(quaternion, cal, thetaPrior)
            % Base-relative tip tangent, expressed in the arm frame.
            if nargin < 3 || ~isfinite(thetaPrior), thetaPrior = 0; end
            est = struct('valid', false, 'theta', NaN, 'phi', NaN, ...
                'thetaObservable', false, 'nonBendingRotation_rad', NaN);
            if isempty(cal) || ~isequal(size(quaternion), [4 2]), return; end
            [B, okB] = ndiSupportFnc.quaternionRotation(quaternion(:,1));
            [T, okT] = ndiSupportFnc.quaternionRotation(quaternion(:,2));
            if ~okB || ~okT, return; end
            A = B*cal.R_baseSensor_from_arm;
            R = A.'*(T*cal.R_tipSensor_from_tip);
            tangent = R(:,3);
            phi = atan2(hypot(tangent(1), tangent(2)), tangent(3));
            observable = hypot(tangent(1), tangent(2)) ...
                > sin(cal.straightAngleTol_rad);
            if observable
                theta = atan2(tangent(2), tangent(1));
            else
                theta = atan2(sin(thetaPrior), cos(thetaPrior));
            end
            c = cos(theta); s = sin(theta);
            cp = cos(phi); sp = sin(phi);
            Rcc = [c*c*cp+s*s, c*s*(cp-1), c*sp; ...
                   c*s*(cp-1), s*s*cp+c*c, s*sp; ...
                   -c*sp, -s*sp, cp];
            extra = Rcc.'*R;
            est.valid = true;
            est.theta = theta;
            est.phi = phi;
            est.thetaObservable = observable;
            est.nonBendingRotation_rad = acos(max(-1, ...
                min(1, (trace(extra)-1)/2)));
        end

        function [R, valid] = quaternionRotation(q)
            % NDI API: [q0; qx; qy; qz], sensor-local to NDI rotation.
            R = nan(3); valid = false;
            if numel(q) ~= 4 || any(~isfinite(q(:))), return; end
            q = double(q(:)); n = norm(q);
            if abs(n-1) > 0.02, return; end
            q = q/n; w = q(1); x = q(2); y = q(3); z = q(4);
            R = [w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y); ...
                 2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x); ...
                 2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z];
            valid = true;
        end

        function lengths = anglesToLengths(theta, phi, r)
            % Constant-curvature, zero common-extension estimate (metres).
            lengths = r*phi*[-cos(theta); ...
                cos(theta)/2-sqrt(3)*sin(theta)/2; ...
                cos(theta)/2+sqrt(3)*sin(theta)/2];
        end
    end
end
