        function [x1,y1,z1,e1,x2,y2,z2,e2,valid,frame,hostTime,overrun,pollStart, ...
                thetaOrientation,phiOrientation,lengthChange_m,orientationValid, ...
                thetaObservable,nonBendingRotation_rad,R_NDI_from_arm, ...
                gravity_arm,baseOrientationValid,ndiDiagnostics] = updatedStepBaseRotation(obj)
            % Pace physical polls, not just simulated time. Never burst to catch up.
            if isempty(obj.clockStart), obj.clockStart=tic; end
            remaining = obj.nextPoll-toc(obj.clockStart);
            if remaining > 0, pause(remaining); end
            pollStart = toc(obj.clockStart);
            [x,y,z,ep,valid,frame,quaternion,rawDiagnostics] = obj.aurora.measureTipPositionAll( ...
                {obj.BaseSensorId,obj.TipSensorId});
            hostTime = toc(obj.clockStart);
            previousFrame=obj.previousFrame;
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
            ndiDiagnostics=ndiFinalizeDiagnostics(rawDiagnostics,frame,previousFrame, ...
                valid,~isempty(obj.calibration),false,pollStart,hostTime);
            R_NDI_from_arm = obj.lastR_NDI_from_arm;
            gravity_arm = obj.lastGravity_arm;
            baseOrientationValid = false;
            if valid(1) && ~isempty(obj.calibration)
                [newR,newGravity,baseOrientationValid] = ...
                    AuroraReader_s2_BaseRotation.armGravityFromBaseQuaternion( ...
                    quaternion(:,1),obj.calibration,obj.ReferenceGravity_arm);
                if baseOrientationValid
                    R_NDI_from_arm = newR;
                    gravity_arm = newGravity;
                    obj.lastR_NDI_from_arm = newR;
                    obj.lastGravity_arm = newGravity;
                end
            end
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
            ndiDiagnostics=ndiFinalizeDiagnostics(rawDiagnostics,frame,previousFrame, ...
                valid,true,orientationValid,pollStart,hostTime);
        end
