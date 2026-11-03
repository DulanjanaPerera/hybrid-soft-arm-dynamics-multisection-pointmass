classdef NdiReplayState < handle
    % State fixture for executing the exact saved/current reader step bodies.
    properties
        aurora
        clockStart
        nextPoll=0
        SamplePeriod=1
        BaseSensorId='0A'
        TipSensorId='0B'
        calibration
        previousFrame=[NaN NaN]
        thetaPrior=0
        lastThetaOrientation=NaN
        lastPhiOrientation=NaN
        lastLengthChange_m=nan(3,1)
        lastNonBendingRotation_rad=NaN
        lastR_NDI_from_arm=nan(3)
        lastGravity_arm=nan(3,1)
        ReferenceGravity_arm=[0;0;-9.81]
        ActuatorRadius=.013
    end
    methods
        function obj=NdiReplayState(driver,cal)
            obj.aurora=driver;obj.calibration=cal;obj.clockStart=tic;
        end
    end
end
