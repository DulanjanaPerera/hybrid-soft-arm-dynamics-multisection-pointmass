% Read XYZ from two Aurora sensors using one BX poll per sample.
% Edit the port and sample count for the current setup. Sensor order follows
% discovered port-handle order; identify base/tip by sensor IDs before use.

serialPort = 'COM12';
baudRate = 921600;
sampleCount = 500;

aurora_device = AuroraDriver_2(serialPort);
trackingStarted = false;
ndiSamples = struct();

try
    aurora_device.openSerialPort();
    aurora_device.init();
    aurora_device.BEEP('2');
    aurora_device.setBaudRate(baudRate);
    aurora_device.detectAndAssignPortHandles();
    if aurora_device.n_port_handles < 2
        error('NDI:TwoSensorsRequired', ...
            'Expected two sensors, found %d.', aurora_device.n_port_handles);
    end
    aurora_device.initPortHandleAll();
    aurora_device.enablePortHandleDynamicAll();

    ndiSamples.sensorIds = {aurora_device.port_handles(1,1).id, ...
        aurora_device.port_handles(1,2).id};
    ndiSamples.hostTime_s = nan(sampleCount,1);
    ndiSamples.position_mm = nan(3,2,sampleCount);
    ndiSamples.indicator = nan(sampleCount,2);
    ndiSamples.frame = nan(sampleCount,2);
    ndiSamples.valid = false(sampleCount,2);

    aurora_device.startTracking();
    trackingStarted = true;
    disp('Start in 2 seconds');
    pause(2);
    aurora_device.BEEP('1');
    disp('Start');

    clockStart = tic;
    for i = 1:sampleCount
        [x,y,z,indicator,valid,frame] = aurora_device.measureTipPositionAll();
        ndiSamples.hostTime_s(i) = toc(clockStart);
        ndiSamples.position_mm(:,:,i) = [x;y;z];
        ndiSamples.indicator(i,:) = indicator;
        ndiSamples.frame(i,:) = frame;
        ndiSamples.valid(i,:) = valid;
        if mod(i,50) == 0 || i == sampleCount
            fprintf('%d/%d samples, valid sensors: [%d %d]\n', ...
                i,sampleCount,valid(1),valid(2));
        end
    end
catch acquisitionError
    cleanupAurora(aurora_device,trackingStarted);
    rethrow(acquisitionError);
end

cleanupAurora(aurora_device,trackingStarted);
ndiSamples.position_m = ndiSamples.position_mm * 1e-3;
clear aurora_device trackingStarted clockStart
disp('End. Data remains in ndiSamples (position_mm and position_m).');

function cleanupAurora(device,trackingStarted)
    if trackingStarted
        try
            device.stopTracking();
        catch stopError
            warning('NDI:StopTrackingFailed', '%s', stopError.message);
        end
    end
    try
        delete(device);
    catch deleteError
        warning('NDI:DeleteFailed', '%s', deleteError.message);
    end
end
