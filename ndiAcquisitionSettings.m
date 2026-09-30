function cfg = ndiAcquisitionSettings()
% Shared acquisition defaults for MATLAB and Simulink.
cfg.Port = 'COM12';
cfg.BaudRate = 921600;
cfg.SamplePeriod = 0.05; % 20 Hz host polling; measured time remains authoritative.
cfg.BaseSensorId = '0A';
cfg.TipSensorId = '0B';
end
