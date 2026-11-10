function ports = yost_list_ports
%YOST_LIST_PORTS List available serial ports without opening a device.
ports = serialportlist('available');
disp(ports);
fprintf('Use dongle ports (for example COM28), not the sensor charging port.\n');
end
