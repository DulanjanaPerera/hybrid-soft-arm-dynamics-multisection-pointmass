classdef AuroraDriver_2 < handle
    properties (Constant)
        COMMAND_FORMAT_1 = 1;
        COMMAND_FORMAT_2 = 2;

        READ_OUT_OF_VOLUME_NOT_ALLOWED = '0001';
        READ_OUT_OF_VOLUME_ALLOWED     = '0801';

        SENSOR_STATUS_VALID    = '01';
        SENSOR_STATUS_MISSING  = '02';
        SENSOR_STATUS_DISABLED = '04';

        BAUD_9600   = '0';
        BAUD_14400  = '1';
        BAUD_19200  = '2';
        BAUD_38400  = '3';
        BAUD_57600  = '4';
        BAUD_115200 = '5';
        BAUD_921600 = '6';
        BAUD_230400 = 'A';

        TT_PRIORITY_STATIC  = 'S';
        TT_PRIORITY_DYNAMIC = 'D';
        TT_PRIORITY_BUTTON  = 'B';

        PHSR_HANDLES_ALL                      = '00';
        PHSR_HANDLES_TO_BE_FREED              = '01';
        PHSR_HANDLES_OCCUPIED                 = '02';
        PHSR_HANDLES_OCCUPIED_AND_INITIALIZED = '03';
        PHSR_HANDLES_ENABLED                  = '04';

        RESET_SOFT = '0';
        RESET_HARD = '1';

        TRACKING_OPTION_NONE                    = '';
        TRACKING_OPTION_FAST_MODE               = '40';
        TRACKING_OPTION_RESET_COUNTER           = '80';
        TRACKING_OPTION_FAST_MODE_RESET_COUNTER = 'C0';
    end

    properties (GetAccess = public, SetAccess = private)
        serial_port;
        n_port_handles;
        port_handles;

        selected_command_format;
        device_init;
    end

    methods (Access = public)
        function obj = AuroraDriver_2(serial_port)
            obj.serial_port = serial(serial_port);     % legacy serial object
            obj.serial_port.Terminator = 'CR';
            set(obj.serial_port,'Timeout',2);          % ADDED
            set(obj.serial_port,'FlowControl','none'); % ADDED
            obj.n_port_handles = 0;
            obj.selected_command_format = obj.COMMAND_FORMAT_2;
            obj.device_init = 0;
        end

        function openSerialPort(obj)
            if strcmp(obj.serial_port.Status, 'closed')
                fopen(obj.serial_port);
            end
        end

        function closeSerialPort(obj)
            if strcmp(obj.serial_port.Status, 'open')
                fclose(obj.serial_port);
            end
        end

        function setBaudRate(obj, baud_rate)
            % Decide whether to enable hardware RTS/CTS for high rates
            use_hw = baud_rate >= 115200;

            switch baud_rate
                case 9600,    baud_rate_code = obj.BAUD_9600;
                case 14400,   baud_rate_code = obj.BAUD_14400;
                case 19200,   baud_rate_code = obj.BAUD_19200;
                case 38400,   baud_rate_code = obj.BAUD_38400;
                case 57600,   baud_rate_code = obj.BAUD_57600;
                case 115200,  baud_rate_code = obj.BAUD_115200;
                case 921600,  baud_rate_code = obj.BAUD_921600;
                case 230400,  baud_rate_code = obj.BAUD_230400;
                otherwise
                    error('Invalid baud rate %d', baud_rate);
            end

            % Ask SCU to switch first (last digit = HW handshaking flag)
            hw_flag = obj.ternary(use_hw, '1', '0');                                     % ADDED
            reply = obj.COMM(baud_rate_code,'0','0','0',hw_flag);                        % CHANGED
            if isempty(reply) || ~strncmpi(reply,'OK',2)
                error('COMM failed: no OKAY reply from SCU (bad baud/port in use).');
            end

            pause(0.2); % >=100 ms after OKAY (NDI guidance)

            % Now host side: set flow control and baud to match
            if use_hw
                set(obj.serial_port,'FlowControl','hardware'); % ADDED
            else
                set(obj.serial_port,'FlowControl','none');     % ADDED
            end
            set(obj.serial_port,'BaudRate', baud_rate);
            pause(0.1);
        end

        function init(obj)
            obj.INIT();
            obj.device_init = 1;
        end

        function startTracking(obj)
            obj.TSTART(obj.TRACKING_OPTION_RESET_COUNTER);
        end

        function stopTracking(obj)
            obj.TSTOP();
        end

        function detectAndAssignPortHandles(obj)
            reply = obj.PHSR(obj.PHSR_HANDLES_ALL);
            if isempty(reply) || length(reply) < 2
                error('PHSR timeout/empty reply – check baud/handshaking.');
            end
            obj.n_port_handles = hex2dec(reply(1:2));
            for i_port_handle = 1:obj.n_port_handles
                s = 3 + 5*(i_port_handle - 1);
                id = reply(s:s+1);
                status = reply(s+2:s+4);
                if i_port_handle == 1
                    obj.port_handles = PortHandle(id, status);
                else
                    obj.port_handles(1,i_port_handle) = PortHandle(id, status);
                end
            end
        end

        function updatePortHandleStatusAll(obj)
            reply = obj.PHSR(obj.PHSR_HANDLES_ALL);
            if isempty(reply) || length(reply) < 2
                error('PHSR timeout/empty reply – check baud/handshaking.');
            end
            n_found_port_handles = hex2dec(reply(1:2));
            for i_found_port_handle = 1:n_found_port_handles
                s = 3 + 5*(i_found_port_handle - 1);
                id = reply(s:s+1);
                status = reply(s+2:s+4);
                for i_port_handle = 1:obj.n_port_handles
                    if strcmp(obj.port_handles(1,i_port_handle).id, id)
                        obj.port_handles(1,i_port_handle).updateStatus(status);
                        break;
                    end
                end
            end
        end

        function initPortHandle(obj, port_handle_id)
            obj.PINIT(port_handle_id);
        end

        function initPortHandleAll(obj)
            for i_port_handle = 1:obj.n_port_handles
                obj.initPortHandle(obj.port_handles(1,i_port_handle).id);
            end
            obj.updatePortHandleStatusAll();
        end

        function enablePortHandleDynamic(obj, port_handle_id)
            obj.PENA(port_handle_id, obj.TT_PRIORITY_DYNAMIC);
        end

        function enablePortHandleDynamicAll(obj)
            for i_port_handle = 1:obj.n_port_handles
                obj.enablePortHandleDynamic(obj.port_handles(1,i_port_handle).id);
            end
            obj.updatePortHandleStatusAll();
        end

        function updateSensorDataAll(obj)
            if obj.device_init == 1
                obj.sendCommand(sprintf('BX %s', obj.READ_OUT_OF_VOLUME_ALLOWED));

                % A handle omitted from this reply must not retain old data.
                for i_port_handle = 1:obj.n_port_handles
                    obj.port_handles(1,i_port_handle).updateSensorStatus( ...
                        obj.SENSOR_STATUS_MISSING);
                end

                start_sequence = obj.readScalar('uint16'); %#ok<NASGU>
                reply_length   = obj.readScalar('uint16'); %#ok<NASGU>
                header_CRC     = obj.readScalar('uint16'); %#ok<NASGU>

                num_handle_reads = obj.readScalar('uint8');

                for i_handle_reads = 1:num_handle_reads
                    handle_id     = dec2hex(obj.readScalar('uint8'), 2);
                    sensor_status = dec2hex(obj.readScalar('uint8'), 2);

                    handle_index = 0;
                    for i_port_handle = 1:obj.n_port_handles
                        if strcmp(obj.port_handles(1,i_port_handle).id, handle_id)
                            handle_index = i_port_handle;
                            break;
                        end
                    end
                    if handle_index > 0
                        obj.port_handles(1,handle_index).updateSensorStatus(sensor_status);
                    end

                    if strcmp(sensor_status, obj.SENSOR_STATUS_DISABLED) == 0
                        if strcmp(sensor_status, obj.SENSOR_STATUS_VALID)
                            q0 = obj.readScalar('float32');
                            qX = obj.readScalar('float32');
                            qY = obj.readScalar('float32');
                            qZ = obj.readScalar('float32');
                            rot = [q0 qX qY qZ];

                            tX = obj.readScalar('float32');
                            tY = obj.readScalar('float32');
                            tZ = obj.readScalar('float32');
                            trans = [tX tY tZ];

                            error_val = obj.readScalar('float32');

                            if handle_index > 0
                                obj.port_handles(1,handle_index).updateTrans(trans);
                                obj.port_handles(1,handle_index).updateRot(rot);
                                obj.port_handles(1,handle_index).updateError(error_val);
                            end
                        end

                        handle_status = dec2hex(obj.readScalar('uint32'), 8);
                        frame_number  = obj.readScalar('uint32');

                        if handle_index > 0
                            obj.port_handles(1,handle_index).updateStatusComplete(handle_status);
                            obj.port_handles(1,handle_index).updateFrameNumber(frame_number);
                        end
                    end
                end

                system_status = obj.readScalar('uint16'); %#ok<NASGU>
                crc           = obj.readScalar('uint16'); %#ok<NASGU>
            end
        end

        function status = readSensorStatus(obj)
            obj.updateSensorDataAll();
            status = obj.port_handles(1,1).sensor_status;
        end

        function [angle, error] = measureTipOrientation(obj)
            obj.updateSensorDataAll();
            rot = obj.port_handles(1,1).rot;
            [~, angle, ~] = quat2angle(rot);
            error = obj.port_handles(1,1).error;
        end

        function [rot1, rot2] = measureTipOrientationAll(obj)
            obj.updateSensorDataAll();
            rot1 = obj.port_handles(1,1).rot;
            rot2 = obj.port_handles(1,2).rot;
        end

        function [x, y, z, error_val, valid, frame] = measureTipPosition(obj)
            obj.updateSensorDataAll();
            if obj.n_port_handles < 1
                error('AuroraDriver_2:NoSensor', 'No sensor port handle was discovered.');
            end
            handle = obj.port_handles(1,1);
            valid = handle.hasUsablePosition();
            T = handle.trans;
            if ~valid, T(:) = NaN; end
            x = T(1); y = T(2); z = T(3);
            error_val = handle.error;
            frame = handle.frame_number;
        end

        function [x, y, z, error_val, valid, frame] = measureTipPositionAll(obj)
            obj.updateSensorDataAll();
            if obj.n_port_handles < 2
                error('AuroraDriver_2:TwoSensorsRequired', ...
                    'Two sensor port handles are required.');
            end
            handle1 = obj.port_handles(1,1);
            handle2 = obj.port_handles(1,2);
            valid = [handle1.hasUsablePosition(), handle2.hasUsablePosition()];
            T1 = handle1.trans;
            T2 = handle2.trans;
            if ~valid(1), T1(:) = NaN; end
            if ~valid(2), T2(:) = NaN; end
            x = [T1(1), T2(1)];
            y = [T1(2), T2(2)];
            z = [T1(3), T2(3)];
            error_val = [handle1.error, handle2.error];
            frame = [handle1.frame_number, handle2.frame_number];
        end

        function error_val = getError(obj)
            obj.updateSensorDataAll();
            error_val = obj.port_handles(1,1).error;
            if ~obj.port_handles(1,1).hasUsablePosition()
                error_val = 99;
            end
        end

        function sensor_available = isSensorAvailable(obj)
            if obj.device_init == 1
                obj.updateSensorDataAll();
                if ~obj.port_handles(1,1).hasUsablePosition()
                    sensor_available = 0;
                else
                    sensor_available = 1;
                end
            else
                sensor_available = 0;
            end
        end
    end

    methods (Access = public)
        function sendCommand(obj, command)
            if obj.selected_command_format == obj.COMMAND_FORMAT_1
                % Not implemented (CRC & ':' separators)
            else
                fprintf(obj.serial_port, command); % appends CR automatically
            end
        end

        function reply = sendCommandAndGetReply(obj, command)
            obj.sendCommand(command);
            t0 = tic;
            while obj.serial_port.BytesAvailable == 0 && toc(t0) < obj.serial_port.Timeout
                pause(0.01);
            end
            if obj.serial_port.BytesAvailable == 0
                reply = ''; % timeout
                return
            end
            reply = fgetl(obj.serial_port);
        end

        function reply = APIREV(obj), reply = obj.sendCommandAndGetReply('APIREV '); end
        function reply = BEEP(obj, n_beep), reply = obj.sendCommandAndGetReply(sprintf('BEEP %s', n_beep)); end

        function [reply_body, error_checking] = BX(obj, reply_option)
            obj.sendCommand(sprintf('BX %s', reply_option));
            start_sequence = fread(obj.serial_port, 1, 'uint16');
            reply_length   = fread(obj.serial_port, 1, 'uint16');
            header_CRC     = fread(obj.serial_port,  1, 'uint16');
            reply_body     = fread(obj.serial_port, reply_length, 'uint8');
            crc            = fread(obj.serial_port,  1, 'uint16');
            error_checking = [start_sequence; header_CRC; crc];
        end

        function reply = COMM(obj, baud_rate, data_bits, parity, stop_bits, hardware_handshaking)
            reply = obj.sendCommandAndGetReply(sprintf('COMM %s%s%s%s%s', baud_rate, data_bits, parity, stop_bits, hardware_handshaking));
        end

        function reply = ECHO(obj, message), reply = obj.sendCommandAndGetReply(sprintf('ECHO %s', message)); end
        function reply = GET(obj, user_parameter_name), reply = obj.sendCommandAndGetReply(sprintf('GET %s', user_parameter_name)); end
        function reply = INIT(obj), reply = obj.sendCommandAndGetReply('INIT '); end
        function reply = LED(obj, port_handle, led_number, state), reply = obj.sendCommandAndGetReply(sprintf('LED %s%s%s', port_handle, led_number, state)); end
        function reply = PDIS(obj, port_handle), reply = obj.sendCommandAndGetReply(sprintf('PDIS %s', port_handle)); end
        function reply = PENA(obj, port_handle, tool_tracking_priority), reply = obj.sendCommandAndGetReply(sprintf('PENA %s%s', port_handle, tool_tracking_priority)); end
        function reply = PHF(obj, port_handle), reply = obj.sendCommandAndGetReply(sprintf('PHF %s', port_handle)); end
        function reply = PHINF(obj, port_handle, reply_option), reply = obj.sendCommandAndGetReply(sprintf('PHINF %s%s', port_handle, reply_option)); end
        function reply = PHSR(obj, reply_option), reply = obj.sendCommandAndGetReply(sprintf('PHSR %s', reply_option)); end
        function reply = PINIT(obj, port_handle), reply = obj.sendCommandAndGetReply(sprintf('PINIT %s', port_handle)); end
        function reply = PPRD(obj, port_handle, srom_device_address), reply = obj.sendCommandAndGetReply(sprintf('PPRD %s%s', port_handle, srom_device_address)); end
        function reply = PPWR(obj, port_handle, srom_device_address, srom_device_data), reply = obj.sendCommandAndGetReply(sprintf('PPWR %s%s%s', port_handle, srom_device_address, srom_device_data)); end
        function reply = PSEL(obj, port_handle, tool_srom_device_id), reply = obj.sendCommandAndGetReply(sprintf('PSEL %s%s', port_handle, tool_srom_device_id)); end
        function reply = PSOUT(obj, port_handle, gpio_1_state, gpio_2_state, gpio_3_state), reply = obj.sendCommandAndGetReply(sprintf('PSOUT %s%s%s%s', port_handle, gpio_1_state, gpio_2_state, gpio_3_state)); end
        function reply = PSRCH(obj, port_handle), reply = obj.sendCommandAndGetReply(sprintf('PSRCH %s', port_handle)); end
        function reply = PURD(obj, port_handle, user_srom_device_address), reply = obj.sendCommandAndGetReply(sprintf('PURD %s%s', port_handle, user_srom_device_address)); end
        function reply = PUWR(obj, port_handle, user_srom_device_address, user_srom_device_data), reply = obj.sendCommandAndGetReply(sprintf('PUWR %s%s%s', port_handle, user_srom_device_address, user_srom_device_data)); end
        function reply = PVWR(obj, port_handle, start_address, tool_definition_data), reply = obj.sendCommandAndGetReply(sprintf('PUWR %s%s%s', port_handle, start_address, tool_definition_data)); end
        function reply = RESET(obj, reset_option), reply = obj.sendCommandAndGetReply(sprintf('RESET %s', reset_option)); end
        function reply = SFLIST(obj, reply_option), reply = obj.sendCommandAndGetReply(sprintf('SFLIST %s', reply_option)); end
        function reply = TSTART(obj, reply_option), reply = obj.sendCommandAndGetReply(sprintf('TSTART %s', reply_option)); end
        function reply = TSTOP(obj), reply = obj.sendCommandAndGetReply('TSTOP '); end
        function reply = TTCFG(obj, port_handle), reply = obj.sendCommandAndGetReply(sprintf('TTCFG %s', port_handle)); end
        function reply = TX(obj, reply_option), reply = obj.sendCommandAndGetReply(sprintf('TX %s', reply_option)); end
        function reply = VER(obj, reply_option), reply = obj.sendCommandAndGetReply(sprintf('VER %s', reply_option)); end
        function reply = VSEL(obj, volume_number), reply = obj.sendCommandAndGetReply(sprintf('VSEL %s', volume_number)); end
    end

    methods
        function delete(obj)
            if strcmp(obj.serial_port.Status, 'open')
                obj.RESET(obj.RESET_SOFT);
                pause(3);
                obj.closeSerialPort();
            end
            delete(obj.serial_port);
        end
    end

    methods (Access = private)
        function value = readScalar(obj, precision)
            [value, count] = fread(obj.serial_port, 1, precision);
            if count ~= 1
                error('AuroraDriver_2:IncompleteBX', ...
                    'Incomplete BX reply while reading %s.', precision);
            end
        end

        function out = ternary(~, cond, a, b)
            if cond, out = a; else, out = b; end
        end
    end
end
