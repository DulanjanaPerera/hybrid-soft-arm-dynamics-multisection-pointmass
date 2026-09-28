classdef AuroraReader_s2_Position_only < matlab.System
    % Stream XYZ measurements from two Aurora sensors into Simulink.

    properties (Nontunable)
        Port = 'COM12'
        BaudRate = 230400
    end

    properties(Access = private)
        aurora
        trackingStarted = false
    end

    methods(Access = protected)
        function setupImpl(obj)
            obj.aurora = AuroraDriver_2(obj.Port);
            try
                obj.aurora.openSerialPort();
                obj.aurora.init();
                obj.aurora.BEEP('2');
                obj.aurora.setBaudRate(obj.BaudRate);
                obj.aurora.detectAndAssignPortHandles();
                if obj.aurora.n_port_handles < 2
                    error('NDI:TwoSensorsRequired', ...
                        'Expected two Aurora sensors, found %d.', ...
                        obj.aurora.n_port_handles);
                end
                obj.aurora.initPortHandleAll();
                obj.aurora.enablePortHandleDynamicAll();
                obj.aurora.startTracking();
                obj.trackingStarted = true;
                obj.aurora.BEEP('1');
            catch setupError
                try
                    releaseImpl(obj);
                catch
                    % Preserve the error that caused setup to fail.
                end
                rethrow(setupError);
            end
        end

        function [x1, y1, z1, e1, x2, y2, z2, e2] = stepImpl(obj)
            [x, y, z, ep] = obj.aurora.measureTipPositionAll();

            x1 = x(1); y1 = y(1); z1 = z(1); e1 = ep(1);
            x2 = x(2); y2 = y(2); z2 = z(2); e2 = ep(2);
        end

        function releaseImpl(obj)
            if isempty(obj.aurora)
                return;
            end
            if obj.trackingStarted
                try
                    obj.aurora.stopTracking();
                catch
                    % Continue closing the serial connection.
                end
                obj.trackingStarted = false;
            end
            try
                obj.aurora.BEEP('3');
            catch
                % A disconnected device should not prevent cleanup.
            end
            delete(obj.aurora);
            obj.aurora = [];
        end

        % Output size
        function [sz1, sz2, sz3, sz4, sz5, sz6, sz7, sz8] = getOutputSizeImpl(~)
            sz1 = [1 1]; sz2 = [1 1]; sz3 = [1 1]; sz4 = [1 1];
            sz5 = [1 1]; sz6 = [1 1]; sz7 = [1 1]; sz8 = [1 1];
        end

        % Output data type
        function [dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt8] = getOutputDataTypeImpl(~)
            dt1 = 'double'; dt2 = 'double'; dt3 = 'double'; dt4 = 'double';
            dt5 = 'double'; dt6 = 'double'; dt7 = 'double'; dt8 = 'double';
        end

        % Output complexity
        function [cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8] = isOutputComplexImpl(~)
            cp1 = false; cp2 = false; cp3 = false; cp4 = false;
            cp5 = false; cp6 = false; cp7 = false; cp8 = false;
        end

        % Output fixed size
        function [fs1, fs2, fs3, fs4, fs5, fs6, fs7, fs8] = isOutputFixedSizeImpl(~)
            fs1 = true; fs2 = true; fs3 = true; fs4 = true;
            fs5 = true; fs6 = true; fs7 = true; fs8 = true;
        end
    end
end
