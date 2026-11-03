classdef NdiMockTransport < handle
    % Offline byte transport; never creates or opens a serial/DAQ device.
    properties
        Status='closed'
        BaudRate=9600
        Timeout=2
        Packet=uint8([])
        Buffer=uint8([])
        Offset=1
        Ascii=''
        BxCount=0
        Commands={}
    end
    properties (Dependent)
        BytesAvailable
    end
    methods
        function value=get.BytesAvailable(obj)
            value=numel(obj.Ascii)+max(0,numel(obj.Buffer)-obj.Offset+1);
        end
        function stage(obj,bytes)
            obj.Packet=uint8(bytes(:)');
        end
        function fprintf(obj,command)
            obj.Commands{end+1}=command;
            if startsWith(command,'BX ')
                obj.BxCount=obj.BxCount+1;obj.Buffer=obj.Packet;obj.Offset=1;
            elseif startsWith(command,'PHSR')
                obj.Ascii='020A0310B031';
            else
                obj.Ascii='OKAY';
            end
        end
        function line=fgetl(obj),line=obj.Ascii;obj.Ascii='';end
        function flushinput(obj),obj.Ascii='';obj.Buffer=uint8([]);obj.Offset=1;end
        function fopen(obj),obj.Status='open';end
        function fclose(obj),obj.Status='closed';end
        function [value,count]=fread(obj,n,precision)
            assert(n==1);
            switch precision
                case 'uint8',width=1;kind='uint8';
                case 'uint16',width=2;kind='uint16';
                case 'uint32',width=4;kind='uint32';
                case 'float32',width=4;kind='single';
                otherwise,error('Unexpected mock precision.');
            end
            if obj.Offset+width-1>numel(obj.Buffer),value=[];count=0;return;end
            bytes=obj.Buffer(obj.Offset:obj.Offset+width-1);obj.Offset=obj.Offset+width;
            value=double(typecast(bytes,kind));count=1;
        end
    end
end
