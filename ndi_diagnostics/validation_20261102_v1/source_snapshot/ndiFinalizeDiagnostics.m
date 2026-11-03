function diagnostics=ndiFinalizeDiagnostics(raw,frame,previousFrame,sampleValid, ...
    calibrationAvailable,orientationValid,pollStart,hostTime)
%NDIFINALIZEDIAGNOSTICS Passive SAME-BX diagnostics; never changes acceptance.
% Schema1:32x2 double, columns base/tip. See ndiDiagnosticsSchema.
%#codegen
diagnostics=raw;
diagnostics(17,:)=double(isfinite(frame));
diagnostics(18,:)=double(isfinite(frame)&frame~=previousFrame);
diagnostics(19,:)=double(all(isfinite(frame))&&frame(1)==frame(2));
diagnostics(20,:)=double(sampleValid);
diagnostics(21,:)=double(calibrationAvailable);
diagnostics(22,:)=double(orientationValid);
diagnostics(31,:)=pollStart;diagnostics(32,:)=hostTime;
for j=1:2
    mask=uint32(0);
    if raw(4,j)==0,mask=bitor(mask,uint32(1));end
    if raw(4,j)==1&&raw(5,j)~=1,mask=bitor(mask,uint32(2));end
    if raw(8,j)==1,mask=bitor(mask,uint32(4));end
    if raw(9,j)==1,mask=bitor(mask,uint32(8));end
    if raw(10,j)==1,mask=bitor(mask,uint32(16));end
    if raw(5,j)==1&&raw(11,j)==0,mask=bitor(mask,uint32(32));end
    if ~isfinite(frame(j)),mask=bitor(mask,uint32(64));end
    if isfinite(frame(j))&&frame(j)==previousFrame(j),mask=bitor(mask,uint32(128));end
    % Mirrors existing pair comparison, including a missing/nonfinite frame.
    if frame(1)~=frame(2),mask=bitor(mask,uint32(256));end
    if all(sampleValid)&&calibrationAvailable&&raw(13,j)==0
        mask=bitor(mask,uint32(512));
    end
    if ~calibrationAvailable,mask=bitor(mask,uint32(1024));end
    if all(sampleValid)&&calibrationAvailable&&all(raw(13,:)==1)&&~orientationValid
        mask=bitor(mask,uint32(2048));
    end
    if ~all(sampleValid),mask=bitor(mask,uint32(4096));end
    diagnostics(23,j)=double(mask);
end
end
