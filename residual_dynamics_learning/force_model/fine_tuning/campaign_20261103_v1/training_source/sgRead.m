function T=sgRead(path,varargin)
%SGREAD Explicit CSV dialect; Windows path columns confuse auto-detection.
T=readtable(path,'Delimiter',',','VariableNamingRule','preserve',varargin{:});
end
