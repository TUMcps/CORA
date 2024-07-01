function [sys,params,spec] = json2cora_linearSys(filename)
% json2cora_linearSys - read a verification benchmark from the JSON format
%   specified in benchmark2json
%
% Syntax:
%    [sys,params,spec] = json2cora_linearSys(filepath)
%
% Inputs:
%    filename - (optional) name of .json file, incl. extension
%
% Outputs:
%    sys - linearSys object
%    params - model parameters
%    spec - safety specifications

% Authors:       Mark Wetzlinger
% Written:       07-June-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read raw data from file
try
    fid = fopen(filename,'r');
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
catch
    throw(CORAerror('CORA:specialError','Could not read file'));
end

% convert JSON to struct
S = jsondecode(str);

% instantiate linear system
sys = linearSys(S.A,S.B,[],S.C);

% instantiate model parameters
params.R0 = interval(S.X0.lowerbound,S.X0.upperbound);
params.U = interval(S.U.lowerbound,S.U.upperbound);
params.tFinal = S.tend;

% instantiate specifications
spec = specification();
for i=1:numel(S.unsafeSet)
    spec(i,1) = specification(...
        interval(S.unsafeSet(i).lowerbound,S.unsafeSet(i).upperbound),...
        'unsafeSet',...
        interval(0,params.tFinal));
end

% ------------------------------ END OF CODE ------------------------------
