function res = benchmark2json(sys,params,spec,varargin)
% benchmark2json - write a verification benchmark to JSON format
%    format for ARCH24:
%    {
%       "version":"1.0",
%       "A": [[8,1,6],[3,5,7],[4,9,2]],
%       "B": [[1],[2],[3]],
%       "C": [[1,2,3]],
%       "X0":{"type":"interval","lowerbound":[-1,-2,-3],"upperbound":[3,4,5]},
%       "U":{"type":"interval","lowerbound":[-0.05],"upperbound":[0.05]},
%       "tend":2,
%       "unsafeSet":[
%         {"type":"interval","lowerbound":[-1],"upperbound":[1]},
%         {"type":"interval","lowerbound":[-5],"upperbound":[-3]}
%         ]
%    }
%
% Syntax:
%    res = benchmark2json(sys,params,spec)
%    res = benchmark2json(sys,params,spec,filename)
%
% Inputs:
%    sys - linearSys object
%    params - model parameters
%    spec - safety specifications
%    filename - (optional) name of .json file, incl. extension
%
% Outputs:
%    res - JSON file successfully generated

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       07-June-2024
% Last update:   28-June-2024 (TL, superior implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default filename unless given
filename = setDefaultValues({'verificationBenchmark.json'},varargin);

% conversion only supported for restrictions used at ARCH24: R0, U must be
% intervals, all specs must be intervals and unsafe sets
if ~isa(params.R0,'interval')
    throw(CORAerror('CORA:notSupported','Only intervals supported for params.R0.'));
end
if ~isa(params.U,'interval')
    throw(CORAerror('CORA:notSupported','Only intervals supported for params.U.'));
end
for s=1:length(spec)
    if ~isa(spec(s).set,'interval') || ~strcmp(spec(s).type,'unsafeSet')
        throw(CORAerror('CORA:notSupported',...
            'All specifications must be unsafe sets and intervals.'));
    end
end

% init struct
S = struct;

% version number
S.version = "1.0";

% dynamics
S.A = sys.A;
S.B = sys.B;
S.C = sys.C;

% initial set
S.X0 = aux_set2struct(params.R0);

% input set
S.U = aux_set2struct(params.U);

% time horizon
S.tend = params.tFinal;

% specifications
S.unsafeSet = aux_spec2struct(spec);

% convert to JSON
jsonstring = jsonencode(S,'PrettyPrint',true);

% write to root CORA folder
try
    fid = fopen([CORAROOT filesep filename],'w');
    fprintf(fid,'%s',jsonstring);
    fclose(fid);
catch ME
    throw(CORAerror('CORA:specialError','Unable to open/write to/close file.'));
end

end


% Auxiliary functions -----------------------------------------------------

function S = aux_set2struct(set)

if isa(set,'interval')
    S.type = 'interval';
    S.lowerbound = set.inf;
    S.upperbound = set.sup;
else
    throw(CORAerror('CORA:notSupported','Only interval supported'));
end

end

function S = aux_spec2struct(spec)

for i=1:numel(spec)
    S(i) = aux_set2struct(spec(i).set);
end

end

% ------------------------------ END OF CODE ------------------------------
