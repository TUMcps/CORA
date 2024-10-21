function handle = writeMatrixFile(M,path,fname,varargin)
% writeMatrixFile - creates an m-file with a given name on a given path
%    that outputs an array of nD matrices; the generated file looks like:
%    "   
%    function [out1,out2,...] = fname(in1,in2,...)
%    
%    out1 = [...];
%    out2 = [...];
%    ...
%
%    end
%    "
%    where the output variables out1, out2, ... arbitrarily depend on the
%    input variables in1, in2, ... 
%
% Syntax:
%    writeMatrixFile(M,path,fname)
%    writeMatrixFile(M,path,fname,...)
%       <'VarNamesIn',varNamesIn> - variable names for input arguments
%       <'VarNamesOut',varNamesOut> - variable names for output arguments
%       <'BracketSubs',bracketSubsOn> - true/false whether bracketSubs
%           should be called
%       <'IntervalArithmetic',intervalOn> - true/false whether output
%           matrix should be converted to an interval
%       <'Sparse',sparseOn> - true/false whether matrix should be sparse
%
% Inputs:
%    M - cell array of matrices
%    path - path where the function should be created
%    fname - function name for the file computing the Jacobian
%    varNamesIn - cell array with variable names for file input arguments
%    varNamesOut - cell array with variable names for file output arguments
%    bracketSubsOn - true/false whether bracketSubs should be called
%
% Outputs:
%    handle - function handle to file
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: derive

% Authors:       Mark Wetzlinger
% Written:       13-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(3,Inf);
% name-value pairs -> number of input arguments is always a multiple of 2+1
if mod(nargin,2) ~= 1
    throw(CORAerror('CORA:oddNumberInputArgs'));
end

inputArgsCheck({{M,'att','cell'};...
                {path,'att','char'};...
                {fname,'att','char'}});

% read input arguments
NVpairs = varargin(1:end);
% check list of name-value pairs
checkNameValuePairs(NVpairs,{'VarNamesIn','VarNamesOut',...
    'IntervalArithmetic','BracketSubs','Sparse'});
% interval arithmetic true/false?
[NVpairs,bracketSubsOn] = readNameValuePair(NVpairs,'BracketSubs',@islogical,false);
% interval arithmetic true/false?
[NVpairs,intervalOn] = readNameValuePair(NVpairs,'IntervalArithmetic',@islogical,false);
% verbose output?
[NVpairs,sparseOn] = readNameValuePair(NVpairs,'Sparse',@islogical,false);

% names for input/output arguments given?
[NVpairs,varNamesIn] = readNameValuePair(NVpairs,'VarNamesIn');
[NVpairs,varNamesOut] = readNameValuePair(NVpairs,'VarNamesOut');
% ...their default values are a bit more complicated
[varNamesIn,varNamesOut] = aux_setDefaultValues(varNamesIn,varNamesOut);

% ensure that matrices and variable names are of equal length
numVars = numel(M);
if numVars ~= numel(varNamesOut)
    throw(CORAerror('CORA:wrongValue','second',...
        'Number of variable names must match number of matrices.'));
end

% open file
fid = fopen([path filesep fname '.m'],'w');

% try-catch to ensure that file will be closed
try
    % write first line
    fprintf(fid, 'function [%s] = %s(%s)\n\n', ...
        strjoin(varNamesOut,','), fname, strjoin(varNamesIn,','));

    % write matrices to file
    arrayfun(@(i) writeMatrix(fid,M{i},varNamesOut{i},...
        'BracketSubs',bracketSubsOn,...
        'IntervalArithmetic',intervalOn,...
        'Sparse',sparseOn),...
        1:numVars);

    % properly end file
    fprintf(fid, 'end\n');

catch ME
    % close file
    fclose(fid);
    rethrow(ME)
end

% close file
fclose(fid);

% output argument
handle = eval(['@', fname]);

end


% Auxiliary functions -----------------------------------------------------

function [varNamesIn,varNamesOut] = aux_setDefaultValues(varNamesIn,varNamesOut)

% by default, the number of input arguments is 1
if isempty(varNamesIn)
    varNamesIn = {'in1'};
end

% by default, the number of output arguments is equal to the number of
% Jacobians
if isempty(varNamesIn)
    varNamesOut = arrayfun(@(i) sprintf('out%i',i),1:numel(J),...
        'UniformOutput',false);
end

end

% ------------------------------ END OF CODE ------------------------------
