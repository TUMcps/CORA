function writeMatrix(fid,M,varName,varargin)
% writeMatrix - writes the content of a symbolic matrix with a custom
%    variable name to a file; for aesthetic reasons, we differentiate
%    between 2D, 3D, and nD cases
%
% Syntax:
%    writeMatrix(fid,M,varName)
%    writeMatrix(fid,M,varName,...)
%
% Inputs:
%    fid - identifier of the file to which the matrix is written
%    M - symbolic nD matrix
%    varName - variable name
%    Name-Value pairs (all optional, arbitrary order):
%       <'BracketSubs',bracketSubsOn> - true/false whether bracketSubs
%           should be called
%       <'Sparse',sparseOn> - true/false whether matrix should be sparse
%       <'IntervalArithmetic',intervalOn> - true/false whether output
%           matrix should be converted to an interval
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: matlabFunction, writeMatrixFile, bracketSubs

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

% read out size of matrix
sz = size(M);

% read input arguments
NVpairs = varargin(1:end);
% check list of name-value pairs
checkNameValuePairs(NVpairs,{'BracketSubs','Sparse','IntervalArithmetic'});
% function handle given?
[NVpairs,bracketSubsOn] = readNameValuePair(NVpairs,'BracketSubs',@islogical,false);
% symbolic function given? (default value determined by size)
[NVpairs,sparseOn] = readNameValuePair(NVpairs,'Sparse',@islogical,numel(sz)>=3);
% interval conversion desired?
[NVpairs,intervalOn] = readNameValuePair(NVpairs,'IntervalArithmetic',@islogical,false);


% different formatting for different dimensions of the matrix
if numel(sz) == 2
    aux_write2D(fid,M,varName,bracketSubsOn,sparseOn,intervalOn);
elseif numel(sz) == 3
    aux_write3D(fid,M,varName,bracketSubsOn,sparseOn,intervalOn);
else
    aux_writenD(fid,M,varName,bracketSubsOn,sparseOn,intervalOn);
end

end


% Auxiliary functions -----------------------------------------------------

function aux_write2D(fid,M,varName,bracketSubsOn,sparseOn,intervalOn)
% 1D, 2D

if bracketSubsOn
    M_char = bracketSubs(char(M));
else
    M_char = char(M);
end
str = replace(M_char,';',[';...' newline]);

% additional text: sparse, then interval
if sparseOn
    str = sprintf('sparse(%s)', str);
end
if intervalOn
    str = sprintf('interval(%s)', str);
end

% write to file
fprintf(fid, '%s = %s;\n\n', varName, str);

end

function aux_write3D(fid,M,varName,bracketSubsOn,sparseOn,intervalOn)
% 3D: read out content of each page

sz = size(M);

% init the matrix by zeros, enclose by sparse/interval if necessary
initStr = sprintf('zeros(%i,%i,%i)', sz);
if sparseOn
    initStr = sprintf('sparse(%s)', initStr);
end
if intervalOn
    initStr = sprintf('interval(%s)', initStr);
end
fprintf(fid, '%s = %s;\n\n', varName, initStr);

% insert value if page is non-zero
for k=1:sz(3)
    % skip line if entries are all-zero (already zero via initialization)
    if all(M(:,:,k) == zeros(sz(1),sz(2)),'all')
        continue
    end

    if bracketSubsOn
        M_char = bracketSubs(char(M(:,:,k)));
    else
        M_char = char(M(:,:,k));
    end
    str = replace(M_char,';',[';...' newline]);
    fprintf(fid, '%s(:,:,%i) = %s;\n\n', varName, k, str);
end

end

function aux_writenD(fid,M,varName,bracketSubsOn,sparseOn,intervalOn)
% nD: use reshape

% sparse and interval currently not supported
if sparseOn || intervalOn
    throw(CORAerror('CORA:notSupported',...
        ['Sparsity and conversion to interval currently not supported '...
        'for matrices larger than 3D.']));
end

if bracketSubsOn
    M_char = bracketSubs(char(M(:)));
else
    M_char = char(M(:));
end
str = sprintf('reshape(%s,[%s]);', M_char, num2str(size(M)));
fprintf(fid, '%s = %s;\n\n', varName, str);

end

% ------------------------------ END OF CODE ------------------------------
