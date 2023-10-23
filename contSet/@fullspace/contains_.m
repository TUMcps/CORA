function res = contains_(fs,S,varargin)
% contains_ - determines if a full-dimensional space contains a set or a
%    point
%    case R^0: only contains R^0, 0 (not representable in MATLAB), and []
%
% Syntax:
%    res = contains_(fs,S)
%    res = contains_(fs,S,type)
%    res = contains_(fs,S,type,tol)
%
% Inputs:
%    fs - fullspace object
%    S - contSet object or numerical vector
%    type - 'exact' or 'approx'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    fs = fullspace(2);
%    p = [1;1];
%    res = contains(fs,p);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (MW, rename contains_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if fs.dimension == 0
    throw(CORAerror('CORA:notSupported',...
        'Containment check for R^0 not supported'));
end

% full-dimensional space contains all other sets, including itself
res = true;

% ------------------------------ END OF CODE ------------------------------
