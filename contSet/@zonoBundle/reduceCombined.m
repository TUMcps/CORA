function zB = reduceCombined(zB,option,varargin)
% reduceCombined - Reduces the order of a zonotope bundle by not reducing
%    each zonotope separately, but in a combined fashion
%
% Syntax:
%    zB = reduceCombined(zB,option,varargin)
%
% Inputs:
%    zB - zonoBundle object
%    option - reduction method selector
%    order - desired order
%    filterLength - ???
%
% Outputs:
%    zB - reduced zonotope bundle
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       21-February-2011
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
[order,filterLength] = setDefaultValues({1,[]},varargin);

% check input arguments
inputArgsCheck({{zB,'att','zonoBundle'};
                {order,'att','numeric','nonnan'};
                {filterLength,'att','numeric','nonnan'}});

%option='methC'
if strcmp(option,'methC')
    [zB,t] = reduceMethC(zB,filterLength);
%wrong argument
else
    throw(CORAerror('CORA:wrongValue','second','methC'));
end


% ------------------------------ END OF CODE ------------------------------
