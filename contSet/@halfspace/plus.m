function [h] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the addition of a vector with a
% halfspace
%
% Syntax:  
%    [h] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - halfspace object or numerical vector
%    summand2 - halfspace object or numerical vector
%
% Outputs:
%    h - halfspace object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      28-August-2013
% Last update:  16-March-2021 (MW, error handling)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing: assign halfspace and summand
if isa(summand1,'halfspace')
    h=summand1;
    summand=summand2;
    
elseif isa(summand2,'halfspace')
    % switch order
    h=summand2;
    summand=summand1;  
end

% error handling
if isempty(h)
    % empty case
    [msg,id] = errEmptySet();
    error(id,msg);
elseif ~isvector(summand)
    % summand not a vector
    [msg,id] = errWrongInput('summand');
    error(id,msg);
elseif dim(h) ~= length(summand)
    % dimension mismatch
    [id,msg] = errDimMismatch();
    error(id,msg);
end
        

% compute Minkowski sum
if isnumeric(summand)
    h.d = h.d + h.c.'*summand;
else
    h = [];
end

%------------- END OF CODE --------------