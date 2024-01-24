function [res,S] = representsa(S,type,varargin)
% representsa - this function overloads the contSet/representsa function
%     for 'numeric' variables
%
% Syntax:
%    res = representsa(S,type)
%    res = representsa(S,type,tol)
%    [res,S] = representsa(S,type)
%    [res,S] = representsa(S,type,tol)
%
% Inputs:
%    S - numeric vector/matrix
%    type - char array
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%    S - converted set
%
% Example: 
%    S = [];
%    representsa_(S,'emptySet')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa, isempty, contSet/isemptyobject

% Authors:       Tobias Ladner
% Written:       16-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[res,S] = representsa_(S,type,varargin{:});

% ------------------------------ END OF CODE ------------------------------
