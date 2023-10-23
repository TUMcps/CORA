function res = ne(trans1,trans2,varargin)
% ne - overloads '~='-operator for transitions
%
% Syntax:
%    res = trans1 ~= trans2
%    res = ne(trans1,trans2)
%    res = ne(trans1,trans2,tol)
%
% Inputs:
%    trans1 - transition object
%    trans2 - transition object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    % guard set
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%
%    % reset function
%    reset1 = struct('A',[1,0;0,-0.75],'c',[0;0]);
%    reset2 = struct('A',[1,0;0,-0.75],'c',[1;0]);
%
%    % transition
%    trans1 = transition(guard,reset1,1);
%    trans2 = transition(guard,reset2,1);
%
%    % comparison
%    res = trans1 ~= trans2
%
% Other m-files required: transition
% Subfunctions: none
% MAT-files required: none
%
% See also: transition/isequal

% Authors:       Mark Wetzlinger
% Written:       10-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(trans1,trans2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
