function res = isequal(reset1,reset2,varargin)
% isequal - checks if two reset functions have equal pre-/post-state and
%    input dimensions
%
% Syntax:
%    res = isequal(reset1,reset2)
%
% Inputs:
%    reset1 - abstractReset object
%    reset2 - abstractReset object
%
% Outputs:
%    res - true/false
%
% Example: 
%    reset1 = abstractReset(2,1,2);
%    reset2 = abstractReset(2,1,3);
%    isequal(reset1,reset1);
%    isequal(reset1,reset2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
narginchk(2,3);
inputArgsCheck({{reset1,'att','abstractReset'};
                {reset2,'att','abstractReset'}});

% compare number of states and inputs
res = reset1.preStateDim == reset2.preStateDim ...
    && reset1.inputDim == reset2.inputDim ...
    && reset1.postStateDim == reset2.postStateDim;

% ------------------------------ END OF CODE ------------------------------
