function linReset_ = lift(linReset,N,M,stateBind,inputBind,id)
% lift - lifts a linear reset function to a higher-dimensional space
%
% Syntax:
%    linReset_ = lift(linReset,N,M,stateBind,inputBind,id)
%
% Inputs:
%    linReset - linearReset object
%    N - dimension of the higher-dimensional state space
%    M - dimension of higher-dimensional input space
%    stateBind - states of the high-dimensional space that correspond to
%                the states of the low-dimensional reset object
%    inputBind - inputs of the high-dimensional space that correspond to
%                the inputs of the low-dimensional reset object
%    id - true/false whether identity reset function should be used for all
%         other states
%
% Outputs:
%    linReset_ - lifted linearReset object
%
% Example:
%    A = [1 2; 0 -1]; B = [2 0 1; -1 0 0]; c = [1; -5];
%    linReset = linearReset(A,B,c);
%    N = 6; stateBind = [2,3];
%    M = 5; inputBind = [2,3,4];
%    id = true;
%    linReset_ = lift(linReset,N,M,stateBind,inputBind,id);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearReset/lift

% Authors:       Mark Wetzlinger
% Written:       07-September-2024
% Last update:   10-October-2024 (MW, support input matrix)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% note: reset function needs to map from R^n -> R^n
if linReset.preStateDim ~= linReset.postStateDim
    throw(CORAerror('CORA:notSupported',...
        ['Projection of reset functions to higher-dimensional spaces '...
        'only supported for R^n -> R^n.']));
end

% init A matrix by identity or zeros
if id
    Aproj = eye(N);
else
    Aproj = zeros(N);
end
Bproj = zeros(N,M);
cProj = zeros(N,1);

% insert state and input mappings A and B, and vector c
Aproj(stateBind,stateBind) = linReset.A;
Bproj(stateBind,inputBind) = linReset.B;
cProj(stateBind) = linReset.c;

% construct resulting reset object
linReset_ = linearReset(Aproj,Bproj,cProj);

% ------------------------------ END OF CODE ------------------------------
