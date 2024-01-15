function I = projectHighDim_(I,N,proj)
% projectHighDim_ - projects an interval set to a higher-dimensional space
%
% Syntax:
%    S = projectHighDim_(S,N,proj)
%
% Inputs:
%    I - interval object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional interval object
%
% Outputs:
%    I - interval object in the higher-dimensional space
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/projectHighDim, interval/project, interval/lift

% Authors:       Tobias Ladner
% Written:       13-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if representsa_(I,'emptySet',eps)
    throw(CORAerror('CORA:notSupported',...
        'Operation projectHighDim is not supported for empty intervals.'));
end

% init bounds
lb = zeros(N,1);
ub = zeros(N,1);

% project
lb(proj) = I.inf;
ub(proj) = I.sup;

% result
I = interval(lb, ub);

% ------------------------------ END OF CODE ------------------------------
