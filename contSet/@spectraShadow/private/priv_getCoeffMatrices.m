function [A0,Ai] = priv_getCoeffMatrices(SpS)
% priv_getCoeffMatrices - returns a 'separated' list of matrices that make
%    up big coefficient matrix A required by the constructor
%
% Syntax:
%    [A0,Ai] = priv_getCoeffMatrices(S)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    A0 - matrix of the constant coefficients
%    Ai - cell array of the other coefficient matrices
%
%    A0 and Ai can then represent A as follows:
%    A = [A0 Ai{1} Ai{2} ... Ai{m}];
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       01-August-2023 
% Last update:   ---    
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

A = SpS.A;
if isempty(A)
    A0 = [];
    Ai = {};
    return
end

k = size(A,1);
m = size(A,2)/k-1;

A0 = A(:,1:k);

Ai = cell([1 m]);
for i=1:m
    Ai{i} = A(:,i*k+1:(i+1)*k);
end

% ------------------------------ END OF CODE ------------------------------
