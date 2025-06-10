function res = contains(matZ,S,varargin)
% contains - checks if the given points are within the set
%
% Syntax:
%    res = contains(matZ,S)
%    res = contains(matZ,S,type)
%    res = contains(matZ,S,type,tol)
%
% Inputs:
%    matZ - matZonotope object
%    S - numeric
%    method - str, 'exact', 'approx'
%    tol - numeric, tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    C = [0 0; 0 0];
%    G(:,:,1) = [1 3; -1 2]; G(:,:,2) = [2 0; 1 -1];
%    matZ = matZonotope(C,G);
%    res = contains(matZ,C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/contains

% Authors:       Tobias Ladner
% Written:       20-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default values
[type,tol] = setDefaultValues({'exact',1e-12},varargin);

% check input arguments
inputArgsCheck({{matZ,'att','matZonotope'},...
                {S,'att','numeric','nonempty'},...
                {type,'str',{'exact','approx'}},...
                {tol,'att','numeric',{'scalar','nonnegative'}}});

% use zonotope/contains
[n,m,p] = size(S);
res = contains(zonotope(matZ),reshape(S,n*m,p));

end

% ------------------------------ END OF CODE ------------------------------
