function [Z] = exactPlus(Z1,Z2,varargin)
% exactPlus - Adds two quadZonotopes by adding all generators of common
% generator factors. It has to be ensured from outside % that the generator 
% factors match
%
% Syntax:  
%    [Z] = exactPlus(Z1,Z2)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%    nrOfGens - limit on the nr of generators that can be added exactly
%
% Outputs:
%    Z - final zonotope
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      30-August-2013
% Last update:  06-September-2013
% Last revision:---

%------------- BEGIN CODE --------------

%obtain matrice
Zmat1 = get(Z1,'Z');
Zmat2 = get(Z2,'Z');

%number of vectors
nrOfVecs1 = length(Zmat1(1,:));
nrOfVecs2 = length(Zmat2(1,:));

if nargin == 2
    maxVecs = min([nrOfVecs1, nrOfVecs2]);
elseif nargin == 3
    maxVecs = min([nrOfVecs1, nrOfVecs2, varargin{1}+1]);
end

Zmat_front = Zmat1(:,1:maxVecs) + Zmat2(:,1:maxVecs);
Zmat_rest_1 = Zmat1(:,maxVecs+1:end);
Zmat_rest_2 = Zmat2(:,maxVecs+1:end);

Z = zonotope([Zmat_front, Zmat_rest_1, Zmat_rest_2]);

%------------- END OF CODE --------------