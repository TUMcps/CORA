function Z = exactPlus(Z1,Z2,varargin)
% exactPlus - adds two zonotopes by adding all generators of common
%    generator factors. Caution: It has to be ensured from outside that the
%    generator factors match!
%
% Syntax:  
%    Z = exactPlus(Z1,Z2)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%    nrOfGens - (optional) limit on the nr of generators that can be
%               added exactly
%
% Outputs:
%    Z - final zonotope
%
% Example: 
%    Z1 = zonotope([0;0],[1 2 -1; 1 -1 3]);
%    Z2 = zonotope([0;0],[2 4 -1; 2 -2 3]);
%    nrOfGens = 2;
%    Zexact = exactPlus(Z1,Z2);
%    Z = Z1 + Z2;
% 
%    plot(Z1); hold on;
%    plot(Z2);
%    plot(Zexact,[1,2],'g');
%    plot(Z,[1,2],'r');
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

%obtain matrix
Zmat1 = Z1.Z;
Zmat2 = Z2.Z;

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