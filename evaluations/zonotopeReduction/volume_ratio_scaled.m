function [vr] = volume_ratio_scaled(Z, Zred)
% volume - Computes the volume ratio of a zonotope Z and Zred
%
% Syntax:  
%    [vol] = volume(Z)
%
% Inputs:
%    Z - zonotope object (original one)
%    Zred - zonotope object (over-approximates Z) 
%
% Outputs:
%    vr - volume ratio, see Alhoff's Diss, page 19, Definition 2.5
%
%

% Author:       Anna Kopetzki
% Written:      28-August-2017
% Last update:  ----
% Last revision:----

%------------- BEGIN CODE --------------

%dimension and nrOfGenerators
G=Z.Z(:,2:end);
[dim,nrOfGen]=size(G);

Gred = Zred.Z(:,2:end);
[dimred, nrOfGenRed] = size(Gred);


W = diag(ones(dim,1));
V = volume(W*Z);
Vred = volume(W*Zred);


vr = (Vred/V)^(1/dim);







%------------- END OF CODE --------------
