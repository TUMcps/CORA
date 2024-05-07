function [zSq,zH] = dependentTerms(matZ,r)
% dependentTerms - computes exact Taylor terms of a matrix zonotope square
%    and a matrix zonotope exponential
%
%    These different tasks are computed in one m-file to save computation
%    time: the for loop has to be executed only once and help functions do
%    not have to be called so often
%
% Syntax:
%    [zSq,zH] = dependentTerms(matZ,r)
%
% Inputs:
%    matZ - matZonotope object
%    r - time step size
%
% Outputs:
%    zSq - exact square matrix
%    zH - exact Taylor terms up to second order
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       24-September-2010
% Last update:   26-April-2024 (TL, major speed using new matZonotope class)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%load data from object structure
C=matZ.C;
G=matZ.G;
n=dim(matZ);
gens=matZ.numgens;


%square computation--------------------------------------------------------
%new center
sqC = C^2*r^2 + 0.5*sum(pagemtimes(G,G),3)*r^2;

%get generators
sqG = cat(3, ...
    (pagemtimes(C,G) + pagemtimes(G,C))*r^2, ...
    0.5*pagemtimes(G,G)*r^2 ...
);
%get indices for 3rd set of generators
if (gens>=2)
    ind = combinator(gens,2,'c');
    sqG = cat(3, ...
            sqG, ...
            ( ...
                pagemtimes( ...
                    G(:,:,ind(:,1)), ...
                    G(:,:,ind(:,2))) ...
                + ...            
                pagemtimes( ...
                    G(:,:,ind(:,2)), ...
                    G(:,:,ind(:,1))) ...
            ) * r^2 ...
    );
end
%--------------------------------------------------------------------------


%H computation-------------------------------------------------------------
%new center
HC = eye(n) + C*r + 0.5*sqC;

%get generators
HG = 0.5 * sqG;
HG(:,:,1:gens) = HG(:,:,1:gens) + G*r;

%--------------------------------------------------------------------------

%write as matrix zonotopes
zSq=matZonotope(sqC,sqG);
zH=matZonotope(HC,HG);

% ------------------------------ END OF CODE ------------------------------
