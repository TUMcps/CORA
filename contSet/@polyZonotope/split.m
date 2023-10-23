function [pZsplit] = split(pZ,varargin)
% split - splits a polynomial zonotope into two or more polynomial
%    zonotopes that jointly enclose that polynomial zonotope
%
% Syntax:
%    [pZsplit] = split(pZ)
%    [pZsplit] = split(pZ, gen)
%
% Description:
%    If only one input is provided, all possible splits of the
%    parallelotope that over-approximates the polynomial zonotope are
%    calculated. If the index of one specific generator is passed as a
%    second input argument, then the split at this specific generator is
%    calculated.
%
% Inputs:
%    pZ - polyZonotope object
%    gen - generator of the over-approximating parallelotope that is
%          splitted
%
% Outputs:
%    pZsplit - cell array of parallotopes represented as polyZonotopes
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 2;0 2 2],[0;0],[1 0 3;0 1 1]);
%    pZsplit = split(pZ);
%    
%    figure; hold on; xlim([-5,5]); ylim([-8,8]);
%    plot(pZsplit{1}{1},[1,2],'FaceColor','b');
%    plot(pZsplit{1}{2},[1,2],'FaceColor','g');
%    plot(pZ,[1,2],'FaceColor','r');
%
%    figure; hold on; xlim([-5,5]); ylim([-8,8]);
%    plot(pZsplit{2}{1},[1,2],'FaceColor','y');
%    plot(pZsplit{2}{2},[1,2],'FaceColor','c');
%    plot(pZ,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: splitLongestGen, splitDepFactor, zonotope/split

% Authors:       Niklas Kochdumper
% Written:       28-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% calculate a parallelotope that over-approximates the polynomial zonotope
zono = zonotope(pZ);
P = reduce(zono,'methC',1);


% Case 1: calculate all possible splits
if nargin == 1
    
    % split all parallelotope generators
    n = length(center(P));
    pZsplit = cell(n,1);    
    
    for iDim = 1:n
        pZsplit{iDim} = aux_splitOneDim(P,iDim); 
    end
    
% Case 2: split at the specified generator
elseif nargin == 2
    
    gen = varargin{1};
    pZsplit = aux_splitOneDim(P,gen);
    
end

end


% Auxiliary functions -----------------------------------------------------

function Zsplit = aux_splitOneDim(Z,splitDim)

    % center and generator matrix
    c = center(Z);
    G = generators(Z);

    % compute centers of splitted parallelotope
    c1 = c-G(:,splitDim)/2;
    c2 = c+G(:,splitDim)/2;

    % compute new set of generators
    Gnew = G;
    Gnew(:,splitDim) = Gnew(:,splitDim)/2;

    % generate splitted parallelpipeds
    Zsplit{1} = polyZonotope(c1,Gnew,[]);
    Zsplit{2} = polyZonotope(c2,Gnew,[]);   
    
end

% ------------------------------ END OF CODE ------------------------------
