function [p] = randPointExtreme(pZ)
% randPointExtreme - generates a random point close to the border of a 
%                    polynomial zonotope
%
% Syntax:  
%    [p] = randPointExtreme(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    p - random extreme point in R^n
%
% Example: 
%    pZ = polyZonotope([0;0], [2 0 1;1 2 1],[],[1 0 1;0 1 3]);
%    
%    p = randPointExtreme(pZ);
%
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(p(1),p(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: randPoint

% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % center
    p = pZ.c;

    % Part 1: dependent generators
    if ~isempty(pZ.G)

        % genrator factors randomly in the interval [-1,1]
        N = size(pZ.expMat,1);
        beta = sign(rand(N,1)*2 - ones(N,1));

        fact = prod(beta.^pZ.expMat,1);
        p = p + pZ.G * fact';
    end

    % Part 2: independent generators
    if ~isempty(pZ.Grest)

        % genrator factors randomly in the interval [-1,1]
        N = size(pZ.Grest,2);
        alpha = sign(rand(N,1)*2 - ones(N,1));

        p = p + pZ.Grest * alpha;
    end

%------------- END OF CODE --------------