function pgon = polygon(pZ,varargin)
% polygon - creates a polygon enclosure of a 2-dimensional polyZonotope
%
% Syntax:  
%    pgon = polygon(pZ)
%    pgon = polygon(pZ,splits)
%
% Inputs:
%    pZ - polyZonotope object
%    splits - number of splits for refinement (optional)
%
% Outputs:
%    pgon - polygon object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    pgon = polygon(pZ,8);
%
%    plot(pgon);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Author:       Niklas Kochdumper
% Written:      08-April-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    splits = 8;

    if nargin > 1 && ~isempty(varargin{1})
       splits = varargin{1}; 
    end
    
    % check if polynomial zonotope is two dimensional
    if dim(pZ) ~= 2
       error('Method "polygon" is only applicable for 2D polyZonotopes!');
    end

    % split the polynomial zonotope multiple times to obtain a better 
    % over-approximation of the real shape
    pZsplit{1} = pZ;
    pZfin = {};

    for i=1:splits
        pZnew = [];
        for j=1:length(pZsplit)
            try
                res = splitLongestGen(pZsplit{j});
                pZnew{end+1} = res{1};
                pZnew{end+1} = res{2};
            catch
                pZfin{end+1} = pZsplit{j};
            end
            
        end
        pZsplit = pZnew;
    end
    
    pZsplit = [pZsplit,pZfin];

    % over-approximate all splitted sets with zonotopes, convert them to 2D
    % polytopes (Matlab build-in) and compute the union
    warOrig = warning;
    warning('off','all');
    pgon = [];

    for i = 1:length(pZsplit)

        % zonotope over-approximation
        zono = zonotope(pZsplit{i});

        % calculate vertices of zonotope
        vert = vertices(zono);

        % tranform to 2D polytope
        pgonTemp = polygon(vert(1,:),vert(2,:));

        % calculate union with previous sets
        pgon = pgonTemp | pgon; 
    end
    
    warning(warOrig);
end

%------------- END OF CODE --------------