function Zbundle = generateRandom(varargin)
% generateRandom - Generates a random zonoBundle
%
% Syntax:  
%    Zbundle = generateRandom()
%    Zbundle = generateRandom(dim,zons)
%
% Inputs:
%    dim - (optional) dimension
%    zons - (optional) number of zonotopes
%
% Outputs:
%    Zbundle - random zonoBundle
%
% Example: 
%    Zbundle = zonoBundle.generateRandom();
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % bounds for parameters
    dim_low = 1;
    dim_up = 5;

    zons_low = 2;
    zons_up = 5;

    % randomly compute parameters
    zons = zons_low + floor(rand(1) * (zons_up - zons_low + 1));
    n = dim_low + floor(rand(1) * (dim_up - dim_low + 1));

    % parse input arguments
    if nargin >= 1 && ~isempty(varargin{1})
       n = varargin{1};
    end

    if nargin >= 2 && ~isempty(varargin{2})
       zons = varargin{2};  
    end

    % construct random zonotope bundle
    Z = cell(zons,1);
    Z{1} = zonotope.generateRandom(n);
    for z=2:zons
        inside = false;
        while ~inside
            inside = true;
            cen = randPoint(Z{1});
            for i=2:z-1
                if ~in(Z{i},cen)
                    inside = false; continue;
                end
            end
        end
        Z{z} = zonotope.generateRandom(n,cen);
    end

    % instantiate interval
    Zbundle = zonoBundle(Z);

%------------- END OF CODE --------------