function ratio = approxVolumeRatio(pZ,varargin)
% approxVolumeRatio - Calculates the approximate ratio of the volumes 
%    between the polynomial zonotope constructed by only the dependent
%    generators of the given polynomial zonotope and the zonotope
%    constructed by the independent generator part of the polynomial
%    zonotope; the ratio is computed as ratio = (V_ind/V_dep)^(1/n)
%
% Syntax:
%    ratio = approxVolumeRatio(pZ)
%    ratio = approxVolumeRatio(pZ,type)
%
% Inputs:
%    pZ - polyZonotope object
%    type - method used to calculate the volume ('interval' or 'pca')
%
% Outputs:
%    ratio - appoximate volume ratio (V_ind/V_dep)^(1/n)
%
% Example:
%    pZ = polyZonotope([0;0],[1 0 1;1 2 -2],[-1 0.3;1.2 0.2],[1 0 1;0 1 2]);
%    ratio = approxVolumeRatio(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/volume, interval/volume

% Authors:       Niklas Kochdumper
% Written:       25-July-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
type = setDefaultValues({'interval'},varargin);

% check input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
                {type,'str',{'interval','pca'}}});

% special cases
if isempty(pZ.GI)
    ratio = 0;
    
elseif isempty(pZ.G)
    ratio = inf;
    
else

    if strcmp(type,'pca')
        % calculate state-space-transformation with pca
        G = [pZ.G -pZ.G pZ.GI -pZ.GI];
        [T,~,~] = svd(G);
    
        % transform the polynomial zonotope to the new state space
        pZ = T'*pZ;
    end

    % over-approximate the independent generators part with an interval
    n = length(pZ.c);
    zono = zonotope([zeros(n,1),pZ.GI]);
    Iind = interval(zono);

    % over-approximate the dependent generators part with an interval
    pZ.GI = [];
    Idep = interval(pZ);
    
    % remove dimensions that are all-zero
    ind1 = find(rad(Idep) == 0);
    ind2 = find(rad(Iind) == 0);
    
    ind = unique([ind1;ind2]);
    ind = setdiff(1:length(Idep),ind);  

    % calculate the volumes of the parallelotopes
    Vind = volume_(Iind(ind));
    Vdep = volume_(Idep(ind));

    % calculate the volume ratio
    ratio = (Vind/Vdep)^(1/n);
end

% ------------------------------ END OF CODE ------------------------------
