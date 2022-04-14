function ratio = approxVolumeRatio(pZ,varargin)
% approxVolumeRatio - Calculate the approximate ratio of the volumes 
%                     between the dependent generator and the independent  
%                     generator part of the polynomial zonotope 
%                     ratio = (V_ind/V_dep)^(1/n)
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

% Author:       Niklas Kochdumper
% Written:      25-July-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
type = 'interval';
if nargin >= 2
   type = varargin{1}; 
end

% special cases
if isempty(pZ.Grest)
    ratio = 0;
    
elseif isempty(pZ.G)
    ratio = inf;
    
else

    if strcmp(type,'pca')
        % calculate state-space-transformation with pca
        G = [pZ.G -pZ.G pZ.Grest -pZ.Grest];
        [T,~,~] = svd(G);
    
        % transform the polynomial zonotope to the new state space
        pZ = T'*pZ;
    end

    % over-approximate the independent generators part with an interval
    n = length(pZ.c);
    zono = zonotope([zeros(n,1),pZ.Grest]);
    Iind = interval(zono);

    % over-approximate the dependent generators part with an interval
    pZ.Grest = [];
    Idep = interval(pZ);
    
    % remove dimensions that are all-zero
    ind1 = find(rad(Idep) == 0);
    ind2 = find(rad(Iind) == 0);
    
    ind = unique([ind1;ind2]);
    ind = setdiff(1:length(Idep),ind);  

    % calculate the volumes of the parallelotopes
    Vind = volume(Iind(ind));
    Vdep = volume(Idep(ind));

    % calculate the volume ratio
    ratio = (Vind/Vdep)^(1/n);
end

%------------- END OF CODE --------------