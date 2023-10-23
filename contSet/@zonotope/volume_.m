function vol = volume_(Z,method,order,varargin)
% volume_ - computes the volume of a zonotope
%
% Syntax:
%    vol = volume_(Z,method,order)
%
% Inputs:
%    Z - zonotope object
%    method - (optional) method for approximation 
%           - 'exact' (default)
%           - 'reduce' for reduced zonotope with order o
%           - 'alamo', see [2]
%    order - (optional) zonotope order for reduction before computation
%
% Outputs:
%    vol - volume
%
% Example: 
%    Z = zonotope([1 -1 0; 0 0 -1]);
%    vol = volume(Z)
%
% References:
%    [1] E. Grover et al. "Determinants and the volumes of parallelotopes 
%        and zonotopes", 2010 
%    [2] Alamo et al. "Bounded error identification of systems with 
%        time-varying parameters", TAC 2006.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Matthias Althoff
% Written:       24-August-2007 
% Last update:   19-July-2010
%                02-September-2019 (incl. approx)
%                04-May-2020 (MW, add vol=0 cases)
%                09-September-2020 (MA, Alamo approx added, reduce changed)
%                27-July-2022 (ME, included batchCombinator)
%                18-August-2022 (MW, include standardized preprocessing)
% Last revision: 27-March-2023 (MW, rename volume_)

% ------------------------------ BEGIN CODE -------------------------------

% dimension and nrOfGenerators
G = Z.G;
[n,nrOfGens] = size(G);

% non full-dimensional set
if nrOfGens < n || rank(G) < n
    vol = 0; return
end

% exact computation
if strcmp(method,'exact')
    % exact calculation

    % number of combinations processed at once
    batch_size = 1e7;
    comb_state = struct;

    accVol = 0;

    while true
        
        [batch, comb_state] = batchCombinator(nrOfGens, int16(n), batch_size, comb_state);

        for j=1:length(batch(:,1))
            try
                currVol = abs(det(G(:,batch(j,:))));
                accVol = accVol + currVol;
            catch
                throw(CORAerror('CORA:specialError',...
                'Parallelotope volume could not be computed.'));            
            end
        end

        if comb_state.done == true
            break;
        end

    end

    % multiply result by factor
    vol=2^n*accVol;

% over-approximative volume using order reduction
elseif strcmp(method,'reduce')
    % reduce zonotope
    Zred = reduce(Z,'pca',order); 
    vol = volume(Zred);
  
% approximation according to [2]    
elseif strcmp(method,'alamo')
    vol = 2^n*sqrt(det(G*G')); 
end


% ------------------------------ END OF CODE ------------------------------
