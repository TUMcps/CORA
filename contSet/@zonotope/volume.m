function vol = volume(Z,type,o)
% volume - computes the volume of a zonotope according to page 40 in [1]
%
% Syntax:  
%    vol = volume(Z,type,o)
%
% Inputs:
%    Z - zonotope object
%    type - type of approximation 
%           - 'none' (default)
%           - 'reduce' for reduced zonotope with order o
%           - 'alamo', see [2]
%
% Outputs:
%    vol - volume
%
% Example: 
%    Z=zonotope([1 -1 0; 0 0 -1]);
%    vol=volume(Z)
%
% References:
%    [1] E. Grover et al. "Determinants and the volumes of parallelotopes 
%        and zonotopes", 2010 
%    [2] Alamo et al. "Bounded error identification of systems with 
%        time-varying parameters§, TAC 2006.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      24-August-2007 
% Last update:  19-July-2010
%               02-September-2019 (incl. approx)
%               04-May-2020 (MW, add vol=0 cases)
%               09-September-2020 (MA, Alamo approx added, reduce changed)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin < 2
    type = 'none';
elseif nargin < 3
    r = 0;   
end

%dimension and nrOfGenerators
G=generators(Z);
[n,nrOfGens]=size(G);

if nrOfGens < n || rank(G) < n
    vol = 0;
    return
end

% exact computation
if strcmp(type,'none')
    % exact calculation

    %possible combinations of n=dim generators from all generators
    comb = combinator(nrOfGens,n,'c');
    nrOfComb=length(comb(:,1));

    accVol = 0;

    for i=1:nrOfComb
        try
            currVol = abs(det(G(:,comb(i,:))));
            accVol = accVol + currVol;
        catch
            currVol=0;
            disp('parallelogram volume could not be computed');
        end
    end

    vol=2^n*accVol;

% over-approximative volume using order reduction
elseif strcmp(type,'reduce')
    % reduce zonotope
    Zred = reduce(Z,'pca',o); 
    vol = volume(Zred);
  
% approximation according to [2]    
elseif strcmp(type,'alamo')
    vol = 2^n*sqrt(det(G*G')); 
end


%------------- END OF CODE --------------