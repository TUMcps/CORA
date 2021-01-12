function vol = volume(Z,varargin)
% volume - Computes the volume of a zonotope according to page 40 in [1]
%
% Syntax:  
%    vol = volume(Z,r)
%
% Inputs:
%    Z - zonotope object
%    r - number of filtered generators in addition to minimal requirement
%           (full rank) for approximative calculation
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
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      24-August-2007 
% Last update:  19-July-2010
%               ??-May-2016 (AK, accumulator)
%               02-Sep-2019 (incl. approx)
%               04-May-2020 (MW, add vol=0 cases)
%               12-January-2021 (MW, clarification of approx case)
% Last revision:---

%------------- BEGIN CODE --------------

r = [];
approx = false;
if nargin == 2
    approx = true;
    if isscalar(varargin{1}) && varargin{1} >= 0 && abs(mod(varargin{1},1)) < eps
        r = varargin{1};
    else
        error("Number of filtered generators invalid!");
    end
end

%dimension and nrOfGenerators
G=generators(Z);
[n,nrOfGens]=size(G);

if nrOfGens < n || rank(G) < n
    vol = 0;
    return
elseif ~isempty(r) && r > nrOfGens - n
    warning(sprintf(['Number of filtered generators matches/exceeds size of generator matrix!\n',...
        '--> default method ''exact computation'' selected']));
    approx = false; r = [];
end

if ~approx
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

    vol = 2^n*accVol;

else
    % approximative calculation

    %prefilter generators
    h = vecnorm(G'*G,1); % same as vecnorm((G'*G)',1);
    % for-loop below is faster if dim low and nrOfGens high:
    % h = zeros(nrOfGens,1);
    % for i=1:nrOfGens
    %     h(i)=norm(G(:,i)'*G,1);
    % end
    [~,indexFiltered] = sort(h);
    Gfiltered = G(:,indexFiltered((end-n-r+1):end));

    %dimension and new nrOfGenerators
    [n,nrOfGens] = size(Gfiltered);

    %possible combinations of n=dim generators from all generators
    comb = combinator(nrOfGens,n,'c');
    nrOfComb = length(comb(:,1));

    parallelogramVol = zeros(nrOfComb,1);
    for i=1:nrOfComb
        parallelogramVol(i) = abs(det(Gfiltered(:,comb(i,:))));
    end

    vol = 2^n*sum(parallelogramVol);
    
end


%------------- END OF CODE --------------