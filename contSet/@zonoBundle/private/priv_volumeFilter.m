function Gfinal = priv_volumeFilter(G,Z,varargin)
% priv_volumeFilter - filters out generators by directly choosing the
%    smallest volume
%
% Syntax:
%    Gfinal = priv_volumeFilter(G,Z)
%    Gfinal = priv_volumeFilter(G,Z,nrOfPics)
%
% Inputs:
%    G - cells of generator matrices
%    Z - original zonotope
%    nrOfPicks - number of parallelotopes that are picked
%
% Outputs:
%    Gfinal - final generator matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       12-September-2008
% Last update:   15-September-2008
%                25-July-2016 (intervalhull replaced by interval)
%                16-March-2019 (sort removed)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);
nrOfPicks = setDefaultValues({1},varargin);

%obtain dimension
dim=length(G{1}(:,1));
vol = zeros(length(G),1);

%determine generators by exact volume minimization:
for i=1:length(G)
    
    %Get transformation matrix P
    P=G{i};

    %check rank of P
    if rank(P)<dim
        vol(i)=inf;
    else    
  
        %compute reduced zonotope
        Ztrans=pinv(P)*Z;
        Zinterval=interval(Ztrans);
        Zred=P*zonotope(Zinterval);

        %compute volume
        vol(i)=volume(Zred);
    end
end

% obtain indices corresponding to the smallest values
[~,index] = mink(vol, nrOfPicks);

Gfinal = cell(nrOfPicks,1);
for i=1:nrOfPicks
    Gfinal{i}=G{index(i)};
end

% ------------------------------ END OF CODE ------------------------------
