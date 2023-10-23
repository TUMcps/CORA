function Znew = splitFirstGen(Z)
% splitFirstGen - splits first generator, which is in direction of the
%    vector field
%
% Syntax:
%    Znew = splitFirstGen(Z)
%
% Inputs:
%    Z - cell array of zonotope objects
%
% Outputs:
%    Znew - cell array of remaining zonotope objects
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       09-October-2008
% Last update:   14-March-2019 (sort removed)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initialize Znew
Znew=[];

%split first generator
for i=1:length(Z)
    %find longest generator
    G=Z{i}.G;
    for j=1:length(G(1,:))
        h(j)=norm(G(:,j)'*G,1);
    end
    [~,index]=max(h); 
    
    %split longest generator
    Ztemp = split(Z{i},index);
    %Ztemp = split(Z{i},1);
    %write to Znew
    counter=length(Znew);
    Znew{counter+1}=Ztemp{1};
    Znew{counter+2}=Ztemp{2};
end

% ------------------------------ END OF CODE ------------------------------
