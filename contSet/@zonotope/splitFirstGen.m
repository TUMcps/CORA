function [Znew] = splitFirstGen(Zdummy,Z)
% splitFirstGen - splits first generator, which is in direction of the
% vector field
%
% Syntax:  
%    [Zrem] = firstSplitGen(Zdummy,Z)
%
% Inputs:
%    Z - cell array of zonotope objects
%
% Outputs:
%    Zrem - cell array of remaining zonotope objects
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 09-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%initialize Znew
Znew=[];

%split first generator
for i=1:length(Z)
    %find longest generator
    G=Z{i}.Z(:,2:end);
    for j=1:length(G(1,:))
        h(j)=norm(G(:,j)'*G,1);
    end
    [value,index]=sort(h); 
    
    %split longest generator
    Ztemp = split(Z{i},index(end));
    %Ztemp = split(Z{i},1);
    %write to Znew
    counter=length(Znew);
    Znew{counter+1}=Ztemp{1};
    Znew{counter+2}=Ztemp{2};
end

%------------- END OF CODE --------------