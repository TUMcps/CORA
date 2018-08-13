function [V] = collect(Vdummy,Vcells)
% collect - collects cell arrays of vertices
%
% Syntax:  
%    [V] = collect(V,Vcells)
%
% Inputs:
%    V - vertices object
%    Vcells - cell array of vertices objects
%
% Outputs:
%    V - vertices object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Matthias Althoff
% Written: 07-October-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%delete empty vertices cells
i=1;
while i<=length(Vcells)
    if isempty(Vcells{i}.V)
        Vcells(i)=[];
    else
        i=i+1;
    end
end

%collect vertices
V=Vcells{1};
for i=2:length(Vcells)
   newPoints=length(Vcells{i}.V);
   V.V(:,(end+1):(end+newPoints))=Vcells{i}.V;
end

%------------- END OF CODE --------------