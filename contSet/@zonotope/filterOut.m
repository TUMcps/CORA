function Zrem = filterOut(Z)
% filterOut - deletes parallelotopes that are covered by other parallelotopes
%
% Syntax:
%    Zrem = filterOut(Z)
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

% Authors:       Matthias Althoff
% Written:       09-October-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initialize Zrem
Zrem=[];

%sort the parallelotopes by volume
for i=1:length(Z)
    vol(i)=volume(Z{i});
end
[val,index]=sort(vol);

%convert to halfspace representation
for i=1:length(Z)
    P{i} = polytope(Z{i});
end

%intersect parallelotopes
for i=1:length(Z)
    ind=index(i);
    Pint=P{ind};
    for j=1:length(Z)
        if j~=ind
            Pint=Pint\P{j};
        end
    end
    %is parallelotope empty?
    [xCheb, RCheb] = chebyball(Pint);
    if RCheb~=-inf
        counter=length(Zrem);
        Zrem{counter+1}=Z{ind};
    else
        disp('canceled!!');
    end        
end


% ------------------------------ END OF CODE ------------------------------
