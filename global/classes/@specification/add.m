function spec = add(spec1,spec2)
% add - joins two specification objects
%
% Syntax:  
%    spec = add(spec1,spec2)
%
% Inputs:
%    spec1 - specification object
%    spec2 - specification object
%
% Outputs:
%    spec - resulting specification object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Author:       Niklas Kochdumper
% Written:      29-May-2020             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% quick exit in case one of the two specification objects is empty
if isempty(spec1)
    spec = spec2; return;
end
if isempty(spec2)
    spec = spec1; return; 
end

% instantiate joint specification object
spec = repelem(spec1(1,1),size(spec1,1)+size(spec2,1),1);
% counter for total number of specifications
cnt = 1;

% loop over all specifications in the first object
for i = 1:size(spec1,1)
    spec(cnt,1) = spec1(i,1);
    cnt = cnt + 1;
end

% loop over all specifications in the second object
for i = 1:size(spec2,1)
    spec(cnt,1) = spec2(i,1);
    cnt = cnt + 1;
end

%------------- END OF CODE --------------