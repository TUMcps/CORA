function res = test_unitvector()
% test_unitvector - unit test function for instantiation of the standard
%    unit vector
%
% Syntax:
%    res = test_unitvector()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       27-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% empty case
resvec(end+1) = isempty(unitvector(0,0));
resvec(end+1) = isempty(unitvector(1,0)); % index is ignored for n=0

% simple cases
resvec(end+1) = isequal(unitvector(1,1),[1]);
resvec(end+1) = isequal(unitvector(1,2),[1;0]);
resvec(end+1) = isequal(unitvector(2,2),[0;1]);
resvec(end+1) = isequal(unitvector(3,4),[0;0;1;0]);
resvec(end+1) = isequal(unitvector(2,5),[0;1;0;0;0]);

% compare with identity
resvec(end+1) = isequal(eye(3),[unitvector(1,3),unitvector(2,3),unitvector(3,3)]);

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
