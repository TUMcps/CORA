function res = test_removeCollinearVertices2D()
% test_removeCollinearVertices2D - unit test function for the removal of
%    collinear points
%
% Syntax:
%    res = test_removeCollinearVertices2D()
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

% Authors:       Mark Wetzlinger
% Written:       21-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% init vertices
V = [-1 0; 0 0; 1 0; 1 1]';
% removal
V_ = removeCollinearVertices2D(V);
res(end+1,1) = compareMatrices(V_,[-1 0; 1 0; 1 1]');

% too few vertices
V = [1; 0];
% removal
V_ = removeCollinearVertices2D(V);
res(end+1,1) = compareMatrices(V,V_);

% with tolerance
V = [-1 0; 0 10*eps; 1 0; 1 1]';
% removal
V_ = removeCollinearVertices2D(V,1e-10);
res(end+1,1) = compareMatrices(V_,[-1 0; 1 0; 1 1]',1e-10);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
