function [val,x] = supportFunc(obj,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a zonoBundle object
%               along a certain direction
%
% Syntax:  
%    [val,x] = supportFunc(obj,dir)
%    [val,x] = supportFunc(obj,dir,type)
%
% Inputs:
%    obj - zonoBundle object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%    x - support vector
%
% Example: 
%    zono1 = zonotope([0 1 2 0;0 1 0 2]);
%    zono2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({zono1,zono2});
%
%    val = supportFunc(zB,[1;1]);
%   
%    figure
%    hold on
%    plot(zono1,[1,2],'b');
%    plot(zono2,[1,2],'b');
%    plot(zB,[1,2],'r');
%    plot(conHyperplane(halfspace([1;1],val),[],[]),[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % parse input arguments
    type = 'upper';
    
    if nargin >= 3 && ~isempty(varargin{1})
        type = varargin{1};
    end

    % initialization
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    
    % loop over all parallel sets
    for i = 1:obj.parallelSets
        
       % get object properties 
       zono = obj.Z{i};
       G = generators(zono);
       c = center(zono);
       [n,m] = size(G);
       
       % construct equality constraint matrices
       Aeq = blkdiag(Aeq,-G);
       beq = [beq;c];
       lb = [lb;-ones(m,1)];
       ub = [ub;ones(m,1)];
 
    end
    
    % add optimal point as an additional variable
    A = [eye(size(Aeq,2));-eye(size(Aeq,2))];
    A = [zeros(size(A,1),n),A];
    
    b = [ub;-lb];
    
    Aeq = [repmat(eye(n),[obj.parallelSets,1]),Aeq];
    
    f = [dir;zeros(length(lb),1)];
    
    % linear program options
    options = optimoptions('linprog','display','off');
    
    % upper or lower bound
    if strcmp(type,'lower')
        
       % solve linear program
       [x,val] = linprog(f',A,b,Aeq,beq,[],[],options);
        
    else
        
       % solve linear program
       [x,val] = linprog(-f',A,b,Aeq,beq,[],[],options);
       val = -val;
       
    end
    
    % truncate support vector
    x = x(1:n);
end

%------------- END OF CODE --------------