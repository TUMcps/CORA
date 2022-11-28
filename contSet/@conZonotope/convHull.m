function cZ = convHull(cZ,varargin)
% convHull - computes the convex hull of a constrained zonotope and one or
%    multiple other set representations
%
% Syntax:  
%    cZ = convHull(cZ,S)
%    cZ = convHull(cZ,S1,...,Sm)
%
% Inputs:
%    cZ - conZonotope object
%    S1,...Sm - contSet objects
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
%
%    Z = [4 2 0 0;4 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ2 = conZonotope(Z,A,b);
%
%    res = convHull(cZ1,cZ2);
%
%    figure; hold on;
%    plot(cZ1,[1,2],'FaceColor','r');
%    plot(cZ2,[1,2],'FaceColor','b');
%    plot(res,[1,2],'g','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/enclose

% Author:        Niklas Kochdumper
% Written:       13-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

%------------- BEGIN CODE --------------

    % use different algorithms for the case with only one or two sets
    if nargin == 2

        S = varargin{1};
        if isemptyobject(S)
            cZ = cZ; return
        end
        if isemptyobject(cZ)
            cZ = S; return
        end
        % find a conZonotope object
        if ~isa(cZ,'conZonotope')
            temp = cZ;
            cZ = S;
            S = temp;
        end

        % handle different classes of the second set
        if isa(S, 'conZonotope')

            cZ = convHullconZonotope(cZ,S);

        elseif isa(S,'zonotope') || isa(S,'interval') || ...
               isa(S,'mptPolytope') || isa(S,'zonoBundle') || ...
               isnumeric(S)

            cZ = convHullconZonotope(cZ,conZonotope(S));
            
        elseif isa(S,'polyZonotope') || isa(S,'conPolyZono')
            
            cZ = convHull(polyZonotope(cZ),S);

        else
            % throw error for given arguments
            throw(CORAerror('CORA:noops',cZ,S));
        end

    elseif all(cellfun(@(x) isa(x,'conZonotope'),varargin,'UniformOutput',true))

        % multiple sets only if all are constrained zonotopes
        list = [{cZ},varargin];
        cZ = convHullMany(list);

    else

        throw(CORAerror('CORA:noops',cZ,varargin(:)));

    end
end


% Auxiliary Functions -----------------------------------------------------

function cZ = convHullconZonotope(cZ1,cZ2)
% compute convex hull of two constrained zonotopes

    % obtain object properties
    c1 = cZ1.Z(:,1);
    c2 = cZ2.Z(:,1);

    G1 = cZ1.Z(:,2:end);
    G2 = cZ2.Z(:,2:end);
    
    n = size(G1,1);
    m1 = size(G1,2);
    m2 = size(G2,2);
    
    A1 = cZ1.A;
    A2 = cZ2.A;
    
    b1 = cZ1.b;
    b2 = cZ2.b;
    
    if isempty(A1)
       A1 = zeros(1,m1);
       b1 = 0;
    end
    
    if isempty(A2)
       A2 = zeros(1,m2);
       b2 = 0;
    end
    
    p1 = size(A1,1);
    p2 = size(A2,1);

    % compute center
    c = 0.5 * (c1 + c2);

    % compute generator matrix
    G = [G1 G2 0.5*(c1-c2) zeros(n,3*m1 + 3*m2)];

    % compute constraint matrix
    A1_ = [blkdiag(A1,A2),0.5*[-b1;b2]];
    A2_ = [-2*eye(2*m1), [-2*eye(m1);-2*eye(m1)]];
    A3_ = [-2*eye(2*m2), [-2*eye(m2);-2*eye(m2)]];
    
    A = blkdiag(A1_,A2_,A3_);
    
    A(p1+p2+1:p1+p2+2*m1,1:m1) = [-2*eye(m1);2*eye(m1)];
    A(p1+p2+2*m1+1:end,m1+1:m1+m2) = [2*eye(m2);-2*eye(m2)];
    A(p1+p2+1:end,m1+m2+1) = [ones(2*m1,1);-ones(2*m2,1)];

    % compute constraint offset
    b = [0.5*b1;0.5*b2;3*ones(2*m1+2*m2,1)];
    
    % remove trivial constraints 0 = 0
    ind = find(sum(abs([A,b]),2) > 0);
    A = A(ind,:);
    b = b(ind);

    % construct resulting constraint zonotope object
    cZ = conZonotope([c,G],A,b);
end


function cZ = convHullMany(list)
% compute convex hull of many constrained zonotopes

    % initialize variables
    n = size(cZ1.Z,1);
    
    A_ = cell(length(list),1);
    b_ = [];
    a = [];
    G_ = [];
    c_ = zeros(n,1);
    
    % loop over all sets
    for i = 1:length(list)
        
        % obtain object properties
        A = list{i}.A;
        b = list{i}.b;
        G = list{i}.Z(:,2:end);
        c = list{i}.Z(:,1);
        
        p = size(A,1);
        m = size(G,2);
        
        % construct constraint matrix and vector
        I = eye(m);
        O = zeros(m);
        O_ = zeros(p,m);
        o = ones(m,1);
        
        if ~isempty(A)
            A_{i} = [-0.5*b A O_ O_ O_;o -2*I -2*I O -2*I;o 2*I O -2*I -2*I];
            b_ = [b_;0.5*b;3*ones(2*m,1)];
        else
            A_{i} = [o -2*I -2*I O -2*I;o 2*I O -2*I -2*I];
            b_ = [b_;3*ones(2*m,1)];
        end
        
        a = [a, 0.5, zeros(1,4*m)];
        
        % construct generator matrix
        G_ = [G_,0.5*c,G,zeros(n,3*m)];
        
        % construct center 
        c_ = c_ + 0.5*c;
        
    end
    
    % construct the resulting conZonotope object
    A = [blkdiag(A_{:});a];
    b = [b_;1-0.5*length(list)];
    
    cZ = conZonotope([c_,G_],A,b);  

end

%------------- END OF CODE --------------