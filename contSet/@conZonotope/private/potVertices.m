function V = potVertices(obj)
% potVertices - calculate all potential vertices of a constrained zonotope
%               object. The points vertices are either real vertices or are
%               located inside the constrained zonotope 
%
% Syntax:  
%    res = potVertices(obj)
%
% Inputs:
%    obj - matrix of size (n,m) containing the m potential vertices
%
% Outputs:
%    res - c-zonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"
%   [2] C. Lara, J. Flores, F. Calderon.
%       "On the Hyperbox - Hyperplane Intersection Problem"

% Author:       Niklas Kochdumper
% Written:      13-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
% Calculate extreme points for the zonotope factors ksi
if isempty(obj.ksi)

    % remove all constraints for which the elimination does not result
    % in an over-approximation (less ksi-dimensions -> speed-up)
    obj = reduceConstraints(obj);
    
    % remove all all-zero columns in the constraint matrix
    temp = sum(abs(obj.A),1);
    ind1 = find(temp == 0);
    ind2 = setdiff(1:size(obj.A,2),ind1);
    
    A_ = obj.A(:,ind2);
    
    % bounding constraints
    n = length(ind2);
    A = [diag(ones(n,1));-diag(ones(n,1))];
    b = ones(2*n,1);
    
    % method for calculation depending on dimension of input
    if size(obj.A,1) ~= 1
        % calculate the extreme points in ksi-space (= polytope vertices)
        try
            ksi_2 = lcon2vert(A,b,A_,obj.b);
        catch
            poly = Polyhedron([A;A_;-A_],[b;obj.b;-obj.b]);
            computeVRep(poly);
            ksi_2 = poly.V;
        end
    elseif n < 5
        % calculate the intersection points via [2] naive method 
        hyperbox = zeros(n,2);
        hyperbox(:,1) = -b(n+1:end,1);
        hyperbox(:,2) = b(1:n,1);
        ksi_2 = boxPlaneIntersectNaive(hyperbox,obj.b,A_);
    else
        % calculate the intersection points via [2] proposed method
        hyperbox = zeros(n,2);
        hyperbox(:,1) = -b(n+1:end,1);
        hyperbox(:,2) = b(1:n,1);
        ksi_2 = boxPlaneIntersect(hyperbox,obj.b,A_);
    end
    
    % combine factors for the potential vertices
    if ~isempty(ind1)
        int = interval(-ones(length(ind1),1),ones(length(ind1),1));
        ksi_1 = vertices(int)';

        ksi = zeros(size(ksi_1,1)*size(ksi_2,1),size(ksi_1,2)+size(ksi_2,2));

        counter = 1;
        m = size(ksi_2,1);

        for i = 1:size(ksi_1,1)
           ksi(counter:counter+m-1,ind1) = ones(m,1) * ksi_1(i,:);
           ksi(counter:counter+m-1,ind2) = ksi_2;
           counter = counter + m;
        end
        
    else
        ksi = ksi_2; 
    end
    
    obj.ksi = ksi';

end
    
% Calculate the corresponding zonotope points (ksi-space -> real space)
V = obj.Z * [ones(1,size(obj.ksi,2)); obj.ksi];

%------------- END OF CODE --------------