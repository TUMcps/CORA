function obj = removeRedundancies(obj,varargin)
% removeRedundancies - remove redundant halfspaces of a polytope
%
% Syntax:  
%    obj = removeRedundancies(obj)
%    obj = removeRedundancies(obj,method)
%
% Inputs:
%    obj - mptPolytope object
%    method - 'all' (remove all redundant halfspaces) or 'aligned' (remove
%             aligned halfspaces)
%
% Outputs:
%    obj - mptPolytope object
%
% Example:
%    A = [1 1; 1 0; 1 1; -1 0; 0 -1; 1 1; 1 2];
%    b = [7; 3; 5; -1; -1; 6; 20];
%
%    poly = mptPolytope(A,b);
%
%    poly1 = removeRedundancies(poly,'aligned');
%    poly2 = removeRedundancies(poly,'all');
%
%    get(poly1,'A')
%    get(poly2,'A')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mptPolytope

% Author:       Niklas Kochdumper
% Written:      15-April-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    method = setDefaultValues({'all'},varargin);

    % check input arguments
    inputArgsCheck({{obj,'att','mptPolytope'};
                    {method,'str',{'all','aligned'}}});
    
    % remove redundancies
    if strcmp(method,'all')
        
        % remove all redundant halfspaces using the MPT-toolbox mehtod
        P = minHRep(obj.P);
        obj.P = P;
    
    elseif strcmp(method,'aligned')
        
        % remove halfspaces with aligned normal vectors
        obj = removeAligned(obj);
        
    elseif strcmp(method,'linprog')
        
        % remove halfspaces with aligned normal vectors
        obj = removeAligned(obj);
        
        % remove halfspaces using linear programming
        obj = removeLinProg(obj);
        
    end
end


% Auxiliary Functions -----------------------------------------------------

function obj = removeAligned(obj)

    % get object properties
    A = obj.P.A;
    b = obj.P.b;
    
    % normalize the normal vectors
    len = sqrt(sum(A.^2,2));
    temp = diag(1./len);
    
    A = temp*A;
    b = temp*b;
    
    A = A(len > 0,:); 
    b = b(len > 0);
    
    % sort the marix rows to detect aligned normal vectors
    [A,ind] = sortrows(A);
    b = b(ind);
    
    % remove all aligned halfspaces
    counter = 1;
    cTemp = 1;
    A_ = zeros(size(A));
    b_ = zeros(size(b));
    
    while counter < size(A,1)
       
        % check if two normal vectors are identical
        if sum(abs(A(counter,:)-A(counter+1,:))) < eps
            
            a_ = A(counter,:);
            
            % determine minimum b
            bmin = min(b(counter),b(counter+1));
            counter = counter + 2;
            
            while counter <= size(A,1)
                
                if sum(abs(A(counter,:)-a_)) < eps
                   bmin = min(bmin,b(counter));
                   counter = counter + 1;
                else
                   break; 
                end
            end
            
            % store unified normal vector
            A_(cTemp,:) = a_;
            b_(cTemp,:) = bmin;
            cTemp = cTemp + 1;
            
        else
            A_(cTemp,:) = A(counter,:);
            b_(cTemp,:) = b(counter,:);
            cTemp = cTemp + 1;
            counter = counter + 1;
        end  
    end
    
    % add last element
    if counter == size(A,1)
        A_(cTemp,:) = A(end,:);
        b_(cTemp,:) = b(end);
        cTemp = cTemp + 1;
    end
    
    % construct final polytope
    obj = mptPolytope(A_(1:cTemp-1,:),b_(1:cTemp-1));
end

function obj = removeLinProg(obj)
% remove redundant halfspaces using linear programming

    A = obj.P.A; b = obj.P.b;
    ind = []; cnt = 1;
    
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    
    while cnt <= size(A,1)
       
        tmp2 = setdiff(1:size(A,1),cnt);
        try
            [~,val] = linprog(-A(cnt,:),A(tmp2,:),b(tmp2),[],[],[],[],options);
        catch
            val = [];
        end
            
        if isempty(val) || ~(-val < b(cnt) + eps)
           cnt = cnt + 1;
        else
            A(cnt,:) = []; b(cnt) = [];
        end
    end

    obj = mptPolytope(A,b);
end

%------------- END OF CODE --------------