function [res,x] = verify(obj,X0,spec)
% verify - automated verifciation for specification on neural networks
%
% Syntax:  
%    [res,x] = verify(obj,X0,spec)
%
% Inputs:
%    obj - object of class neuralNetworkOld
%    X0 - initial set represented as an object of class interval
%    spec - specification represented as an object of class specification
%
% Outputs:
%    res - true if specification is satisfied, false if not
%    x - counterexample in terms of an initial point violating the specs.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetworkOld

% Author:       Niklas Kochdumper
% Written:      23-November-2021             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = [];

    % try to falsify the specification first to detect easy cases fast
    [temp,x] = fastInitialFalsification(obj,X0,spec);
    if ~temp
       res = false; return; 
    end
    
    % loop until specification can be verified or falsified
    listX0 = {struct('set',X0,'specs',(1:size(spec,1)))};
    
    while true
        
        listX0_ = {}; res_ = zeros(length(listX0_),1); 
        x_ = cell(length(listX0_),1);
        
        % loop over all sets that are not verified yet
        for i = 1:length(listX0)
            
            % compute image of the neural network
            R = evaluate(obj,conZonotope(listX0{i}.set));
            
            % check if specifications are satisfied
            ind = checkSpec(spec(listX0{i}.specs),R);
            
            if isempty(ind)
               continue; 
            end
            
            % try to falsify specification using the most critical point
            ind = listX0{i}.specs(ind);
            [temp,x_{i}] = falsificationCriticalPoint(obj,X0,R,spec(ind));

            if ~temp
                res_(i) = -1; continue;
            else
                % refine interval by excluding points already verified
                list = contractInitialSet(listX0{i}.set,R,spec,ind);
                listX0_ = [listX0_; list];
            end
        end
        
        % update list of open sets
        ind = find(res_ == -1);
        if ~isempty(ind)
            res = false; x = x_{ind(1)}; return;
        else
            if isempty(listX0_)
               res = true; return; 
            else
               listX0 = listX0_; 
            end
        end
    end
end


% Auxiliary Functions -----------------------------------------------------

function [res,x] = fastInitialFalsification(obj,X0,spec)
% fast falsification method to detect easy cases fast

    res = true; x = [];

    % map center and point on the boundary through the network
    Z = zonotope(X0); c = center(Z); G = generators(Z);
    points = [c, c+G, c-G];
    
    output = evaluate(obj,points);
    
    % compute robustness for each output value
    r = robustness(spec,output);
    
    if min(r) < 0
       x = points(:,r < 0); x = x(:,1); res = false; return; 
    end
    
    % approximate robustness values by a linear function
    temp = pinv([ones(size(points,2),1),points'])*r';
    a = temp(2:end);
    
    % determine the point that is expected to violate the specification the
    % most based on the linear estimation
    fac = sign(a);
    x0 = center(X0) + diag(rad(X0))*fac;
    val = evaluate(obj,x0);
    if ~check(spec,val)
       x = x0; res = false; 
    end
end

function ind = checkSpec(spec,R)
% determine indices of all specifications that are violated

    ind = [];
    
    for i = 1:size(spec,1)
       if ~check(spec(i),R)
          ind = [ind; i];
       end
    end
end

function list = contractInitialSet(X0,R,spec,ind)
% contract the initial set by excluding states that are already guaranteed
% to satisfy the specification

    % initialization
    l = infimum(X0); u = supremum(X0); c = (l+u)/2; G = diag((u-l)/2);
    n = length(c); list = []; splitting = false;
    
    % loop over all specifications
    for i = 1:length(ind)
        
       if strcmp(spec(ind(i)).type,'unsafeSet')
           
           % contract interval based on intersection with unsafe set
           temp = R & spec(ind(i)).set;
           cZ = conZonotope(c,[G,zeros(n,size(temp.A,2)-n)],temp.A,temp.b);
           int = interval(cZ);
           
           % check if contracted, otherwise split initial set
           if volume(int(rad(int)>0)) > 0.5*volume(X0(rad(X0)>0))
               splitting = true;
           else
               list = combineSets(list,int,ind(i));
           end
           
       elseif strcmp(spec(ind(i)).type,'safeSet')
           
           poly = spec(ind(i)).set;
           
           % loop over all halfspace constraints
           for j = 1:size(poly.P.A,1)
              
              % check if set is fully in current halfspace or not
              hs = halfspace(-poly.P.A(j,:),-poly.P.b(j));
              
              if isIntersecting(hs,R)
                  
                   % contract interval based on intersection with unsafe 
                   % set = the outside of the safe set
                   temp = R & hs;
                   cZ = conZonotope(c,[G,zeros(n,size(temp.A,2)-n)], ...
                                                            temp.A,temp.b);
                   int = interval(cZ);
                   
                   % check if contracted, otherwise split initial set
                   if volume(int(rad(int)>0)) > 0.5*volume(X0(rad(X0)>0))
                       splitting = true; break;
                   else
                       list = combineSets(list,int,ind(i));
                   end
              end
           end
       end
       
       % split the initial set if it could not be contracted based on the 
       % specifications
       if splitting
           [~,index] = max(sum(R.A(:,1:n).^2,1));
           temp = split(X0,index);
           list = {struct('set',temp{1},'specs',ind); ...
                                        struct('set',temp{2},'specs',ind)};
       end
    end
end

function list = combineSets(list,set,ind)
% combine intervals that are close together to avoid exponential growth of
% the number of sets

    found = false;
    
    for i = 1:length(list)
       if contains(list{i}.set,center(set)) || contains(set,center(list{i}.set))
          list{i}.set = list{i}.set | set;
          list{i}.specs = unique([list{i}.specs, ind]);
          found = true;
       end
    end
    
    if ~found
       list = [list; {struct('set',set,'specs',ind)}]; 
    end
end

function [res,x] = falsificationCriticalPoint(obj,X0,R,spec)
% try to falsify the specification using the most critical point from the
% initial set extracted from the reachable set

    res = true; x = [];
    Z = zonotope(X0); c = center(Z); G = generators(Z); n = length(c);
    
    for i = 1:size(spec,1)
        
        poly = spec(i).set;
    
        if strcmp(spec(i).type,'unsafeSet')

            % determine halfspace constraint that is violated most
            for j = 1:size(poly.P.A,1)
               [val,~,a] = supportFunc(R,poly.P.A(j,:)','lower'); 
               if val < poly.P.b(j)
                   x = c + G*a(1:n);
                   val = evaluate(obj,x);
                   if ~check(spec(i),val)
                      res = false; return; 
                   end
               end
            end

        elseif strcmp(spec(i).type,'safeSet')

            % determine halfspace constraint that is violated most
            for j = 1:size(poly.P.A,1)
               [val,~,a] = supportFunc(R,poly.P.A(j,:)'); 
               if val > poly.P.b(j)
                   x = c + G*a(1:n);
                   val = evaluate(obj,x);
                   if ~check(spec(i),val)
                      res = false; return; 
                   end
               end
            end
        end
    end
end

%------------- END OF CODE --------------