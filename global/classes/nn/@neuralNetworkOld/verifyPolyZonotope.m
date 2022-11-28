function [res,x] = verifyPolyZonotope(obj,X0,spec,varargin)
% verifyPolyZonotope - automated verifciation for specification on neural
%    networks
%
% Syntax:  
%    [res,x] = verifyPolyZonotope(obj,X0,spec)
%    [res,x] = verifyPolyZonotope(obj,X0,spec,options)
%
% Inputs:
%    obj - object of class neuralNetworkOld
%    X0 - initial set represented as an object of class interval
%    spec - specificaton represented as an object of class specification
%    options - struct containing algorithm settings
%       .nrGen:         maximum number of generators
%       .layerQuad:     number of layers that are evaluated using a 
%                       quadratic approximation
%
% Outputs:
%    res - true if spec. is satisfied, false if not
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

    res = []; x = [];
    
    % parse input arguments
    nrGen = []; layerQuad = 0;
    if nargin > 3
        options = varargin{1};
        if isfield(options,'nrGen')
           nrGen = options.nrGen; 
        end
        if isfield(options,'layerQuad')
           layerQuad = options.layerQuad; 
        end
    end
    
    % split neural network into one part that is evaluated using a
    % quadratic approximation, and one part using a linear approximation
    if layerQuad > 0
       n = layerQuad;
       NNquad = neuralNetworkOld(obj.W(1:n),obj.b(1:n),obj.actFun(1:n));
       NNlin = neuralNetworkOld(obj.W(n+1:end),obj.b(n+1:end), ...
                                                      obj.actFun(n+1:end));
    end
    
    % loop until specification can be verified or falsified
    listX0 = {struct('set',X0,'specs',spec)};
    
    while true
        
        listX0_ = {}; res_ = zeros(length(listX0_),1); 
        x_ = cell(length(listX0_),1);
        
        % loop over all sets that are not verified yet
        for i = 1:length(listX0)
            
            % compute image of the neural network
            pZ = polyZonotope(listX0{i}.set);
            
            if layerQuad > 0
                R = evaluate(NNquad,pZ,'approx','quad',nrGen);
                R = evaluate(NNlin,R,'approx');
            else
                R = evaluate(obj,pZ,'approx','quad',nrGen);
            end
                                                        
            % enclose polynomial zonotope by a zonotope
            R = zonotopeEnclosure(R,dim(X0));
            
            % check if specifications are satisfied
            ind = checkSpec(listX0{i}.specs,R);
            
            if isempty(ind)
               continue; 
            end
            
            % try to falsify specification using the most critical point
            [temp,x_{i}] = ...
                 falsificationCriticalPoint(obj,X0,R,listX0{i}.specs(ind));

            if ~temp
                res_(i) = -1; continue;
            else
                % split the initial set
                list = splitInitialSet(listX0{i},R,ind);
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

function ind = checkSpec(spec,R)
% determine indices of all specifications that are violated

    ind = [];
    
    for i = 1:size(spec,1)
        
       if strcmp(spec(i,1).type,'safeSet') || isa(spec(i,1).set,'halfspace')
           res = check(spec(i,1),R);
       else
           resTemp = isIntersecting(spec(i,1).set,R,'approx');
           if resTemp
              res = check(spec(i,1),R); 
           else
              res = true;
           end
       end

       if ~res
          ind = [ind; i];
       end
    end
end

function list = splitInitialSet(X0,R,ind)
% split the initial set along the most promising dimension

    n = length(X0.set);
    [val,index] = max(sum(R.Z(:,2:n+1).^2,1));
    if val == 0
      [~,index] = max(rad(X0.set)); 
    end
    temp = split(X0.set,index);
    list = {struct('set',temp{1},'specs',X0.specs(ind)); ...
                             struct('set',temp{2},'specs',X0.specs(ind))};
end

function [res,x] = falsificationCriticalPoint(obj,X0,R,spec)
% try to falsify the specification using the most critical point from the
% initial set extracted from the reachable set

    res = true; x = [];
    c = center(X0); G = diag(rad(X0)); n = length(c);
    
    for i = 1:size(spec,1)
        
        poly = spec(i).set;
        
        if isa(poly,'mptPolytope')
            A = poly.P.A; b = poly.P.b;
        elseif isa(poly,'halfspace')
            A = poly.c'; b = poly.d;
        end
    
        if strcmp(spec(i).type,'unsafeSet')

            % determine halfspace constraint that is violated most
            for j = 1:size(A,1)
               [val,~,a] = supportFunc(R,A(j,:)','lower'); 
               if val < b(j)
                   x = c + G*a(1:n);
                   val = evaluate(obj,x);
                   if ~check(spec(i),val)
                      res = false; return; 
                   end
               end
            end

        elseif strcmp(spec(i).type,'safeSet')

            % determine halfspace constraint that is violated most
            for j = 1:size(A,1)
               [val,~,a] = supportFunc(R,A(j,:)'); 
               if val > b(j)
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

function Z = zonotopeEnclosure(pZ,n)
% enclose a polynomial zonotope by a zonotope, where we explicitely keep
% the dependencies between initial and final states

    Z = zonotope(pZ); 
    c = center(Z); G = generators(Z); 
    G_ = zeros(size(G,1),n); index = [];

    for i = 1:n
        ind = find(pZ.id == i);
        if ~isempty(ind)
           ind_ = find(pZ.expMat(ind,:) == 1);
           for j = 1:length(ind_)
              if sum(pZ.expMat(:,ind_(j))) == 1
                  G_(:,i) = pZ.G(:,ind_(j));
                  index = [index; ind_(j)];
              end
           end
        end
    end
    
    Z = zonotope(c,[G_ G(:,setdiff(1:size(G,2),index))]);
end

%------------- END OF CODE --------------