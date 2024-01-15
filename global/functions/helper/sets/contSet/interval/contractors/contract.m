function res = contract(f,dom,varargin)
% contract - contracts a interval domain to tightly enclose nonlinear
%            constraints
%
% Syntax:
%    res = contract(f,dom)
%    res = contract(f,dom,alg)
%    res = contract(f,dom,alg,iter)
%    res = contract(f,dom,alg,iter,splits)
%    res = contract(f,dom,alg,iter,splits,jacHan)
%
% Inputs:
%    f - function handle for the constraint f(x) = 0
%    dom - initial domain (class: interval)
%    alg - algorithm used for contraction. The available algorithms are 
%          'forwardBackward' (see Sec. 4.2.4 in [1]),  'linearize' 
%          (see Sec. 4.3.4 in [1]), 'polynomial' (see [2]), 'interval', and 
%          'all' (all contractors together)
%    iter - number of iteration (integer > 0 or 'fixpoint')
%    splits - number of recursive splits (integer > 0)
%    jacHan - handle for Jacobian
%
% Outputs:
%    res - contracted domain (class: interval)
%
% Example: 
%    f = @(x) x(1)^2 + x(2)^2 - 4;
%    dom = interval([1;1],[3;3]);
%   
%    res = contract(f,dom);
%
%    figure; hold on;
%    xlim([-3,3]); ylim([-3,3]);
%    plot(dom,[1,2],'r');
%    plot(res,[1,2],'g');
%    syms x1 x2
%    ls = levelSet(f([x1;x2]),[x1;x2],'==');
%    plot(ls,[1,2],'b');   
%
% References:
%    [1] L. Jaulin et al. "Applied Interval Analysis", 2006
%    [2] G. Trombettoni et al. "A Box-Consistency Contractor Based on 
%        Extremal Functions", 2010
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contractForwardBackward, contractParallelLinearization

% Authors:       Niklas Kochdumper
% Written:       04-November-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % set default values
    [alg,iter,splits,jacHan] = setDefaultValues({'forwardBackward',...
        1,[],[]},varargin);

    inputArgsCheck({{alg,'str',{'forwardBackward','linearize',...
        'polynomial','interval','all'}}});
    % check input arguments
    if ischar(iter)
        inputArgsCheck({{iter,'str','fixpoint'}});
        iter = 10000;
    else
        inputArgsCheck({{iter,'att','numeric',{'integer','nonnegative'}}});
    end
    
    % precompute jacobian matrix
    if isempty(jacHan) && ismember(alg,{'linearize','interval','all'})

        % construct symbolic variables
        n = length(dom);
        vars = sym('x',[n,1]);
        fSym = f(vars);

        % compute jacobian matrix
        jac = jacobian(fSym,vars);
        jacHan = matlabFunction(jac,'Vars',{vars});
    end

    % splitting of intervals considered or not
    if isempty(splits)                                  % no splitting
    
        % iteratively contract the domain
        dom_ = dom;

        for i = 1:iter

            % contract the domain using the selected algorithm
            if strcmp(alg,'forwardBackward')
                dom = contractForwardBackward(f,dom);
            elseif strcmp(alg,'linearize')
                dom = contractParallelLinearization(f,dom,jacHan);
            elseif strcmp(alg,'polynomial')
                dom = contractPolyBoxRevise(f,dom);
            elseif strcmp(alg,'interval')
                dom = contractInterval(f,dom,jacHan);
            elseif strcmp(alg,'all')
                dom = contractForwardBackward(f,dom);
                if ~representsa_(dom,'emptySet',eps)
                    dom = contractInterval(f,dom,jacHan);
                    if ~representsa_(dom,'emptySet',eps)
                        dom = contractParallelLinearization(f,dom,jacHan);
                        if ~representsa_(dom,'emptySet',eps)
                           dom = contractPolyBoxRevise(f,dom); 
                        end
                    end
                end
            end
            
            % check if set is empty
            if representsa_(dom,'emptySet',eps)
               res = [];
               return;
            end

            % check for convergence
            if all(abs(infimum(dom)-infimum(dom_)) < eps) && ...
               all(abs(supremum(dom)-supremum(dom_)) < eps)
                break;
            else
                dom_ = dom;
            end
        end

        res = dom;
        
    else                                                % splitting
       
        % initialization
        list = {dom};
        
        % loop over the number of recursive splits
        for i = 1:splits
            
            list_ = cell(2*length(list),1);
            counter = 1;
            
            % loop over all sets in the list
            for j = 1:length(list)
                
                % determine the best dimension to split and split the 
                % domain along the determined dimension
                domSplit = aux_bestSplit(f,list{j});
                
                % loop over all splitted domains
                for k = 1:length(domSplit)
                    
                    % check if the domain is empty
                    temp = f(domSplit{k});
                    
                    p = zeros(length(temp),1);
                    
                    if ~contains(temp,p)
                       continue; 
                    end
                
                    % contract the splitted domain
                    domTemp = contract(f,domSplit{k},alg,iter,[],jacHan);
                
                    % update the queue
                    if ~representsa_(domTemp,'emptySet',eps)
                        list_{counter} = domTemp;
                        counter = counter + 1;  
                    end
                end
            end
            
            list = list_(1:counter-1);
        end
        
        % unite all the contracted intervals
        if ~isempty(list)
            res = list{1};

            for i = 2:length(list)
               res = res | list{i}; 
            end
        else
            res = [];
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function list = aux_bestSplit(f,dom)
% determine the dimension along which it is best to split the domain

    n = length(dom);
    vol = zeros(n,1);

    % loop over all dimensions
    for i = 1:n
       
        % split the domain at the current dimension
        dom_ = split(dom,i);
        
        % evaluate constraint function for the one of the splitted sets
        val = f(dom_{1});
        
        % evaluate the quality of the split
        vol(i) = volume(val);
    end
    
    % determine best dimension to split
    [~,ind] = min(vol);
    
    % split the domain at the determined dimension
    list = split(dom,ind);
end

% ------------------------------ END OF CODE ------------------------------
