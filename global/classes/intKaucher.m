classdef intKaucher
% intKaucher - implementation of Kaucher arithmetic [1]. For demonstration
%    purposes, we consider Example 1 in [2]
%
% Syntax:
%    obj = intKaucher(ll,rl)
%
% Inputs:
%    ll - left limit
%    rl - right limit
%
% Outputs:
%    obj - generated object
%
% Example:
%    % function f
%    f = @(x) x^2 - x;
%
%    % compute gradient
%    syms x;
%    df = gradient(f(x));
%    df = matlabFunction(df);
%
%    % compute bounds for gradient
%    I = interval(2,3);
%    c = center(I);
%    gr = df(I);
% 
%    % compute inner-approximation of the range
%    x = intKaucher(3,2);
%    gr = intKaucher(infimum(gr),supremum(gr));
%    
%    res = f(c) + gr*(x - c);
%
% References:
%    [1] E. Kaucher "Interval Analysis in the Extended Interval Space", 
%        Springer 1980 
%    [2] E. Goubault et al. "Forward Inner-Approximated Reachability of 
%        Non-Linear Continuous Systems", HSCC 2017
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Niklas Kochdumper
% Written:       18-October-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    inf = [];
    sup = [];
end
    
methods 
    
    % class constructor
    function obj = intKaucher(ll,rl)
        if any(size(ll) ~= size(rl))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Input arguments must have identical dimensions!'));
        end
        obj.inf = ll;
        obj.sup = rl;
    end
    
    
    % Mathematical Operators ----------------------------------------------
    
    function res = plus(obj1,obj2)
       
        % switch the objects so that the first one is the interval
        if ~isa(obj1,'intKaucher')
            temp = obj1;
            obj1 = obj2;
            obj2 = temp;
        end
        
        % initialize resulting object
        res = obj1;
        
        % compute the sum
        if isa(obj2,'intKaucher')
            res.inf = res.inf + obj2.inf;
            res.sup = res.sup + obj2.sup;
        elseif isnumeric(obj2)
            res.inf = res.inf + obj2;
            res.sup = res.sup + obj2;
        else
            throw(CORAerror('CORA:notSupported','Desired sum not supported.'));
        end
    end
    
    function res = minus(obj1,obj2)
       res = obj1 + (-obj2);
    end
    
    function res = times(obj1,obj2)
       
        if isa(obj1,'intKaucher')
           if isa(obj2,'intKaucher')
               
               res = obj1;
               
               % different cases of inervals
               if aux_inP(obj1) && aux_inDualZ(obj2)
                    res.inf = obj1.inf * obj2.inf;
                    res.sup = obj1.inf * obj2.sup;
               elseif aux_inMinusP(obj1) && aux_inDualZ(obj2)
                    res.inf = obj1.sup * obj2.sup;
                    res.sup = obj1.sup * obj2.inf;
               elseif aux_inZ(obj1) && aux_inDualZ(obj2)
                    res.inf = 0;
                    res.sup = 0;
               else
                    throw(CORAerror('CORA:notSupported',...
                        'Desired multiplication not supported.'));
               end

           elseif isnumeric(obj2)
               if obj2 < 0
                  res = (-obj1)*(-obj2); 
               else
                  res = obj1;
                  res.inf = res.inf * obj2;
                  res.sup = res.sup * obj2;
               end
           else
                throw(CORAerror('CORA:notSupported',...
                    'Desired multiplication not supported.'));
           end
        else
            if isnumeric(obj1)
               res = obj2 * obj1; 
            else
                throw(CORAerror('CORA:notSupported',...
                    'Desired multiplication not supported.'));
            end
        end     
    end
    
    function res = uminus(obj)
       res = obj;
       res.inf = -obj.sup;
       res.sup = -obj.inf;
    end
    
    function res = mtimes(obj1,obj2)
        
        % two scalars
        if isscalar(obj1) && isscalar(obj2)
           res = obj1.*obj2; 
        
        % scalar and matrix
        elseif isscalar(obj1)
           [n1,n2] = size(obj2);
           res = obj2;
           for i = 1:n1
              for j = 1:n2
                  res(i,j) = obj1.*obj2(i,j);
              end
           end
           
        % scalar and matrix 
        elseif isscalar(obj2)
           [n1,n2] = size(obj1);
           res = obj1;
           for i = 1:n1
              for j = 1:n2
                  res(i,j) = obj1(i,j).*obj2;
              end
           end
           
        % matrix and matrix   
        else
            % initialize result
            infRes = zeros(size(obj1,1),size(obj2,2));
            supRes = zeros(size(obj1,1),size(obj2,2));
            
            % get object properties
            inf1 = obj1.inf;
            sup1 = obj1.sup;
            inf2 = obj2.inf;
            sup2 = obj2.sup;
            
            % matrix multiplication
            for i = 1:size(obj1,1)
                for j = 1:size(obj2,2)
                    temp = intKaucher(0,0);
                    for k = 1:size(obj1,2)
                        temp = temp + intKaucher(inf1(i,k),sup1(i,k)).* ...
                                      intKaucher(inf2(k,j),sup2(k,j));
                    end
                    infRes(i,j) = temp.inf;
                    supRes(i,j) = temp.sup;
                end
            end
            
            res = intKaucher(infRes,supRes);
        end
    end
    
    
    % Additional Methods --------------------------------------------------
    
    function res = interval(obj)
       res = interval(obj.inf,obj.sup); 
    end
    
    function res = prop(obj)
        res = intKaucher(min(obj.inf,obj.sup),max(obj.inf,obj.sup));
    end
    
    function res = isProp(obj)
        res = all(obj.inf >= obj.sup); 
    end
    
    function res = isscalar(obj)
        res = isscalar(obj.inf) && isscalar(obj.sup);
    end
    
    function varargout = size(obj, varargin)
        if nargin > 1
            [varargout{1:nargout}] = size(obj.inf,varargin{1});
        else
            [varargout{1:nargout}] = size(obj.inf);
        end
    end
    
    function display(obj)
        
        % determine size of interval
        [rows, columns] = size(obj.inf);

        name = [inputname(1), ' = '];
        disp(name)
        for i = 1:rows
            str = [];
            % display one row
            for j = 1:columns
                newStr = sprintf('[%0.5f,%0.5f]',full(obj.inf(i,j)),full(obj.sup(i,j)));
                str = [str,' ',newStr];
            end
            disp(str);
        end
    end
    
    function res = subsref(obj, S)

        res = obj;

        if strcmp(S.type,'()')
            
            % only one index specified
            if length(S.subs)==1
                res.inf=obj.inf(S.subs{1});
                res.sup=obj.sup(S.subs{1});
                
            % two indices specified
            elseif length(S.subs)==2
                if strcmp(S.subs{1},':')
                    column=S.subs{2};
                    res.inf=obj.inf(:,column);
                    res.sup=obj.sup(:,column);
                elseif strcmp(S.subs{2},':')
                    row=S.subs{1};
                    res.inf=obj.inf(row,:);
                    res.sup=obj.sup(row,:);
                elseif isnumeric(S.subs{1}) && isnumeric(S.subs{1})
                    row=S.subs{1};
                    column=S.subs{2};
                    res.inf=obj.inf(row,column);
                    res.sup=obj.sup(row,column);
                end
            end
        end
    end
    
    function obj = subsasgn(obj, S, value)

        % check if value is an interval
        if ~isa(value,'intKaucher')
            value = interval(value,value);
        end

        if strcmp(S.type,'()')
            
            % only one index specified
            if length(S.subs)==1
                obj.inf(S.subs{1}) = value.inf;
                obj.sup(S.subs{1}) = value.sup;
                
            % two indices specified
            elseif length(S.subs)==2
                if strcmp(S.subs{1},':')
                    column=S.subs{2};
                    obj.inf(:,column) = value.inf;
                    obj.sup(:,column) = value.sup;
                elseif strcmp(S.subs{2},':')
                    row=S.subs{1};
                    obj.inf(row,:) = value.inf;
                    obj.sup(row,:) = value.sup; 
                elseif isnumeric(S.subs{1}) && isnumeric(S.subs{1})
                    row=S.subs{1};
                    column=S.subs{2};
                    obj.inf(row,column) = value.inf;
                    obj.sup(row,column) = value.sup;
                end
            end
        end
    end

end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_inP(obj)
    res = obj.inf >= 0 && obj.sup >= 0;
end

function res = aux_inMinusP(obj)
    res = obj.inf <= 0 && obj.sup <= 0;
end

function res = aux_inZ(obj)
    res = obj.inf <= 0 && obj.sup >= 0;
end

function res = aux_inDualZ(obj)
    res = obj.inf >= 0 && obj.sup <= 0;
end

% ------------------------------ END OF CODE ------------------------------
