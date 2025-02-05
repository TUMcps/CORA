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
    
    function IK_plus = plus(IK1,IK2)
       
        % switch the objects so that the first one is the interval
        if ~isa(IK1,'intKaucher')
            IK1_copy = IK1;
            IK1 = IK2;
            IK2 = IK1_copy;
        end
        
        % initialize resulting object
        IK_plus = IK1;
        
        % compute the sum
        if isa(IK2,'intKaucher')
            IK_plus.inf = IK_plus.inf + IK2.inf;
            IK_plus.sup = IK_plus.sup + IK2.sup;
        elseif isnumeric(IK2)
            IK_plus.inf = IK_plus.inf + IK2;
            IK_plus.sup = IK_plus.sup + IK2;
        else
            throw(CORAerror('CORA:notSupported','Desired sum not supported.'));
        end
    end
    
    function IK_minus = minus(IK1,IK2)
        IK_minus = IK1 + (-IK2);
    end
    
    function IK_times = times(IK1,IK2)
       
        if isa(IK1,'intKaucher')
           if isa(IK2,'intKaucher')
               
               IK_times = IK1;
               
               % different cases of inervals
               if aux_inP(IK1) && aux_inDualZ(IK2)
                    IK_times.inf = IK1.inf * IK2.inf;
                    IK_times.sup = IK1.inf * IK2.sup;
               elseif aux_inMinusP(IK1) && aux_inDualZ(IK2)
                    IK_times.inf = IK1.sup * IK2.sup;
                    IK_times.sup = IK1.sup * IK2.inf;
               elseif aux_inZ(IK1) && aux_inDualZ(IK2)
                    IK_times.inf = 0;
                    IK_times.sup = 0;
               else
                    throw(CORAerror('CORA:notSupported',...
                        'Desired multiplication not supported.'));
               end

               % IK2 is numeric
           elseif isnumeric(IK2)
               if IK2 < 0
                  IK_times = (-IK1)*(-IK2); 
               else
                  IK_times = IK1;
                  IK_times.inf = IK_times.inf * IK2;
                  IK_times.sup = IK_times.sup * IK2;
               end
           else
                throw(CORAerror('CORA:notSupported',...
                    'Desired multiplication not supported.'));
           end
        else % IK1 is not intKaucher
            if isnumeric(IK1)
               IK_times = IK2 * IK1; 
            else
                throw(CORAerror('CORA:notSupported',...
                    'Desired multiplication not supported.'));
            end
        end     
    end
    
    function res = uminus(IK)
       res = IK;
       res.inf = -IK.sup;
       res.sup = -IK.inf;
    end
    
    function res = mtimes(factor1,factor2)
        
        % two scalars
        if isscalar(factor1) && isscalar(factor2)
           res = factor1.*factor2; 
        
        % scalar and matrix
        elseif isscalar(factor1)
           [n1,n2] = size(factor2);
           res = factor2;
           for i = 1:n1
              for j = 1:n2
                  res(i,j) = factor1.*factor2(i,j);
              end
           end
           
        % scalar and matrix 
        elseif isscalar(factor2)
           [n1,n2] = size(factor1);
           res = factor1;
           for i = 1:n1
              for j = 1:n2
                  res(i,j) = factor1(i,j).*factor2;
              end
           end
           
        % matrix and matrix   
        else
            % initialize result
            infRes = zeros(size(factor1,1),size(factor2,2));
            supRes = zeros(size(factor1,1),size(factor2,2));
            
            % get object properties
            inf1 = factor1.inf;
            sup1 = factor1.sup;
            inf2 = factor2.inf;
            sup2 = factor2.sup;
            
            % matrix multiplication
            for i = 1:size(factor1,1)
                for j = 1:size(factor2,2)
                    temp = intKaucher(0,0);
                    for k = 1:size(factor1,2)
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
    
    function I = interval(IK)
        I = interval(IK.inf,IK.sup); 
    end
    
    function res = prop(IK)
        % makes a Kaucher interval proper (inf < sup)
        res = intKaucher(min(IK.inf,IK.sup),max(IK.inf,IK.sup));
    end
    
    function res = isProp(IK)
        % checks if a Kaucher interval is proper (inf < sup)
        res = all(IK.inf >= IK.sup); 
    end
    
    function res = isscalar(IK)
        res = isscalar(IK.inf) && isscalar(IK.sup);
    end
    
    function varargout = size(IK, varargin)
        if nargin > 1
            [varargout{1:nargout}] = size(IK.inf,varargin{1});
        else
            [varargout{1:nargout}] = size(IK.inf);
        end
    end
    
    function display(IK)
        
        % determine size of interval
        [rows, columns] = size(IK.inf);

        name = [inputname(1), ' = '];
        disp(name)
        for i = 1:rows
            str = [];
            % display one row
            for j = 1:columns
                newStr = sprintf('[%0.5f,%0.5f]',full(IK.inf(i,j)),full(IK.sup(i,j)));
                str = [str,' ',newStr];
            end
            disp(str);
        end
    end
    
    function IK_sub = subsref(IK, S)

        IK_sub = IK;

        if strcmp(S.type,'()')
            
            % only one index specified
            if length(S.subs)==1
                IK_sub.inf=IK.inf(S.subs{1});
                IK_sub.sup=IK.sup(S.subs{1});
                
            % two indices specified
            elseif length(S.subs)==2
                if strcmp(S.subs{1},':')
                    column=S.subs{2};
                    IK_sub.inf=IK.inf(:,column);
                    IK_sub.sup=IK.sup(:,column);
                elseif strcmp(S.subs{2},':')
                    row=S.subs{1};
                    IK_sub.inf=IK.inf(row,:);
                    IK_sub.sup=IK.sup(row,:);
                elseif isnumeric(S.subs{1}) && isnumeric(S.subs{2})
                    row=S.subs{1};
                    column=S.subs{2};
                    IK_sub.inf=IK.inf(row,column);
                    IK_sub.sup=IK.sup(row,column);
                end
            end
        end
    end
    
    function IK = subsasgn(IK, S, value)

        % check if value is an interval
        if ~isa(value,'intKaucher')
            value = interval(value,value);
        end

        if strcmp(S.type,'()')
            
            % only one index specified
            if length(S.subs)==1
                IK.inf(S.subs{1}) = value.inf;
                IK.sup(S.subs{1}) = value.sup;
                
            % two indices specified
            elseif length(S.subs)==2
                if strcmp(S.subs{1},':')
                    column=S.subs{2};
                    IK.inf(:,column) = value.inf;
                    IK.sup(:,column) = value.sup;
                elseif strcmp(S.subs{2},':')
                    row=S.subs{1};
                    IK.inf(row,:) = value.inf;
                    IK.sup(row,:) = value.sup; 
                elseif isnumeric(S.subs{1}) && isnumeric(S.subs{2})
                    row=S.subs{1};
                    column=S.subs{2};
                    IK.inf(row,column) = value.inf;
                    IK.sup(row,column) = value.sup;
                end
            end
        end
    end

end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_inP(IK)
    res = IK.inf >= 0 && IK.sup >= 0;
end

function res = aux_inMinusP(IK)
    res = IK.inf <= 0 && IK.sup <= 0;
end

function res = aux_inZ(IK)
    res = IK.inf <= 0 && IK.sup >= 0;
end

function res = aux_inDualZ(IK)
    res = IK.inf >= 0 && IK.sup <= 0;
end

% ------------------------------ END OF CODE ------------------------------
