classdef syntaxTree
% interval class 
%
% Syntax:  
%    object constructor: Obj = interval(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - left limit
%    input2 - right limit
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  polytope

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      19-June-2015
% Last update:  18-November-2015
%               26-January-2016
%               15-July-2017 (NK)
%               04-November-2019 (Zhuoling Li)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    nodes = {};
    operator = [];
    value = [];
    funHan = [];
    id = [];
end
    
methods
    
    % class constructor
    function obj = syntaxTree(value,id,varargin)
        obj.value = value;
        obj.id = id;
        if nargin > 3
           obj.operator = varargin{1};
           obj.funHan = varargin{2};
           obj.nodes = varargin{3};
        end
    end
    
    % parsing methods
    function res = sin(obj)
        fHan = @(x,y) sin_(x,y);
        res = syntaxTree(sin(obj.value),[],'sin',fHan,{obj});
    end
    
    function res = cos(obj)
        fHan = @(x,y) cos_(x,y);
        res = syntaxTree(cos(obj.value),[],'cos',fHan,{obj});
    end
    
    function res = power(obj,exp)
        if ~all(size(obj) == [1,1])
            res = obj;
            if all(size(exp) == [1,1])
               exp = repmat(exp,size(obj)); 
            end
            for i = 1:size(obj,1)
                for j = 1:size(obj,2)
                    fHan = @(x,y) power_(x,exp(i,j),y); 
                    res(i,j) = syntaxTree(obj(i,j).value^exp(i,j),[], ...
                                          'power',fHan,{obj(i,j)});
                end
            end
        else
            fHan = @(x,y) power_(x,exp,y); 
            res = syntaxTree(obj.value^exp,[],'power',fHan,{obj});
        end
    end
    
    function res = mpower(obj,exp)
        fHan = @(x,y) power_(x,exp,y); 
        res = syntaxTree(obj.value^exp,[],'power',fHan,{obj});
    end
    
    function res = tan(obj)
        fHan = @(x,y) tan_(x,y);
        res = syntaxTree(tan(obj.value),[],'tan',fHan,{obj});
    end
    
    function res = exp(obj)
        fHan = @(x,y) exp_(x,y);
        res = syntaxTree(exp(obj.value),[],'exp',fHan,{obj});
    end
    
    function res = log(obj)
        fHan = @(x,y) log_(x,y);
        res = syntaxTree(log(obj.value),[],'log',fHan,{obj});
    end
    
    function res = sqrt(obj)
        fHan = @(x,y) sqrt_(x,y);
        res = syntaxTree(sqrt(obj.value),[],'sqrt',fHan,{obj});
    end
    
    function res = asin(obj)
        fHan = @(x,y) asin_(x,y);
        res = syntaxTree(asin(obj.value),[],'asin',fHan,{obj});
    end
    
    function res = acos(obj)
        fHan = @(x,y) acos_(x,y);
        res = syntaxTree(acos(obj.value),[],'acos',fHan,{obj});
    end
    
    function res = atan(obj)
        fHan = @(x,y) atan_(x,y);
        res = syntaxTree(atan(obj.value),[],'atan',fHan,{obj});
    end
    
    function res = sinh(obj)
        fHan = @(x,y) sinh_(x,y);
        res = syntaxTree(sinh(obj.value),[],'sinh',fHan,{obj});
    end
    
    function res = cosh(obj)
        fHan = @(x,y) cosh_(x,y);
        res = syntaxTree(cosh(obj.value),[],'cosh',fHan,{obj});
    end
    
    function res = tanh(obj)
        fHan = @(x,y) tanh_(x,y);
        res = syntaxTree(tanh(obj.value),[],'tanh',fHan,{obj});
    end
    
    function res = asinh(obj)
        fHan = @(x,y) asinh_(x,y);
        res = syntaxTree(asinh(obj.value),[],'asinh',fHan,{obj});
    end
    
    function res = acosh(obj)
        fHan = @(x,y) acosh_(x,y);
        res = syntaxTree(acosh(obj.value),[],'acosh',fHan,{obj});
    end
    
    function res = atanh(obj)
        fHan = @(x,y) atanh_(x,y);
        res = syntaxTree(atanh(obj.value),[],'atanh',fHan,{obj});
    end
    
    function res = plus(obj1,obj2)
        if isscalar(obj1) && isscalar(obj2)
            fHan = @(x,y,z) plus_(x,y,z);
            if isa(obj1,'syntaxTree')
               if isa(obj2,'syntaxTree')
                  res = syntaxTree(obj1.value + obj2.value,[],'+',fHan,{obj1,obj2});
               else
                  res = syntaxTree(obj1.value + obj2,[],'+',fHan,{obj1,obj2});
               end
            else
                if isa(obj2,'syntaxTree')
                   res = syntaxTree(obj1 + obj2.value,[],'+',fHan,{obj1,obj2}); 
                end
            end
        else
            if isscalar(obj1)
                obj1 = repmat(obj1,size(obj2));
            elseif isscalar(obj2)
                obj2 = repmat(obj2,size(obj1)); 
            end
            if ~all(size(obj1) ==  size(obj2))
                error('Dimensions must agree!');
            end
            for i = 1:size(obj1,1)
                for j = 1:size(obj1,2)
                    res(i,j) = obj1(i,j) + obj2(i,j);
                end
            end
        end
    end
    
    function res = minus(obj1,obj2)
        fHan = @(x,y,z) minus_(x,y,z);
        if isa(obj1,'syntaxTree')
           if isa(obj2,'syntaxTree')
              res = syntaxTree(obj1.value - obj2.value,[],'-',fHan,{obj1,obj2});
           else
              res = syntaxTree(obj1.value - obj2,[],'-',fHan,{obj1,obj2});
           end
        else
            if isa(obj2,'syntaxTree')
               res = syntaxTree(obj1 - obj2.value,[],'-',fHan,{obj1,obj2}); 
            end
        end
    end
    
    function res = uminus(obj)
        fHan = @(x,y,z) uminus_(x,y);
        res = syntaxTree(-obj.value,[],'uminus',fHan,{obj});
    end
    
    function res = times(obj1,obj2)
        if isscalar(obj1) && isscalar(obj2)
            fHan = @(x,y,z) times_(x,y,z);
            if isa(obj1,'syntaxTree')
               if isa(obj2,'syntaxTree')
                  res = syntaxTree(obj1.value * obj2.value,[],'*',fHan,{obj1,obj2});
               else
                  res = syntaxTree(obj1.value * obj2,[],'*',fHan,{obj1,obj2});
               end
            else
                if isa(obj2,'syntaxTree')
                   res = syntaxTree(obj1 * obj2.value,[],'*',fHan,{obj1,obj2}); 
                end
            end
        elseif ~isscalar(obj1) && ~isscalar(obj2)
            if ~all(size(obj1) == size(obj2))
               error('Dimensions must agree!'); 
            end
            for i = 1:size(obj1,1)
                for j = 1:size(obj1,2)
                    res(i,j) = obj1(i,j) .* obj2(i,j);
                end
            end
        else
            res = obj1 * obj2;
        end
    end
    
    function res = mtimes(obj1,obj2)
        if isscalar(obj1) && isscalar(obj2)
            res = obj1 .* obj2;
        elseif isscalar(obj1)
            for i = 1:size(obj2,1)
               for j = 1:size(obj2,2) 
                   res(i,j) = obj1 .* obj2(i,j);
               end
            end
        elseif isscalar(obj2)
            for i = 1:size(obj1,1)
               for j = 1:size(obj1,2) 
                   res(i,j) = obj1(i,j) .* obj2;
               end
            end
        else
            if size(obj1,2) ~= size(obj2,1) 
               error('Dimensions must agree!'); 
            end
            for i = 1:size(obj1,1)
                for j = 1:size(obj2,2)
                    res(i,j) = sum(obj1(i,:).*obj2(:,j)');
                end
            end
        end
    end
    
    function res = prod(obj,varargin)
        n = 1;
        if nargin <= 1
            if size(obj,1) == 1
                n = 2;
            end
        else
            n = varargin{1};
        end
        S.type = '()';
        if n == 1
            S.subs = {1,':'};
        elseif n == 2
            S.subs = {':', 1};
        else
            error('Wrong input');
        end
        res = subsref(obj,S);
        for i = 2:size(obj, n)
            S.subs = {i,':'};
            res = res .* subsref(obj,S);
        end
    end
    
    function res = sum(obj,varargin)
        n = 1;
        if nargin <= 1
            if size(obj,1) == 1
                n = 2;
            end
        else
            n = varargin{1};
        end
        S.type = '()';
        if n == 1
            S.subs = {1,':'};
        elseif n == 2
            S.subs = {':', 1};
        else
            error('Wrong input');
        end
        res = subsref(obj,S);
        for i = 2:size(obj, n)
            res = res + obj(i);
        end
    end
    
    
    % backpropagation methods
    function res = backpropagation(obj,value,int)
        
        % check if node is base variable or not
        if ~isempty(obj.id)
            
            if ~isIntersecting(value,int(obj.id)) || isempty(value)
                [msg,ID] = errEmptySet();
                error(ID,msg);              
            else
                res = int;
                res(obj.id) = value & int(obj.id);
            end
            
        else
            
            % binary vs. unary operators
            if length(obj.nodes) == 1
                
               % tighten interval 
               val_ = obj.funHan(value,obj.nodes{1}.value);
               
               % backpropagation for the node
               res = backpropagation(obj.nodes{1},val_,int);
               
            else
                
                if isa(obj.nodes{1},'syntaxTree')
                    if isa(obj.nodes{2},'syntaxTree')
                        
                       % tighten interval
                       [val1_,val2_] = obj.funHan(value,obj.nodes{1}.value,obj.nodes{2}.value); 

                       % backpropagation for the node
                       res1 = backpropagation(obj.nodes{1},val1_,int);
                       res2 = backpropagation(obj.nodes{2},val2_,int);

                       res = res1 & res2;
                       
                    else
                        
                        % tighten interval
                       val_ = obj.funHan(value,obj.nodes{1}.value,obj.nodes{2}); 

                       % backpropagation for the node
                       res = backpropagation(obj.nodes{1},val_,int); 
                    end
                    
                else
                        
                   % tighten interval
                   [~,val_] = obj.funHan(value,obj.nodes{1},obj.nodes{2}.value); 

                   % backpropagation for the node
                   res = backpropagation(obj.nodes{2},val_,int);
                end
            end
            
        end
    end
end
end



% Auxiliary Functions -----------------------------------------------------

function res = sin_(int,intPrev)
        
    % intersection with valid domain
    int = int & interval(-1,1);

    % apply inverse function
    resTemp1 = asin(int);  
    resTemp2 = -resTemp1 + pi;

    % compute updated upper bound
    ind1 = ceil((supremum(intPrev)-supremum(resTemp1))/(2*pi));    
    ind2 = ceil((supremum(intPrev)-supremum(resTemp2))/(2*pi));

    if in(resTemp1 + ind1*2*pi,supremum(intPrev)) || ...
       in(resTemp2 + ind2*2*pi,supremum(intPrev))
   
        supNew = supremum(intPrev);
        
    else
        
        diff1 = supremum(intPrev) - (supremum(resTemp1) + (ind1-1)*2*pi);
        diff2 = supremum(intPrev) - (supremum(resTemp2) + (ind2-1)*2*pi);
        
        if diff1 < diff2
           supNew = supremum(resTemp1) + (ind1-1)*2*pi;
        else
           supNew = supremum(resTemp2) + (ind2-1)*2*pi;
        end        
    end
    
    % compute updated lower bound
    ind1 = floor((infimum(intPrev)-infimum(resTemp1))/(2*pi));
    ind2 = floor((infimum(intPrev)-infimum(resTemp2))/(2*pi));

    if in(resTemp1 + ind1*2*pi,infimum(intPrev)) || ...
       in(resTemp2 + ind2*2*pi,infimum(intPrev))
   
        infiNew = infimum(intPrev);
        
    else
        
        diff1 = infimum(resTemp1) + (ind1+1)*2*pi - infimum(intPrev);
        diff2 = infimum(resTemp2) + (ind2+1)*2*pi - infimum(intPrev);
        
        if diff1 < diff2
           infiNew = infimum(resTemp1) + (ind1+1)*2*pi;
        else
           infiNew = infimum(resTemp2) + (ind2+1)*2*pi;
        end        
    end
    
    % intersect with previous value
    res = intPrev & interval(infiNew,supNew);
end

function res = cos_(int,intPrev)

    % intersection with valid domain
    int = int & interval(-1,1);

    % apply intervese function
    resTemp1 = acos(int);
    resTemp2 = -resTemp1;
   
    % compute updated upper bound
    ind1 = ceil((supremum(intPrev)-supremum(resTemp1))/(2*pi));    
    ind2 = ceil((supremum(intPrev)-supremum(resTemp2))/(2*pi));

    if in(resTemp1 + ind1*2*pi,supremum(intPrev)) || ...
       in(resTemp2 + ind2*2*pi,supremum(intPrev))
   
        supNew = supremum(intPrev);
        
    else
        
        diff1 = supremum(intPrev) - (supremum(resTemp1) + (ind1-1)*2*pi);
        diff2 = supremum(intPrev) - (supremum(resTemp2) + (ind2-1)*2*pi);
        
        if diff1 < diff2
           supNew = supremum(resTemp1) + (ind1-1)*2*pi;
        else
           supNew = supremum(resTemp2) + (ind2-1)*2*pi;
        end        
    end
    
    % compute updated lower bound
    ind1 = floor((infimum(intPrev)-infimum(resTemp1))/(2*pi));
    ind2 = floor((infimum(intPrev)-infimum(resTemp2))/(2*pi));

    if in(resTemp1 + ind1*2*pi,infimum(intPrev)) || ...
       in(resTemp2 + ind2*2*pi,infimum(intPrev))
   
        infiNew = infimum(intPrev);
        
    else
        
        diff1 = infimum(resTemp1) + (ind1+1)*2*pi - infimum(intPrev);
        diff2 = infimum(resTemp2) + (ind2+1)*2*pi - infimum(intPrev);
        
        if diff1 < diff2
           infiNew = infimum(resTemp1) + (ind1+1)*2*pi;
        else
           infiNew = infimum(resTemp2) + (ind2+1)*2*pi;
        end        
    end
    
    % intersect with previous value
    res = intPrev & interval(infiNew,supNew);

end

function res = power_(int,exp,intPrev)
       
    if exp == 0
        res = intPrev;
        return;

    % even vs. uneven exponent
    elseif mod(exp,2) == 0

        % intersection with valid domain
        if infimum(int) < 0
            if supremum(int) < 0
               error('Value outside the valid domain!'); 
            else
               int = interval(0,supremum(int));
            end
        end

        % apply inverse function
        temp1 = nthroot(infimum(int),exp);
        temp2 = nthroot(supremum(int),exp);
        
        % select correct solution
        if supremum(intPrev) > 0 && infimum(intPrev) < 0
           m = max(temp1,temp2);
           res = interval(-m,m);
        else
           if supremum(intPrev) > 0
               res = interval(min(temp1,temp2),max(temp1,temp2));
           else
               res = interval(-max(temp1,temp2),-min(temp1,temp2));
           end
        end

    else

        % apply inverse function
        temp1 = nthroot(infimum(int),exp);
        temp2 = nthroot(supremum(int),exp);

        res = interval(min(temp1,temp2),max(temp1,temp2));
    end    

    % intersect with previous value
    res = intPrev & res;
end

function [res1,res2] = plus_(int,intPrev1,intPrev2)

    if isa(intPrev1,'interval')
        if isa(intPrev2,'interval')
            res1 = intPrev1 & (int - intPrev2);
            res2 = intPrev2 & (int - intPrev1);
        else
            res1 = intPrev1 & (int - intPrev2);
            res2 = [];
        end
    else
        res1 = [];
        res2 = intPrev2 & (int - intPrev1);
    end
end

function [res1,res2] = minus_(int,intPrev1,intPrev2)

    if isa(intPrev1,'interval')
        if isa(intPrev2,'interval')
            res1 = intPrev1 & (int + intPrev2);
            res2 = intPrev2 & (intPrev1-int);
        else
            res1 = intPrev1 & (int + intPrev2);
            res2 = [];
        end
    else
        res1 = [];
        res2 = intPrev2 & (intPrev1-int);
    end
end

function res = uminus_(int,intPrev)
    res = -int & intPrev;
end

function [res1,res2] = times_(int,intPrev1,intPrev2)

    if isa(intPrev1,'interval')
        if isa(intPrev2,'interval')
            if ~in(intPrev2,0)
                res1 = intPrev1 & (int/intPrev2);
            else
                res1 = intPrev1;
            end
            if ~in(intPrev1,0)
                res2 = intPrev2 & (int/intPrev1);
            else
                res2 = intPrev2; 
            end
        else
            res1 = intPrev1 & (int/intPrev2);
            res2 = [];
        end
    else
        res1 = [];
        if intPrev1 == 0
            res2 = intPrev2;
        else
            res2 = intPrev2 & (int/intPrev1);
        end
    end
end

%------- From here updated by Zhuoling Li -----------
function res = tan_(int,intPrev)

    % intersection with valid domain
    int = int & interval(-1,1);

    % apply intervese function
    resTemp1 = atan(int);
    resTemp2 = resTemp1 + pi;
   
    % compute updated upper bound
    ind1 = ceil((supremum(intPrev)-supremum(resTemp1))/(2*pi));    
    ind2 = ceil((supremum(intPrev)-supremum(resTemp2))/(2*pi));

    if in(resTemp1 + ind1*2*pi,supremum(intPrev)) || ...
       in(resTemp2 + ind2*2*pi,supremum(intPrev))
   
        supNew = supremum(intPrev);
        
    else
        
        diff1 = supremum(intPrev) - (supremum(resTemp1) + (ind1-1)*2*pi);
        diff2 = supremum(intPrev) - (supremum(resTemp2) + (ind2-1)*2*pi);
        
        if diff1 < diff2
           supNew = supremum(resTemp1) + (ind1-1)*2*pi;
        else
           supNew = supremum(resTemp2) + (ind2-1)*2*pi;
        end        
    end
    
    % compute updated lower bound
    ind1 = floor((infimum(intPrev)-infimum(resTemp1))/(2*pi));
    ind2 = floor((infimum(intPrev)-infimum(resTemp2))/(2*pi));

    if in(resTemp1 + ind1*2*pi,infimum(intPrev)) || ...
       in(resTemp2 + ind2*2*pi,infimum(intPrev))
   
        infiNew = infimum(intPrev);
        
    else
        
        diff1 = infimum(resTemp1) + (ind1+1)*2*pi - infimum(intPrev);
        diff2 = infimum(resTemp2) + (ind2+1)*2*pi - infimum(intPrev);
        
        if diff1 < diff2
           infiNew = infimum(resTemp1) + (ind1+1)*2*pi;
        else
           infiNew = infimum(resTemp2) + (ind2+1)*2*pi;
        end        
    end
    
    % intersect with previous value
    res = intPrev & interval(infiNew,supNew);

end

function res = exp_(int,intPrev)
       
        % intersection with valid domain
        if infimum(int) < 0
            if supremum(int) < 0
               error('Value outside the valid domain!'); 
            else
               int = interval(0,supremum(int));
            end
        end

        % apply inverse function
        temp1 = log(infimum(int));
        temp2 = log(supremum(int));

        res = interval(min(temp1,temp2),max(temp1,temp2));
        
    % intersect with previous value
    res = intPrev & res;
end

function res = log_(int,intPrev)

        % apply inverse function
        temp1 = exp(infimum(int));
        temp2 = exp(supremum(int));

        res = interval(min(temp1,temp2),max(temp1,temp2));
        
    % intersect with previous value
    res = intPrev & res;
end

function res = sqrt_(int,intPrev)

        % intersection with valid domain
        if infimum(int) < 0
            if supremum(int) < 0
               error('Value outside the valid domain!'); 
            else
               int = interval(0,supremum(int));
            end
        end

        % apply inverse function
        temp1 = power(infimum(int),2);
        temp2 = power(supremum(int),2);

        res = interval(min(temp1,temp2),max(temp1,temp2));
        
    % intersect with previous value
    res = intPrev & res;
end

function res = asin_(int,intPrev)

    % intersection with valid domain
    if infimum(int) < -pi/2
         if supremum(int) < -pi/2
            error('Value outside the valid domain!'); 
         else
            int = interval(-pi/2,pi/2);
         end
    end
    
    % apply inverse function
    temp1 = sin(infimum(int));
    temp2 = sin(supremum(int));

    res = interval(min(temp1,temp2),max(temp1,temp2));        
    
    % intersect with previous value
    res = intPrev & res;

end

function res = acos_(int,intPrev)

    % intersection with valid domain
        if infimum(int) < 0
            if supremum(int) < 0
               error('Value outside the valid domain!'); 
            else
               int = interval(0,supremum(int));
            end
        end

    % apply inverse function
    temp1 = cos(infimum(int));
    temp2 = cos(supremum(int));

    res = interval(min(temp1,temp2),max(temp1,temp2));
        
    % intersect with previous value
    res = intPrev & res;

end

function res = atan_(int,intPrev)

    % intersection with valid domain
    if infimum(int) < -pi/2
         if supremum(int) < -pi/2
            error('Value outside the valid domain!'); 
         else
            int = interval(-pi/2,pi/2);
         end
    end

    % apply inverse function
    temp1 = tan(infimum(int));
    temp2 = tan(supremum(int));

    res = interval(min(temp1,temp2),max(temp1,temp2));
            
    % intersect with previous value
    res = intPrev & res;

end

function res = sinh_(int,intPrev)

    % intersection with valid domain
    int = int & interval(-1,1);

    % apply intervese function
    temp1 = asinh(infimum(int));
    temp2 = asinh(supremum(int));
    
    res = interval(min(temp1,temp2),max(temp1,temp2));
        
    % intersect with previous value
    res = intPrev & res;

end

function res = cosh_(int,intPrev)

   % intersection with valid domain
        if infimum(int) < 1
            if supremum(int) < 1
               error('Value outside the valid domain!'); 
            else
               int = interval(1,supremum(int));
            end
        end

    % apply intervese function
    temp1 = acosh(infimum(int));
    temp2 = acosh(supremum(int));
   
    if supremum(intPrev) > 0 && infimum(intPrev) < 0
           m = max(temp1,temp2);
           res = interval(-m,m);
      else
           if supremum(intPrev) > 0
               res = interval(min(temp1,temp2),max(temp1,temp2));
           else
               res = interval(-max(temp1,temp2),-min(temp1,temp2));
           end
     end
    
    % intersect with previous value
    res = intPrev & res;

end

function res = tanh_(int,intPrev)

   % intersection with valid domain
    int = int & interval(-1,1);

    % apply intervese function
    temp1 = atanh(infimum(int));
    temp2 = atanh(supremum(int));
   
    res = interval(min(temp1,temp2),max(temp1,temp2));
    
    % intersect with previous value
    res = intPrev & res;

end

function res = asinh_(int,intPrev)
    
    % apply inverse function
    temp1 = sinh(infimum(int));
    temp2 = sinh(supremum(int));

    res = interval(min(temp1,temp2),max(temp1,temp2));        
    
    % intersect with previous value
    res = intPrev & res;

end

function res = acosh_(int,intPrev)
    
    % intersection with valid domain
        if infimum(int) < 0
            if supremum(int) < 0
               error('Value outside the valid domain!'); 
            else
               int = interval(0,supremum(int));
            end
        end
        
    % apply inverse function
    temp1 = cosh(infimum(int));
    temp2 = cosh(supremum(int));

    res = interval(min(temp1,temp2),max(temp1,temp2));
            
    % intersect with previous value
    res = intPrev & res;

end

function res = atanh_(int,intPrev)
        
    % apply inverse function
    temp1 = tanh(infimum(int));
    temp2 = tanh(supremum(int));

    res = interval(min(temp1,temp2),max(temp1,temp2));
        
    % intersect with previous value
    res = intPrev & res;

end