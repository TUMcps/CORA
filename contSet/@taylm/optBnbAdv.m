function int = optBnbAdv(obj)
% optBnbAdv - branch and bound optimization with re-expansion of the taylor
%             models at the centers of the subdomains
%
% Syntax:  
%    res = optBnbAdv( obj )
%
% Inputs:
%    obj - taylm
%
% Outputs:
%    res - interval
%
% Example:
%   f= @(x) 1 + x.^5 - x.^4;
%   x = interval(0,1);
%   t = taylm(x,10,'x');
%   T = f(t);
%
%   intReal = interval(0.91808,1)
%
%   int = interval(T)
%   intBnb = interval(T,'bnb')
%   intBnbAdv = interval(T,'bnbAdv')
%   intLinQuad = interval(T,'linQuad')
%
%   X = 0:0.01:1;
%   Y = f(X);
%   plot(X,Y);
%   xlim([0,1]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, taylm/interval, optBnb, optLinQuad
%
% References: 
%   [1] M. Althoff et al. "Implementation of Taylor models in CORA 2018
%       (Tool Presentation)"

% Author:       Niklas Kochdumper
% Written:      14-April-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % Implementation of Alogrithm 1 from reference paper [1]
    % (advanced version, with reexpansion)

    % calculate initial bounds with standarad interval arithmetic
    int = interval(obj,'int');

    % initialize all variables
    unitVec = ones(length(obj.names_of_var),1);
    dom{1} = interval(-unitVec,unitVec);
    bounds{1} = int;
    
    upperBound = inf;
    lowerBound = -inf;

    % loop until the desired precision is achieved
    while true

       % get indices of the subdomains that contain the upper and the lower bounds 
       temp = cell2mat(cellfun(@(x) supremum(x),bounds,'UniformOutput',false));
       [valMax,indMax] = max(temp);
       
       temp = cell2mat(cellfun(@(x) infimum(x),bounds,'UniformOutput',false));
       [valMin,indMin] = min(temp);
       
       % check for convergence 
       if abs(upperBound-valMax) <= 2*obj.eps && abs(lowerBound-valMin) <= 2*obj.eps
          break;
       else
          upperBound = valMax;
          lowerBound = valMin;
       end
       
       % minimum and maximum in the same subdomain
       if indMax == indMin
           
           % split the selected subdomain in half
           [dom1,dom2,int1,int2] = halveDomain(obj,dom,indMax);
           
           % remove the old subdomain from the lists
           bounds{indMax} = [];
           bounds = bounds(~cellfun('isempty',bounds)); 
           
           dom{indMax} = [];
           dom = dom(~cellfun('isempty',dom)); 
           
           % add the splitted subdomains to the lists
           bounds(end+1:end+2) = {int1,int2};
           dom(end+1:end+2) = {dom1,dom2}; 
           
       % minimum and maximum are in different subdomains    
       else
           
           % split the selected subdomains in half
           [dom1,dom2,int1,int2] = halveDomain(obj,dom,indMax);
           [dom3,dom4,int3,int4] = halveDomain(obj,dom,indMin);
           
           % remove the old subdomain from the lists
           bounds{indMax} = [];
           bounds{indMin} = [];
           bounds = bounds(~cellfun('isempty',bounds)); 
           
           dom{indMax} = [];
           dom{indMin} = [];
           dom = dom(~cellfun('isempty',dom)); 
           
           % add the splitted subdomains to the lists
           bounds(end+1:end+4) = {int1,int2,int3,int4};
           dom(end+1:end+4) = {dom1,dom2,dom3,dom4};
                  
       end
    end
    
    % construct the overall boundaries
    int = interval(valMin,valMax);
end


function [dom1,dom2,int1,int2] = halveDomain(obj,dom,ind)

   % determine coordinate with the larges radius
   radius = rad(dom{ind});
   [~,indCoord] = max(radius);

   % split the subdomain along the coordinate with largest radius
   infi = infimum(dom{ind});
   sup = supremum(dom{ind});

   temp = sup;
   temp(indCoord) = infi(indCoord) + radius(indCoord);
   dom1 = interval(infi,temp);

   temp = infi;
   temp(indCoord) = sup(indCoord) - radius(indCoord);
   dom2 = interval(temp,sup);

   % re-expand the taylor models at the center of the halved
   % domains
   t1 = reexpand(obj,dom1);
   t2 = reexpand(obj,dom2);

   % calculate the updated bounds on the new subdomains
   int1 = interval(t1,'int');
   int2 = interval(t2,'int');

end


%------------- END OF CODE --------------