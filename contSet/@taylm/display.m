function display(obj)
% display - Displays the Taylor expression
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    ---
%
% Example: 
%    a = taylm(interval(-1, 1), 3);
%    display(a);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      31-March-2016
%               18-July-2017 (DG) Multivariable polynomial pack is added
% Last update:  11-April-2018 (NK) Sort variables in the polynomial part
% Last revision:---

%------------- BEGIN CODE --------------

name = [inputname(1), ' = '];
disp(name)

% determine the shape of a matrix
[mi,mj] = size(obj);

for i = 1:mi
    rowStr = [];
    for j = 1:mj
        % get a polynomial part; show ratinal numbers as decimals with 5
        % digits
        poly = displayPoly(obj(i,j));
        
        % get an interval part
        remainder = sprintf('[%0.5f,%0.5f]',infimum(obj(i,j).remainder),supremum(obj(i,j).remainder));
        
        % add both parts
        str = [poly, ' + ', remainder];
        
        % add to a row
        rowStr = sprintf('%s \t %s',rowStr ,str);
       
    end
    disp(rowStr)
end

end

function str = displayPoly(obj)

    % get coefficients
    c = obj.coefficients;
    
    % get monomials
    degs = obj.monomials(:,2:end);
    
    % get var names
    names = obj.names_of_var;
    
    % sort var names alphabetically
    [namesSort,ind] = sort(names);
    if size(namesSort,2) ~= 1
       namesSort = transpose(namesSort); 
    end
    
    % adapt exponent matrix to the sorted variable names
    degs = degs(:,ind);
    
    % sort the exponet matrix according to the polynomial order of the
    % terms
    [~,ind] = sortrows([sum(degs,2),fliplr(degs)]);
    degs = degs(ind,:);
    c = c(ind);
    
    % transform the var names to syms
    v = transpose(sym([namesSort]));
    % make a syms expression
    res = c.*prod(repmat(v,[size(c,1) 1]).^degs,2);
    % transform to character array
    temp = num2cell(res);
    charArr = cellfun(@(x) char(vpa(x,5)),temp,'UniformOutput',false);
    
    % join the single monomials
    str = '';
    for i = 1:length(charArr)
       if i == 1
           str = charArr{i};
       else
           if strcmp(charArr{i}(1),'-')
              str = [str,' - ',charArr{i}(2:end)];  
           else
              str = [str,' + ' charArr{i}];
           end
       end
    end


end

%------------- END OF CODE --------------