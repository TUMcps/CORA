function [count,out] = numberOfInputs(f,inpArgs)
% numberOfInputs - computes the number of inputs of a function handle
%
% Syntax:  
%    [count,out] = numberOfInputs(f,inpArgs)
%
% Inputs:
%    f - function handle
%    inpArgs - number of input arguments for the function 
%
% Outputs:
%    count - vector storing the length of each input argument
%    out - length of the output argument
%
% Example:
%    f = @(x,u) [x(1)*x(5)^2; sin(x(3)) + u(2)];
%    numberOfInputs(f,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys

% Author:       Niklas Kochdumper
% Written:      03-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    if inpArgs > 26
       error('That many input arguments are not supported!'); 
    end

    % upper bound
    bound = 100000;
    maxVal = 10;
    s = 'abcdefghijklmnopqrstuvwxyz';
    input = cell(inpArgs,1);
    
    while true
       
        found = 0;
        comb = combinator(maxVal,inpArgs);
        
        if size(comb,1) > bound
           error('Could not determine length of input arguments!'); 
        end
        
        for i = 1:size(comb,1)
            for j = 1:inpArgs
                input{j} = sym(s(j),[comb(i,j),1]);
            end

            % pass symbolic inputs to function handle
            try
                output = f(input{:}); 
                found = 1;
                break;
            catch
                continue;
            end
        end
        
        if found
            break; 
        else
            maxVal = 2*maxVal;
        end
    end
    
    
    % get all symbolic variables that are contained in the output
    vars = symvar(output);
    
    % get required length for each input argument
    count = zeros(inpArgs,1);
    
    for i = 1:length(vars)
        
        % split variable into string and number
        temp = char(vars(i));
        str = temp(1);
        num = str2double(temp(2:end));
        
        % update counter
        ind = find(s == str);
        count(ind) = max(count(ind),num);
    end
    
    out = length(output);

%------------- END OF CODE --------------