function res = evalNthTensor(T,x,order)
% evalNthTensor - evaluates the taylor term that corresponds to the
%                 specified tensor
%
% Syntax:  
%    T = evalNthTensor(f,x,order)
%
% Inputs:
%    T - tensor (symbolic or numberic)
%    x - variable values for which the taylor-term is evaluated
%    order - order of the tensor T
%
% Outputs:
%    res - value of the taylor term
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Niklas Kochdumper
% Written:      08-February-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % initialize the variable that stores the resulting values of the
    % taylor term
    res = repmat(x(1),[length(T),1]);

    % different initialization of the algorithm depending on whether the
    % tensor order is odd or even
    if mod(order,2) == 1        % odd tensor order
        
        % loop over all system dimensions
        for i = 1:length(T)
            temp = T{i};
            
            % first-order is a special case, since the derivative there is
            % stored as a matrix instead of a cell array
            if order == 1
                res(i) = temp * x;
            else
                % calculate the value of the term with a recursive function
                res(i) = x(1) * evalQuadratic(temp{1},x);

                for j = 2:length(x)
                    res(i) = res(i) + x(j) * evalQuadratic(temp{j},x);
                end
            end
        end  
        
    else                        % even tensor order
        
        % loop over all system dimensions
        for i = 1:length(T)
            % call of the recursive function
            res(i) = evalQuadratic(T{i},x);
        end
    end 
    
    % multiply by factorial factor to obtain the final result
    res = 1/factorial(order) * res;

end

function res = evalQuadratic(T,x)
% recursive function that evaluates the value of the taylor term that
% corresponds to the tensor T at the point x

    if iscell(T)    % next recursion level
        
        % exploit symmetry in the tensors due to Schwarz's theorem to
        % speed up the computations
        H = repmat(x(1),[length(x),length(x)]);
        for i = 1:length(x)
           H(i,i) = evalQuadratic(T{i,i},x);
           for j = i+1:length(x)
               temp = evalQuadratic(T{i,j},x);
               H(i,j) = temp;
               H(j,i) = temp;
           end
        end
       
        res = transpose(x)* H * x;
        
    else            % end of the recursion
        res = transpose(x)* T * x;
    end
end

%------------- END OF CODE --------------