function T = generateNthTensor(f,vars,order,varargin)
% generateNthTensor - generates the N-th order tensor for function f
%
% Syntax:
%    T = generateNthTensor(f,vars,order)
%    T = generateNthTensor(f,vars,order,Tprev)
%
% Inputs:
%    f - symbolic function
%    vars - symbolic variables of the function
%    order - order of the tensor that is generated
%    Tprev - tensor for order-1 (faster computation if specified)
%
% Outputs:
%    T - resulting symbolic tensor
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: createHigherOrderTensorFiles, evalNthTensor

% Authors:       Niklas Kochdumper
% Written:       08-February-2018
% Last update:   02-February-2021 (MW, different handling of varargin)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
if nargin < 3
    throw(CORAerror('CORA:notEnoughInputArgs',3));
elseif nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs',4));
elseif nargin == 3
    Tprev = []; % previous tensor not provided
elseif nargin == 4
    Tprev = varargin{1};
    % catch user inputs that are not supported
    if order == 1
        throw(CORAerror('CORA:wrongValue','third',...
            'The computation from previous tensor not supported for order = 1.'));
    end
end

    % initialize tensor 
    T = cell(length(f),1);

    % different algorithms depending on whether or not the previous tensor
    % is provided by the user
    if isempty(Tprev)                  % previous tensor not provided
       
        % different initialization depending on whether the tensor order is
        % odd or even 
        if mod(order,2) == 1        % odd tensor order

            % loop over all system dimensions
            for i = 1:length(f)

                first = jacobian(f(i),vars);

                % first order tensor is a special case since the derivative
                % is not stored in a cell array
                if order == 1
                   T{i} = first; 
                else
                   T{i} = cell(length(first),1);
                   for j = 1:length(first)
                      % call of the recursive function
                      T{i}{j} = aux_hessianRecursive(first(j),vars,order-1); 
                   end
                end
            end

        else                        % even tensor order

            % loop over all system dimensions
            for i = 1:length(f)
               % call of the recursive function
               T{i} = aux_hessianRecursive(f(i),vars,order); 
            end
        end 
        
    else                            % previous tensor provided
        
        % use tensor for order-1 to calculate the current tensor
        Tprev = varargin{1};
        
        % different initialization depending on whether the tensor order is
        % odd or even
        if mod(order,2) == 1        % odd tensor order

            % loop over all sysem dimensions
            for i = 1:length(f)

               T{i} = cell(length(vars),1);
               
               for j = 1:length(vars)
                  % call of the recursive function
                  T{i}{j} = aux_hessianFromPrevious(Tprev{i},vars(j)); 
               end
            end

        else                        % even tensor order
            
            % loop over all system dimensions
            for i = 1:length(f)
               
               % second-order tensor is a special case, since it is derived 
               % from the first-order tensor, which is not stored as a cell
               % array
               if order == 2
                    % fill the quadratic matrix with the derivatives of the
                    % first-order tensor
                    T{i} = repmat(Tprev{1}(1),[length(vars),length(vars)]);
                    for k = 1:length(vars)
                       T{i}(k,k) = diff(Tprev{i}(k),vars(k));
                       for j = k+1:length(vars)
                           temp = diff(Tprev{i}(k),vars(j));
                           T{i}(k,j) = temp;
                           T{i}(j,k) = temp;
                       end
                    end
               else
                    % fill the quadratic matrix with derivatives computed
                    % from the previous tensor by the call to the recursive
                    % function
                    T{i} = cell(length(vars));
                    for k = 1:length(vars)
                       T{i}{k,k} = aux_hessianFromPrevious(Tprev{i}{k},vars(k));
                       for j = k+1:length(vars)
                           temp = aux_hessianFromPrevious(Tprev{i}{k},vars(j));
                           T{i}{k,j} = temp;
                           T{i}{j,k} = temp;
                       end
                    end 
               end
            end
        end
    end
    
end


% Auxiliary functions -----------------------------------------------------

function H = aux_hessianRecursive(f,vars,order)
% recursive function that calculates the tensor of the specified order for
% function f

    d = hessian(f,vars);

    if order == 2       % end of recursion
       H = d; 
    else                % next level of the recursion
       H = cell(length(vars));
       for i = 1:length(vars)
           % exploit symmetry in the tensors due to Schwarz's theorem to
           % speed up the computations
           H{i,i} = aux_hessianRecursive(d(i,i),vars,order-2);
           for j = i+1:length(vars)
               temp = aux_hessianRecursive(d(i,j),vars,order-2);
               H{i,j} = temp;
               H{j,i} = temp;
           end
       end
    end
end

function H = aux_hessianFromPrevious(fprev,var)
% recursive function the derivative of the tensor "fprev" with respect to
% the variable var

    if iscell(fprev)    % next level of the recursion
       H = cell(size(fprev));
       for i = 1:size(H,1)
           % exploit symmetry in the tensors due to Schwarz's theorem to
           % speed up the computations
           H{i,i} = aux_hessianFromPrevious(fprev{i,i},var);
           for j = i+1:size(H,1)
               temp = aux_hessianFromPrevious(fprev{i,j},var);
               H{i,j} = temp;
               H{j,i} = temp;
           end
       end
    else            % end of recursion
       H = diff(fprev,var); 
    end
end

% ------------------------------ END OF CODE ------------------------------
