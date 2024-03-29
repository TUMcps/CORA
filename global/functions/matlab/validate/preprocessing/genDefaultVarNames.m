function names = genDefaultVarNames(mat,names,inputname)
% genDefaultVarNames - generate default variable names for an taylm-,
%    affine- or zoo-object 
%
% Syntax:
%    names = genDefaultVarNames(mat,names,inputname)
%
% Inputs:
%    mat - matrix containting the values of the variables (size of the
%          matrix determines how many variables are generated)
%    names - variable names that got passed to the function as an input
%            argument
%    inputname - name of the variable that contained the variable values
%   
% Outputs:
%    names - cell array containing the default variable names
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, affine, zoo

% Authors:       Niklas Kochdumper
% Written:       10-April-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    
   % no names provided -> create default names
   if isempty(names)

       text = 'x';
       if ~isempty(inputname)
          text = inputname; 
       end

       % input interval is a scalar
       if isscalar(mat)
           names = text;

       % input interval is a matrix => add indices to name 
       else                  
           names = aux_addMatIndicesToName(text,mat);
       end

   % names provided -> check for correctness
   else
       
       % single string instead of cell-array => add matrix indices
       if ~isscalar(mat) && ~iscell(names)
          names = aux_addMatIndicesToName(names,mat); 
       end

       if iscell(names) && any(size(names) ~= size(mat))
           throw(CORAerror('CORA:wrongValue','second',...
               ' has to be a cell array with the same size as input argument ''int''!'));
       end
   end
end


% Auxiliary functions -----------------------------------------------------

function nameMat = aux_addMatIndicesToName(name,mat)
% add the matrix indices to the names of the variables

   x = 1:size(mat,1);
   y = 1:size(mat,2);
   [X,Y] = meshgrid(x,y);
   Z = 10*X'+Y';

   temp1 = cellfun(@(x) num2str(x),num2cell(Z),'UniformOutput',false);
   temp2 = repmat({name},size(mat));
   nameMat = cellfun(@(x,y) strcat(x,y),temp2,temp1,'UniformOutput',false);

end

% ------------------------------ END OF CODE ------------------------------
