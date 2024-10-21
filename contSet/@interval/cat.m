function I = cat(dim,varargin)
% cat - Overloaded 'cat()' operator for intervals, concatenates 
%    given intervals along the given dimension
%
%
% Syntax:
%    I = cat(dim,I1,I2,...)
%
% Inputs:
%    dim - dimension to concatenate along
%    I1,I2,... - interval object
%
% Outputs:
%    I - interval object
%
% Example:
%    I1 = interval([-0.5;0.3]);
%    I2 = interval([2;3]);
%    I = cat(2,I1,I2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cat, interval/cartProd

% Authors:       Tobias Ladner
% Written:       18-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
Is = varargin;
inputArgsCheck({{dim,'att','numeric',{'scalar','positive'}}})
if ~all(cellfun(@(I) isa(I,'interval'), Is,'UniformOutput',true))
    throw(CORAerror('CORA:wrongValue','>= 2','Given objects must be intervals.'))
end

% read inf/sup
infs = cellfun(@(I) I.inf, Is,'UniformOutput',false);
sups = cellfun(@(I) I.sup, Is,'UniformOutput',false);

% call builtin cat
I = interval(cat(dim,infs{:}), cat(dim,sups{:}));

end

% ------------------------------ END OF CODE ------------------------------
