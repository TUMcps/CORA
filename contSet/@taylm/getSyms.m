function res = getSyms( obj )
% getSyms( obj ) - returns a polynomial in a sym form
%
% Syntax:  
%    res = getSyms( obj )
%
% Inputs:
%    obj - a Taylor model
%
% Outputs:
%    res - sym 
%
% Example: 
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Dmitry Grebenyuk
% Written:      02-October-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

	res = arrayfun(@(a) s_getSyms(a), obj, 'UniformOutput', 0);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end

function res = s_getSyms( obj )
    
    % get coefficients
    c = obj.coefficients;
    % get monomials
    degs = obj.monomials(:, 2:end);
    % get var names
    names = obj.names_of_var;
    % transform the var names to syms
    v = sym([names]);
    % make a syms expression
    res = sum(c.*prod(repmat(v,[size(c,1) 1]).^degs,2));
    
end

