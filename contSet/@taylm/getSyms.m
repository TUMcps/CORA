function res = getSyms(tay)
% getSyms - returns a polynomial in a sym form
%
% Syntax:
%    res = getSyms(tay)
%
% Inputs:
%    tay - a Taylor model
%
% Outputs:
%    res - sym
%
% Example:
%    syms x y
%    func = sin(x+y) + exp(-x) + x*y;
%    tay = taylm(func,interval([1;3],[2;4]),4);
%    getSyms(tay)
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Dmitry Grebenyuk
% Written:       02-October-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

	res = arrayfun(@(a) aux_s_getSyms(a),tay,'UniformOutput',false);
    A = cat(1, res{:});
    res = reshape(A, size(res));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_s_getSyms(tay)
    % get coefficients
    c = tay.coefficients;
    % get monomials
    degs = tay.monomials(:, 2:end);
    % get var names
    names = tay.names_of_var;
    % transform the var names to syms
    v = sym([names]);
    % make a syms expression
    res = sum(c.*prod(repmat(v,[size(c,1) 1]).^degs,2));
end

% ------------------------------ END OF CODE ------------------------------
