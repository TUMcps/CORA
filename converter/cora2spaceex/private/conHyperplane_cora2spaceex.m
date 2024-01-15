function eqs = conHyperplane_cora2spaceex(obj)
% conHyperplane_cora2spaceex - generate a string that decribes the object
%                              in SpaceEx format
%
% Syntax:
%    eqs = conHyperplane_cora2spaceex(obj)
%
% Inputs:
%    obj - conHyperplane object
%
% Outputs:
%    eqs - string describing the set in SpaceEx format
%
% Example:
%    c = [1 2];
%    d = 3;
%    A = [1 0];
%    b = -2;
%    ch = conHyperplane(c,d,A,b);
%
%    conHyperplane_cora2spaceex(ch)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cora2spaceex

% Authors:       Niklas Kochdumper
% Written:       18-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % hyperplane equation c*x = d
    c      = obj.a';
    d      = obj.b;

    dim = size(c,1);
    x = sym('x',[dim,1]);

    eq = c'*x;

    eqs = [char(eq),' == ',num2str(d)];

    % inequality constraints C*x < d
    if ~isempty(obj.C)
        C = obj.C*x;
        for idx_guard = 1: size(obj.C,1)
            Ax = char(C(idx_guard));
            bx = num2str(obj.d(idx_guard));
            eq_c = [Ax, ' <= ' , bx];
            eqs = [eqs,' & ',newline,eq_c];
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
