function res = test_stlInterval_representsa
% test_stlInterval_representsa - unit test function of representsa
%
% Syntax:
%    res = test_stlInterval_representsa
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

e = stlInterval();
o = stlInterval(0);
p = stlInterval(1);
cc = stlInterval(0,1,true,true);
co = stlInterval(0,1,true,false);
oc = stlInterval(0,1,false,true);
oo = stlInterval(0,1,false,false);

% emptiness
assert(representsa(e,'emptySet'));
assert(~representsa(p,'emptySet'));
assert(~representsa(cc,'emptySet'));
assert(~representsa(co,'emptySet'));
assert(~representsa(oc,'emptySet'));
assert(~representsa(oo,'emptySet'));

% origin
assert(~representsa(e,'origin'));
assert(representsa(o,'origin'));
assert(~representsa(p,'origin'));
assert(~representsa(cc,'origin'));
assert(~representsa(co,'origin'));
assert(~representsa(oc,'origin'));
assert(~representsa(oo,'origin'));

% point
assert(~representsa(e,'point'));
assert(representsa(o,'point'));
assert(representsa(p,'point'));
assert(~representsa(cc,'point'));
assert(~representsa(co,'point'));
assert(~representsa(oc,'point'));
assert(~representsa(oo,'point'));

% interval
assert(representsa(e,'interval'));
assert(representsa(o,'interval'));
assert(representsa(p,'interval'));
assert(representsa(cc,'interval'));
assert(~representsa(co,'interval'));
assert(~representsa(oc,'interval'));
assert(~representsa(oo,'interval'));

% zonotope
assert(representsa(e,'zonotope'));
assert(representsa(o,'zonotope'));
assert(representsa(p,'zonotope'));
assert(representsa(cc,'zonotope'));
assert(~representsa(co,'zonotope'));
assert(~representsa(oc,'zonotope'));
assert(~representsa(oo,'zonotope'));

res = true;

% ------------------------------ END OF CODE ------------------------------
