function res = test_fourValued
% test_fourValued - unit test function of fourValued
%
% Syntax:
%    res = test_fourValued
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
% See also: -

% Authors:       Florian Lercher
% Written:       21-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% abbreviations for enum members
tt = fourValued.True;
uu = fourValued.Unknown;
ff = fourValued.False;
ii = fourValued.Inconclusive;

% conjunction
assert(isequal(tt & tt, tt));
assert(isequal(tt & uu, uu));
assert(isequal(tt & ff, ff));
assert(isequal(tt & ii, ii));
assert(isequal(uu & tt, uu));
assert(isequal(uu & uu, uu));
assert(isequal(uu & ff, ff));
assert(isequal(uu & ii, ii));
assert(isequal(ff & tt, ff));
assert(isequal(ff & uu, ff));
assert(isequal(ff & ff, ff));
assert(isequal(ff & ii, ff));
assert(isequal(ii & tt, ii));
assert(isequal(ii & uu, ii));
assert(isequal(ii & ff, ff));
assert(isequal(ii & ii, ii));

% disjunction
assert(isequal(tt | tt, tt));
assert(isequal(tt | uu, tt));
assert(isequal(tt | ff, tt));
assert(isequal(tt | ii, tt));
assert(isequal(uu | tt, tt));
assert(isequal(uu | uu, uu));
assert(isequal(uu | ff, uu));
assert(isequal(uu | ii, ii));
assert(isequal(ff | tt, tt));
assert(isequal(ff | uu, uu));
assert(isequal(ff | ff, ff));
assert(isequal(ff | ii, ii));
assert(isequal(ii | tt, tt));
assert(isequal(ii | uu, ii));
assert(isequal(ii | ff, ii));
assert(isequal(ii | ii, ii));

% negation
assert(isequal(~tt, ff));
assert(isequal(~uu, uu));
assert(isequal(~ff, tt));
assert(isequal(~ii, ii));

% fromBool
assert(isequal(fourValued.fromBool(true), tt));
assert(isequal(fourValued.fromBool(false), ff));

% fromKleene
assert(isequal(fourValued.fromKleene(kleene.True), tt));
assert(isequal(fourValued.fromKleene(kleene.Unknown), uu));
assert(isequal(fourValued.fromKleene(kleene.False), ff));

% arrays
arr1 = [tt,ff,uu,ii];
arr2 = [ff,tt,ii,uu];
assert(isequal(arr1 & arr2, [ff,ff,ii,ii]));
assert(isequal(arr1 | arr2, [tt,tt,ii,ii]));
assert(isequal(~arr1, [ff,tt,uu,ii]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
