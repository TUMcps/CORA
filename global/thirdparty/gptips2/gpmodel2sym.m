function [symObj,cellsymObj] = gpmodel2sym(gp,ID,fastMode,useAlias)
%GPMODEL2SYM Create a simplified Symbolic Math object for a multigene symbolic regression model.
%
%   SYMOBJ = GPMODEL2SYM(GP,ID) simplifies and gets the symbolic regression
%   model SYMOBJ with mumeric identifier ID in the GPTIPS datastructure GP.
%
%   SYMOBJ = GPMODEL2SYM(GP,'best') gets the best model of the run (as
%   evaluated on the training data).
%
%   SYMOBJ = GPMODEL2SYM(GP,'valbest') gets the model that performed best
%   on the validation data (if this data exists).
%
%   SYMOBJ = GPMODEL2SYM(GP,'testbest') gets the model that performed best
%   on the test data (if this data exists).
%
%   [SYMOBJ,SYMGENES] = GPMODEL2SYM(GP,ID) gets the overall symbolic model
%   SYMOBJ as well as a cell row array of symbolic objects SYMGENES
%   containing the individual genes (and bias term) that comprise the
%   overall model. The 1st element of SYMGENES is the bias term. The
%   (n+1)th element of SYMGENES is the nth gene of the model.
%
%   [SYMOBJ,SYMGENES] = GPMODEL2SYM(GP,ID,FASTMODE) does the same but uses
%   a slightly faster simplification method (MuPAD symbolic:simplify
%   instead of symbolic:Simplify)
%   
%   [SYMOBJ,SYMGENES] = GPMODEL2SYM(GP,ID,FASTMODE,USEALIAS) does the same
%   but uses aliased variable names (if they exist in
%   gp.nodes.inputs.names) in the returned SYM objects when USEALIAS is
%   TRUE (default = FALSE)
%
%   SYMOBJ = GPMODEL2SYM(GP,GPMODEL) operates on the GPMODEL struct
%   representing a multigene regression model, i.e. the struct returned by
%   the functions GPMODEL2STRUCT or GENES2GPMODEL.
%
%   Remarks:
%
%   As of GPTIPS 2, SYM objects are created using the full default
%   precision of the Symbolic Math toolbox and numbers are represented
%   (except for ERCs) by a ratio of 2 large integers. E.g. consider the
%   following model SYM.
%
%   (5589981434161623*x53)/(8796093022208*x30*(x35 + x46))
%
%   To display SYM objects to N significant digits, use the VPA and CHAR
%   commands, e.g.
%
%   VPA(SYMOBJ,N)
%
%   For instance, using VPA(SYMOBJ,4) on the above example gives:
%
%   (635.5*x53)/(x30*(x35 + x46))
%
%   You can also use PRETTY(VPA(SYMOBJ,4)) to get a more nicely formatted
%   display of the model SYM object (this is similar to the command line
%   output of the GPPRETTY function). This seems glitchy in some versions
%   of MATLAB however.
%
%   Further remarks:
%
%   This assumes each population member is a multigene regression model, 
%   i.e. created using the fitness function REGRESSMULTI_FITFUN - or a 
%   variant thereof with filename beginning REGRESSMULTI.
%
%   (c) Dominic Searson 2009-2015
%
%   GPTIPS 2
%
%   See also GPPRETTY, GPMODEL2MFILE, GPMODELREPORT, GPMODEL2STRUCT,
%   GPMODEL2FUNC

symObj = [];cellsymObj = [];

if ~gp.info.toolbox.symbolic
    error('The Symbolic Math Toolbox is required to use this function.');
end

if nargin < 2
    disp('Basic usage is SYMOBJ = GPMODEL2SYM(GP,ID)');
    disp('or SYMOBJ = GPMODEL2SYM(GP,''best'')');
    disp('or SYMOBJ = GPMODEL2SYM(GP,''valbest'')');
    return;
end

if nargin < 3|| isempty(fastMode)
    fastMode = false;
end

if nargin < 4|| isempty(useAlias)
    useAlias = false;
end

if isnumeric(ID) && ( ID < 1 || ID > numel(gp.pop) )
    error('The supplied numerical model indentifier ID is not valid.');
end

if ~strncmpi('regressmulti', func2str(gp.fitness.fitfun),12)
    error('GPMODEL2SYM may only be used for multigene symbolic regression problems.');
end

try
    [symObj, cellsymObj] = gppretty(gp,ID,[],true,fastMode,useAlias);
catch
    error(['Could not get symbolic object(s) for this model. ' lasterr]);
end