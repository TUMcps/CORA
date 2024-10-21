function refexpr = gpreformat(gp,expr,useAlias)
%GPREFORMAT Reformats encoded trees so that the Symbolic Math toolbox can process them properly.
%
%   REFEXPR = GPREFORMAT(GP,EXPR) reformats EXPR which must be an encoded
%   tree expression string - e.g. g(a(m(x2),x1)) - or a linear single
%   dimensional cell array of such encoded string expressions. The output
%   REFEXPR is a decoded expression, e.g. sin((sqrt(x2))*(x1)) or a cell
%   array of such expressions.
%
%   Remarks:
%
%   Genes are accessed as elements of the cell arrays that comprise the
%   population. The jth gene of the ith individual is accessed as
%   GP.POP{I}{J}.
%
%   As of GPTIPS 2, additional processing is performed to allow:
%  (1) better simplification of 'power', 'add3', 'mult3' and 'negexp' nodes
%  (2) processing of any user supplied input variables aliases in the
%   gp.nodes.inputs.names cell array.
%
%   Example:
%
%   To reformat the 1st gene of population member 12 use:
%
%   GPREFORMAT(GP,GP.POP{12}{1})
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPPRETTY, GPMODEL2SYM, GPMODEL2STRUCT, TREE2EVALSTR

if nargin < 2
    disp('Usage is GPREFORMAT(GP,EXPR) where EXPR is a tree expression or a cell array of such expresssions.');
    return;
end

%if useAlias is set to false then variable names will be substituted for
%aliases where they exist
if nargin < 3 || isempty(useAlias)
   useAlias = true; 
end

if ~iscell(expr)
    array_in = {expr};
else
    array_in = expr;
end

len = length(array_in);
refexpr = cell(1,len);

%temporarily add power, times, exp to function set (if not present) to
%allow improved model simplication
gp = addPower(gp);
gp = addTimes(gp);
gp = addexp(gp);

for j=1:len
    
    tempstr = array_in{j};
    
    %skip processing if single input node
    if ~ (tempstr(1) == 'x')
        
        %convert cube,square and sqrt to prefix power
        tempstr = func2power(tempstr,gp);
        
        %convert mult3 to expanded prefix expression using 'times'
        tempstr = mult3_2times(tempstr,gp);
        
        %convert add3 to expanded prefix expression using 'plus'
        tempstr = add3_2plus(tempstr,gp);
        
        %convert 'neg(x)' to 'times(-1,x)
        tempstr = neg2_minus1times(tempstr,gp);
        
        %convert 'negexp(x)' to 'exp(times(-1,x))'
        tempstr = negexp_2exptimes(tempstr,gp);
        
        %Convert certain prefix nodes to infix format
        tempstr = pref2inf(tempstr,gp);
        
        %locate function identifiers of certain infix nodes and replace with
        %+,- etc
        for i=1:length(gp.nodes.functions.afid)
            
            func = gp.nodes.functions.active_name_UC{i};
            id = gp.nodes.functions.afid(i);
            
            if strcmpi(func,'times')
                tempstr = strrep(tempstr,id,'*');
            elseif strcmpi(func,'minus')
                tempstr = strrep(tempstr,id,'-');
            elseif strcmpi(func,'plus')
                tempstr = strrep(tempstr,id,'+');
            elseif strcmpi(func,'rdivide')
                tempstr = strrep(tempstr,id,'/');
            elseif strcmpi(func,'power')
                tempstr = strrep(tempstr,id,'^');
            end
        end
        
        %replace remaining function identifiers with function names
        tempstr = tree2evalstr({tempstr},gp);
        tempstr = tempstr{1};
        
        %remove square brackets
        tempstr = strrep(tempstr,'[','');
        tempstr = strrep(tempstr,']','');
    else %single x node
        tempstr = ['(' tempstr ')'];
    end
    
    %if user defined variable aliases are supplied then replace the x1, x2
    %etc. with the plain versions of the aliases (unless supressed by
    %useAlias function parameter)
    if useAlias && ~isempty(gp.nodes.inputs.names)
        varInds = find(~cellfun(@isempty,gp.nodes.inputs.namesPlain));
        for i = 1:numel(varInds)
            tempstr = strrep(tempstr,['(x' int2str(varInds(i)) ')'],['(' strtrim(gp.nodes.inputs.namesPlain{varInds(i)}) ')']);
        end
    end
    refexpr{j} = tempstr;
end

if ~iscell(expr)
    refexpr = refexpr{1};
end

function gp = addPower(gp)
%check if power is in active function set, if not then add it temporarily
%(marked as afid = @)
powerloc = find(strcmp(gp.nodes.functions.active_name_UC,'POWER'));

if isempty(powerloc)
    gp.nodes.functions.afid = [gp.nodes.functions.afid '@'];
    gp.nodes.functions.active_name_UC{gp.nodes.functions.num_active + 1} = 'POWER';
    gp.nodes.functions.num_active = gp.nodes.functions.num_active + 1;
    gp.nodes.functions.arity_argt0(gp.nodes.functions.num_active) = 2;
    gp.nodes.functions.afid_argt0(gp.nodes.functions.num_active) = '@';
end

function gp = addexp(gp)
%check if 'exp' is in function set, if not then add temporarily (marked as
%adif = %
exploc = find(strcmp(gp.nodes.functions.active_name_UC,'EXP'));

if isempty(exploc)
    gp.nodes.functions.afid=[gp.nodes.functions.afid '%'];
    gp.nodes.functions.active_name_UC{gp.nodes.functions.num_active +1} = 'EXP';
    gp.nodes.functions.num_active = gp.nodes.functions.num_active + 1;
    gp.nodes.functions.arity_argt0(gp.nodes.functions.num_active) = 1;
    gp.nodes.functions.afid_argt0(gp.nodes.functions.num_active)='%';
end

function gp=addTimes(gp)
%check if times is in active function set, if not then add it temporarily
%(marked as afid = &)
timesloc = find(strcmp(gp.nodes.functions.active_name_UC,'TIMES'));

if isempty(timesloc)
    gp.nodes.functions.afid=[gp.nodes.functions.afid '&'];
    gp.nodes.functions.active_name_UC{gp.nodes.functions.num_active+1} = 'TIMES';
    gp.nodes.functions.num_active = gp.nodes.functions.num_active + 1;
    gp.nodes.functions.arity_argt0(gp.nodes.functions.num_active) = 2;
    gp.nodes.functions.afid_argt0(gp.nodes.functions.num_active) = '&';
end

function strOut = mult3_2times(strIn,gp)
%Converts mult3(x1,x2,x3) to times(times(x1,x2),x3) to allow better
%integration & simplification by Symbolic Math toolbox

strOut = strIn;

%locate mult3 and times nodes in active set
mult3loc = find(strcmp(gp.nodes.functions.active_name_UC,'MULT3'));
if isempty(mult3loc)
    return
end

timesloc = find(strcmp(gp.nodes.functions.active_name_UC,'TIMES'));
if isempty(timesloc)
    return
end

%get afids of mult3 and times nodes
mult3node = gp.nodes.functions.afid(mult3loc);
timesnode = gp.nodes.functions.afid(timesloc);

mult3inds = strfind(strIn,mult3node);
num2replace = numel(mult3inds);

%loop until all replaced
for i=1:num2replace;
    
    mult3inds = strfind(strOut,mult3node);
    node = mult3inds(1);
    
    %extract subtree with mult3 as its root
    [main,subtree] = extract(node,strOut);
    
    %get index of the rootnode of the 1st argument of subtree
    argInds = picknode(subtree,6,gp);
    argInd1 = argInds(2); %2nd node is 1st argumnet
    
    %extract the 1st argument of the subtree
    [r1,arg1] = extract(argInd1,subtree);
    
    %get index of the rootnode of the 2nd argument of subtree
    argInds = picknode(r1,6,gp);
    argInd2 = argInds(2);
    
    %extract the 2nd argument of the subtree
    [r2,arg2] = extract(argInd2,r1);
    
    %get index of the rootnode of the 3rd argument of subtree
    argInds = picknode(r2,6,gp);
    argInd3 = argInds(2);
    [~,arg3] = extract(argInd3,r2);
    
    %wrap up the 3 args of mult3 as times(times(arg1,arg2),arg3)
    newtree = [timesnode '(' timesnode '(' arg1 ',' arg2 '),' arg3 ')'];
    strOut = strrep(main,'$',newtree);
    
end

function strOut = neg2_minus1times(strIn,gp)
%converts neg(x1) to times([-1],x) to allow better integration &
%simplification by Symbolic Math toolbox
strOut = strIn;

%locate neg and times in active set
negloc = find(strcmp(gp.nodes.functions.active_name_UC,'NEG'));
if isempty(negloc)
    return
end

timesloc = find(strcmp(gp.nodes.functions.active_name_UC,'TIMES'));
if isempty(timesloc)
    return
end

%get afids of neg and times nodes
negnode = gp.nodes.functions.afid(negloc);
timesnode = gp.nodes.functions.afid(timesloc);

neginds = strfind(strIn,negnode);
num2replace = numel(neginds);

%loop until all replaced
for i=1:num2replace
    
    neginds = strfind(strOut,negnode);
    node = neginds(1);
    
    %extract subtree with neg as its root
    [main,subtree] = extract(node,strOut);
    
    %get index of the rootnode of the 1st argument of subtree
    argInds = picknode(subtree,6,gp);
    argInd1 = argInds(2); %2nd node is 1st argument of neg
    
    %extract the 1st argument of the neg subtree
    [~,arg1] = extract(argInd1,subtree);
    
    %wrap up argument of neg as times([-1],arg1)
    newtree = [timesnode '([-1],' arg1 ')'];
    strOut = strrep(main,'$',newtree);
end

function strOut = negexp_2exptimes(strIn,gp)
%converts negexp(x1) to exp(minus([-1],x)) to allow better integration &
%simplification by Symbolic Math toolbox

strOut = strIn;

%locate negexp, times & exp nodes in active set
negexploc = find(strcmp(gp.nodes.functions.active_name_UC,'NEGEXP'));
if isempty(negexploc)
    return
end

exploc = find(strcmp(gp.nodes.functions.active_name_UC,'EXP'));
if isempty(exploc)
    return
end

timesloc = find(strcmp(gp.nodes.functions.active_name_UC,'TIMES'));
if isempty(timesloc)
    return
end

%get afids of negexp, times and exp nodes
negexpnode = gp.nodes.functions.afid(negexploc);
timesnode = gp.nodes.functions.afid(timesloc);
expnode = gp.nodes.functions.afid(exploc);

negexpinds = strfind(strIn,negexpnode);
num2replace = numel(negexpinds);

%loop until all replaced
for i=1:num2replace
    
    negexpinds = strfind(strOut,negexpnode);
    node = negexpinds(1);
    
    %extract subtree with negexp as its root
    [main,subtree] = extract(node,strOut);
    
    %get index of the rootnode of the 1st argument of subtree
    argInds = picknode(subtree,6,gp);
    argInd1 = argInds(2); %2nd node is 1st argument of negexp
    
    %extract the 1st argument of the negexp subtree
    [~,arg1] = extract(argInd1,subtree);
    
    %wrap up argument of negexp as exp(times([-1],arg1))
    newtree = [expnode '(' timesnode '([-1],' arg1 '))'];
    strOut = strrep(main,'$',newtree);
end

function strOut = add3_2plus(strIn,gp)
%Converts add3(x1,x2,x3) to plus(plus(x1,x2),x3) to allow better
%integration & simplification by Symbolic Math toolbox

strOut = strIn;

%locate add3 and plus nodes in active set
add3loc = find(strcmp(gp.nodes.functions.active_name_UC,'ADD3'));
if isempty(add3loc)
    return
end

plusloc = find(strcmp(gp.nodes.functions.active_name_UC,'PLUS'));
if isempty(plusloc)
    return
end

%get afids of add3 and plus nodes
add3node = gp.nodes.functions.afid(add3loc);
plusnode = gp.nodes.functions.afid(plusloc);

add3inds = strfind(strIn,add3node);
num2replace = numel(add3inds);

%loop until all replaced
for i=1:num2replace;
    
    add3inds = strfind(strOut,add3node);
    node = add3inds(1);
    
    %extract subtree with add3 as its root
    [main,subtree] = extract(node,strOut);
    
    %get index of the rootnode of the 1st argument of subtree
    argInds = picknode(subtree,6,gp);
    argInd1 = argInds(2); %2nd node is 1st argument
    
    %extract the 1st argument of the subtree
    [r1,arg1] = extract(argInd1,subtree);
    
    %get index of the rootnode of the 2nd argument of subtree
    argInds = picknode(r1,6,gp);
    argInd2 = argInds(2);
    
    %extract the 2nd argument of the subtree
    [r2,arg2] = extract(argInd2,r1);
    
    %get index of the rootnode of the 3rd argument of subtree
    argInds = picknode(r2,6,gp);
    argInd3 = argInds(2);
    [~,arg3] = extract(argInd3,r2);
    
    %wrap up the 3 args of add3 as plus(plus(arg1,arg2),arg3)
    newtree = [plusnode '(' plusnode '(' arg1 ',' arg2 '),' arg3 ')'];
    strOut = strrep(main,'$',newtree);
    
end

function strOut = func2power(strIn,gp)
%Converts cube(x) to power(x,3) etc.

%locate power node in active set
powerloc = find(strcmp(gp.nodes.functions.active_name_UC,'POWER'));

if isempty(powerloc)
    return
end

powernode = gp.nodes.functions.afid(powerloc);
inds = picknode(strIn,7,gp);

if inds(1) == -1
    strOut = strIn;
    return;
end
numNodes = numel(inds);

for i=1:numNodes
    inds = picknode(strIn,7,gp);
    ournode = strIn(inds(1));
    ournodeloc = strfind(gp.nodes.functions.afid,ournode);
    funcname = gp.nodes.functions.active_name_UC{ournodeloc};
    
    [main,sub_tree] = extract(inds(1),strIn);
    
    %get index of 2nd node in 'sub_tree' (i.e. first (and only) logical
    %argument)
    keynode_ind = picknode(sub_tree,6,gp); %
    keynode_ind = keynode_ind(2);
    
    %extract the 1st argument expression
    [~,arg1] = extract(keynode_ind,sub_tree);
    
    %wrap this argument into power function
    if strcmpi(funcname,'cube')
        newtree = [powernode '(' arg1 ',[3])'];
    elseif strcmpi(funcname,'square')
        newtree = [powernode '(' arg1 ',[2])'];
    end
    
    strIn = strrep(main,'$',newtree);
end
strOut = strIn;