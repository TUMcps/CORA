function [sendup,status] = pref2inf(expr,gp)
%PREF2INF Recursively extract arguments from a prefix expression and convert to infix where possible.
%
%   [SENDUP,STATUS] = PREF2INF(EXPR,GP)
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
%
%   See also PICKNODE, EXTRACT, GPREFORMAT, TREE2EVALSTR

%get indices of all nodes in current expression
ind = picknode(expr,6,gp);
if isempty(ind)
    sendup = expr; %exits if current expression has no further extractable args
    status = 0;
    return
end

%get first node
ind = ind(1);

%get number of arguments of node
try  
    args = gp.nodes.functions.arity_argt0(strfind(gp.nodes.functions.afid,expr(ind)));
    if isempty(args)
        % if has zero arity then exit
        sendup = expr;
        status = 0;
        return;
    end
catch
    %also exit if error
    sendup = expr;
    status = 0;
    return;
end

if args == 1

    %get subtree rooted in this node
    [~,subtree] = extract(ind,expr);
    node = subtree(1);

    %get index of 2nd node in 'sub_tree' (i.e. first logical argument)
    keynodeInd = picknode(subtree,6,gp);
    keynodeInd = keynodeInd(2);

    %extract the 1st argument expression
    [~,arg1] = extract(keynodeInd,subtree);

    [rec1,rec1_status] = pref2inf(arg1,gp);

    if rec1_status==1
        arg1 = rec1;
    end

    sendup = [node '(' arg1 ')'];
    status = 1;

elseif args > 2
    
    %get subtree rooted in this node
    [~,subtree] = extract(ind,expr);
    node = subtree(1);

    %get index of 2nd node in 'sub_tree' (i.e. first logical argument)
    keynodeInd = picknode(subtree,6,gp); %
    keynodeInd = keynodeInd(2);

    %extract the 1st argument expression
    [remainder,arg1] = extract(keynodeInd,subtree);

    %extract the second argument expression from the remainder
    %find the 1st node after $ in the remainder as this will be the
    %keynode of the 2nd argument we wish to extract
    keynodeInd = picknode(remainder,6,gp);
    tokenInd = strfind(remainder,'$');

    keynodeLoc = find(keynodeInd > tokenInd);
    keynodeInd = keynodeInd(keynodeLoc);
    keynode_ind_1 = keynodeInd(1);
    [remainder2,arg2] = extract(keynode_ind_1,remainder);

    %extract the third argument expression from the remainder
    %find the 1st node after $ in remainder2 as this will be the keynode of
    %the 3nd argument we wish to extract
    keynodeInd = picknode(remainder2,6,gp);
    tokenInd = strfind(remainder2,'$');

    keynodeLoc = find(keynodeInd > max(tokenInd));
    keynodeInd = keynodeInd(keynodeLoc);
    keynode_ind_1 = keynodeInd(1);
    [remainder3,arg3] = extract(keynode_ind_1,remainder2);

    [rec1,rec1_status] = pref2inf(arg1,gp);
    [rec2,rec2_status] = pref2inf(arg2,gp);
    [rec3,rec3_status] = pref2inf(arg3,gp);

    if rec1_status==1
        arg1 = rec1;
    end

    if rec2_status==1
        arg2 = rec2;
    end

    if rec3_status==1
        arg3 = rec3;
    end

    sendup = [node '(' arg1 ',' arg2 ',' arg3 ')'];
    status = 1;

    if args > 3

        %extract the fourth argument expression from the remainder
        %find the 1st node after $ in remainder3 as this will be the
        %keynode of the 4nd argument we wish to extract
        keynodeInd = picknode(remainder3,6,gp);
        tokenInd = strfind(remainder3,'$');

        keynodeLoc = find(keynodeInd > max(tokenInd));
        keynodeInd = keynodeInd(keynodeLoc);
        keynode_ind_1 = keynodeInd(1);
        [~,arg4] = extract(keynode_ind_1,remainder3);
        [rec4,rec4_status] = pref2inf(arg4,gp);

        if rec4_status==1
            arg4 = rec4;
        end

        sendup = [node '(' arg1 ',' arg2 ',' arg3 ',' arg4 ')'];
        status = 1;
    end

else %must have exactly 2 args

    ind = picknode(expr,6,gp);
    %get subtree rooted in this node
    [maintree,subtree] = extract(ind,expr);
    node = subtree(1);

    %get index of 2nd node in 'sub_tree' (i.e. first logical argument)
    keynodeInd = picknode(subtree,6,gp); %
    keynodeInd = keynodeInd(2);

    %extract the 1st argument expression
    [remainder,arg1] = extract(keynodeInd,subtree);

    %extract the second argument expression from the remainder
    %find the 1st node after $ in the remainder as this will be the
    %keynode of the 2nd argument we wish to extract
    keynodeInd = picknode(remainder,6,gp);
    tokenInd = strfind(remainder,'$');

    keynodeLoc = find(keynodeInd>tokenInd);
    keynodeInd = keynodeInd(keynodeLoc);
    keynode_ind_1 = keynodeInd(1);
    [~,arg2] = extract(keynode_ind_1,remainder);

    %process arguments of these arguments if any exist
    [rec1,rec1_status] = pref2inf(arg1,gp);
    [rec2,rec2_status] = pref2inf(arg2,gp);

    if rec1_status==1
        arg1 = rec1;
    end

    if rec2_status==1
        arg2 = rec2;
    end

    %If the node is of infix type (for Matlab symbolic purposes)
    % then send up the results differently
    afid_ind = strfind(gp.nodes.functions.afid, node);
    nodename = gp.nodes.functions.active_name_UC{afid_ind};

    if strcmpi(nodename,'times') || strcmpi(nodename,'minus') ||  ... 
            strcmpi(nodename,'plus') || strcmpi(nodename,'rdivide') || strcmpi(nodename,'power')

        sendup = ['(' arg1 ')' node '(' arg2 ')'];
    else
        sendup = [node '(' arg1 ',' arg2 ')'];
    end
    sendup = strrep(maintree,'$',sendup);

    status = 1; %i.e. ok
end