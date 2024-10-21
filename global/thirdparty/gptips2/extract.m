function [mainTree,subTree] = extract(index,parentExpr)
%EXTRACT Extract a subtree from an encoded tree expression.
%
%   [MAINTREE,SUBTREE] = EXTRACT(INDEX,PARENTEXPR)
%
%   Input args:
%   PARENTEXPR (the parent string expression)
%   INDEX (the index in PARENTEXPR of the root node of the subtree to be
%   extracted)
%
%   Output args:
%   MAINTREE (PARENTEXPR with the removed subtree replaced by '$')
%   SUBTREE  (the extracted subtree)
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also PICKNODE, GETDEPTH, GETCOMPLEXITY

cnode = parentExpr(index);
iplus = index + 1;
iminus = index - 1;
endpos = numel(parentExpr);

if cnode == 'x'  %extracting an input terminal (x1, x2 etc.)
    
    section = parentExpr(iplus:endpos);
    inp_comma_ind = strfind(section,',');
    inp_brack_ind = strfind(section,')');
    
    %if none found then string must consist of single input
    if isempty(inp_brack_ind) && isempty(inp_comma_ind)
        mainTree = '$';
        subTree = parentExpr;
    else
        inp_ind = sort([inp_brack_ind inp_comma_ind]);
        final_ind = inp_ind(1) + index;
        subTree = parentExpr(index:final_ind-1);
        mainTree = [parentExpr(1:iminus) '$' parentExpr(final_ind:endpos)];
    end
    
    return
    
elseif cnode == '[' %ERC
    
    cl_sbr = strfind(parentExpr(iplus:endpos),']');
    final_ind = cl_sbr(1)+index;
    subTree = parentExpr(index:final_ind);
    mainTree = [parentExpr(1:iminus) '$' parentExpr(final_ind+1:endpos)];
    return
    
elseif cnode=='?' %ERC token
    subTree = cnode;
    mainTree = parentExpr;
    mainTree(index) = '$';
    return
    
else %otherwise extract a tree with a function node as root
    
    %subtree defined when number open brackets=number of closed brackets
    search_seg = parentExpr(index:endpos);
    
    %get indices of open brackets
    op = strfind(search_seg,'(');
    cl = strfind(search_seg,')');
    
    %compare indices to determine point where num_open=num_closed
    tr_op = op(2:end);
    l_tr_op = numel(tr_op);
    
    hibvec = tr_op-cl(1:l_tr_op);
    
    cl_ind = find(hibvec > 0);
    
    if isempty(cl_ind)
        j = cl(numel(op));
    else
        cl_ind = cl_ind(1);
        j = cl(cl_ind);
    end
    
    subTree = search_seg(1:j);
    mainTree = [parentExpr(1:iminus) '$'  parentExpr(j+index:endpos)];
    return
end