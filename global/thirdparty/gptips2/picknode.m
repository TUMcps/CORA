function position = picknode(expr,nodetype,gp)
%PICKNODE Select a node (or nodes) of specified type from an encoded GP expression and return its position.
%
%   POSITION = PICKNODE(EXPR,NODETYPE,GP) returns the POSITION where EXPR
%   is the symbolic expression, NODETYPE is the specified node type and GP
%   is the GPTIPS data struct.
%
%   Remarks:
%
%   For NODETYPEs 0, 3 or 4 this function returns the string POSITION of
%   the node or -1 if no appropriate node can be found. If NODETYPE is 5 or
%   6 then POSITION is a vector of sorted node indices. If NODETYPE = 7
%   then POSITION is an unsorted vector of node indices.
%
%   Set NODETYPE argument as follows:
%
%   0 = any node
%   1 = unused
%   2 = unused
%   3 = constant (ERC) selection only
%   4 = input selection only
%   5 = indices of nodes with a builtin MATLAB infix representation,
%   sorted left to right (offline use)
%   6 = indices of all nodes, sorted left to right (offline use)
%   7 = gets indices of cube and square functions (offline use)
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
%
%   See also EXTRACT, MUTATE, CROSSOVER

%mask negative constants
x = strrep(expr,'[-','[#');    

if nodetype == 0 %pick any node
    
    %get indices of all function node, constant and input node locations
    %NB char(97) ='a' and char(122) = 'z' and double('[') = 91
    xd = double(x);
    ind = find((xd <= 122 & xd >= 97) | xd==91);
   
elseif nodetype == 3  %just constants
    
    ind = strfind(x,'[');
    
elseif nodetype == 4 %just inputs
    
    ind = strfind(x,'x');
    
elseif nodetype == 5 %special option, only selects plus, minus, times,rdivide,power (intended for offline use)
    
    idcount = 0; ind = [];
    for i=1:numel(gp.nodes.functions.name)
        
        if gp.nodes.functions.active(i)
            idcount = idcount+1;
            func = gp.nodes.functions.name{i};
            
            if strcmpi(func,'times') || strcmpi(func,'plus') ...
                    || strcmpi(func,'minus') || strcmpi(func,'rdivide') || strcmpi(func,'power')
                id = gp.nodes.functions.afid(idcount);
                %look for id in string and concatenate to indices vector
                ind = [ind strfind(x,id)];
            end
            
        end
    end
    position = sort(ind);
    return
    
elseif nodetype == 6 %offline use (all nodes, sorted left to right)
    
    %get indices of function and input locations
    temp_ind = find(double(x)<=122 & double(x)>=97);
    
    %locations of '@','%' and '&' which may have been added by
    %gpreformat.m
    
    %special case of post-run added power node
    temp_ind = [temp_ind strfind(x,'@')];
    
    %special case of post-run added times node
    temp_ind = [temp_ind strfind(x,'&')];
    
    %special case of post-run added exp node
    temp_ind = [temp_ind strfind(x,'%')];
    
    %add indices of constant locations
    ind = [temp_ind strfind(x,'[')];
    
    position = sort(ind);
    return
    
elseif nodetype == 7 %gets all square and cube node indices
    
    idcount = 0;ind = [];
    for i=1:numel(gp.nodes.functions.name)
        
        if gp.nodes.functions.active(i)
            idcount = idcount+1;
            func = gp.nodes.functions.name{i};
            
            if strcmpi(func,'square') || strcmpi(func,'cube')
                id = gp.nodes.functions.afid(idcount);
                %look for id in string and concatenate to indices vector
                ind = [ind strfind(x,id)];
            end
            
        end
    end
    if isempty(ind)
        position = -1;
    else
        position = ind;
    end
    return;
    
end

%count legal nodes
num_nodes = numel(ind);

%if none legal, return -1
if ~num_nodes
    position = -1;
    return;
else
    %select a random node from those indexed
    fun_choice = ceil(rand*num_nodes);
    position = ind(fun_choice);
end