function exprMut = mutate(expr,gp)
%MUTATE Mutate an encoded symbolic tree expression.
%
%   EXPRMUT = MUTATE(EXPR,GP) mutates the encoded symbolic expression EXPR
%   using the parameters in GP to produce expression EXPRXMUT.
%
%   Remarks:
%
%   This function uses the GP.OPERATORS.MUTATION.MUTATE_PAR field to
%   determine probabilistically which type of mutation event to use.
%
%   The types are:
%
%   mutate_type=1 (Default) Ordinary Koza style subtree mutation
%   mutate_type=2 Mutate a non-constant terminal to another non-constant
%   terminal.
%   mutate_type=3 Constant perturbation. Generate small number using
%   Gaussian and add to existing constant.
%   mutate_type=4 Zero constant. Sets the selected constant to zero.
%   mutate_type=5 Randomise constant. Substitutes a random new value for
%   the selected constant.
%   mutate_type=6 Unity constant. Sets the selected constant to 1.
%
%   The probabilities of each mutate_type being selected are stored in the
%   GP.OPERATORS.MUTATION.MUTATE_PAR field which must be a row vector of
%   length 6.
%
%   I.e.
%
%   GP.OPERATORS.MUTATION.MUTATE_PAR = [p_1 p_2 p_3 p_4 p_5 p_6]
%
%              where p_1 = probability of mutate_type 1 being used
%                    p_2=     "         "    "        2  "      "
%
%              and so on.
%
%   Example:
%
%   If GP.OPERATORS.MUTATION.MUTATE_PAR = [0.5 0.25 0.25 0 0] then approx.
%   1/2 of mutation events will be ordinary subtree mutations, 1/4 will be
%   Gaussian perturbations of a constant and 1/4 will be setting a constant
%   to zero.
%
%   Note:
%
%   If a mutation of a certain type cannot be performed (e.g. there are no
%   constants in the tree) then a default subtree mutation is performed
%   instead.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also CROSSOVER

%pick a mutation event type based on the weights in mutate_par
rnum = rand;

while true %loop until a mutation occurs
    
    %ordinary subtree mutation
    if rnum <= gp.operators.mutation.cumsum_mutate_par(1)
        
        %first get the string position of the node
        position = picknode(expr,0,gp);
        
        %remove the logical subtree for the selected node
        %and return the main tree with '$' replacing the extracted subtree
        maintree = extract(position,expr);
        
        %generate a new tree to replace the old
        %Pick a depth between 1 and gp.treedef.max_mutate_depth;
        depth = ceil(rand*gp.treedef.max_mutate_depth);
        newtree = treegen(gp,depth);
        
        %replace the '$' with the new tree
        exprMut = strrep(maintree,'$',newtree);
        
        break;
        
        %substitute terminals
    elseif rnum <= gp.operators.mutation.cumsum_mutate_par(2)
        
        position = picknode(expr,4,gp);
        
        if position == -1
            rnum = 0;
            continue;
        else
            maintree = extract(position,expr); %extract the constant
            exprMut = strrep(maintree,'$',['x' ...
                sprintf('%d',ceil(rand*gp.nodes.inputs.num_inp))]); %replace it with random terminal
            break;
        end
        
        %constant Gaussian perturbation
    elseif rnum <= gp.operators.mutation.cumsum_mutate_par(3)
        
        position = picknode(expr,3,gp);
        
        if position == -1
            rnum = 0;
            continue;
        else
            
            [maintree,subtree] = extract(position,expr);
            
            %process constant
            oldconst = sscanf(subtree(2:end-1),'%f');
            newconst = oldconst + (randn*gp.operators.mutation.gaussian.std_dev);
            
            %use appropriate string formatting for new constant
            if newconst == fix(newconst)
                newtree = ['[' sprintf('%.0f',newconst) ']'];
            else
                newtree = ['[' sprintf(gp.nodes.const.format,newconst) ']'];
            end
            
            %replace '$' with the new tree
            exprMut = strrep(maintree,'$',newtree);
            break
            
        end
        
        %make constant zero
    elseif rnum <= gp.operators.mutation.cumsum_mutate_par(4)
        
        position = picknode(expr,3,gp);
        
        if position == -1
            rnum = 0;
            continue;
        else
            maintree = extract(position,expr); %extract the constant
            exprMut = strrep(maintree,'$','[0]'); %replace it with zero
            break
        end
        
        %randomise constant
    elseif rnum <= gp.operators.mutation.cumsum_mutate_par(5)
        position = picknode(expr,3,gp);
        
        if position == -1
            rnum = 0;
            continue;
        else
            maintree = extract(position,expr); %extract the constant
            
            %generate a constant in the range
            const = rand*(gp.nodes.const.range(2) - gp.nodes.const.range(1)) + gp.nodes.const.range(1);
            
            %convert the constant to square bracketed string form:
            if const == fix(const)
                arg = ['[' sprintf('%.0f',const) ']'];
            else
                arg = ['[' sprintf(gp.nodes.const.format,const) ']'];
            end
            
            exprMut = strrep(maintree,'$',arg); %insert new constant
            break
        end
        
        %set a constant to 1
    elseif rnum <= gp.operators.mutation.cumsum_mutate_par(6)
        
        position = picknode(expr,3,gp);
        
        if position == -1
            rnum = 0;
            continue;
        else
            maintree = extract(position,expr); %extract the constant
            exprMut = strrep(maintree,'$','[1]'); %replace it with one
            break
        end
    end
end