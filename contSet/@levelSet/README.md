# Level Set

A level set can be defined in multiple ways:

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BLS%7D:=%5C%7Bx%5Cin%5Cmathbb%7BR%7D%5En%5C,%7C%5C,f(x)%5Cleq%200%5C%7D,"/>
</p>
<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BLS%7D:=%5C%7Bx%5Cin%5Cmathbb%7BR%7D%5En%5C,%7C%5C,f(x)%3C0%5C%7D,"/>
</p>
<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BLS%7D:=%5C%7Bx%5Cin%5Cmathbb%7BR%7D%5En%5C,%7C%5C,f(x)=0%5C%7D."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{LS} := \{ x \in \mathbb{R}^n \, | \, f(x) \leq 0 \}
\mathcal{LS} := \{ x \in \mathbb{R}^n \, | \, f(x) < 0 \}
\mathcal{LS} := \{ x \in \mathbb{R}^n \, | \, f(x) = 0 \}
-->

In CORA, level sets are instantiated by

    O = levelSet(eq,vars,op);

where ``eq`` is a symbolic MATLAB function, ``vars`` are the symbolic variables, and ``op`` is the comparison operator.

Example:

    % initialize level set
    vars = sym('x',[2;1]);
    eq = 1/vars(1)^2 - vars(2);
    ls = levelSet(eq,vars,'==');
    % plot level set
    plot(ls);

More information in Section 2.2.2.3 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help levelSet

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>
