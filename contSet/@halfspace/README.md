# Halfspace

A halfspace is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BHS%7D:=%5C%7Bx%5Cin%5Cmathbb%7BR%7D%5En%5C,%7C%5C,c%5E%5Ctop%20x%5Cleq%20d%5C%7D,%5Cquad%20c%5Cin%5Cmathbb%7BR%7D%5En,d%5Cin%5Cmathbb%7BR%7D."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{HS} := \{ x \in \mathbb{R}^n \, | \, c^\top x \leq d \}, \quad c \in \mathbb{R}^n, d \in \mathbb{R} .
-->

Halfspaces are convex sets.

In CORA, halfspaces are instantiated by

    hs = halfspace(c,d);

where ``c`` is the normal vector of the halfspace and ``d`` is the offset of the normal vector.

Example:

    % initialize halfspace
    hs = halfspace([1 1],1);
    % plot halfspace
    plot(hs);

More information in Section 2.2.2.2 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help halfspace

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>