# Constrained Hyperplane

A constrained hyperplane is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BCH%7D:=%5C%7Bx%5Cin%5Cmathbb%7BR%7D%5En%5C,%7C%5C,c%5E%5Ctop%20x=d,Ax%5Cleq%20b%5C%7D,%5Cquad%20c%5Cin%5Cmathbb%7BR%7D%5En,d%5Cin%5Cmathbb%7BR%7D,A%5Cin%5Cmathbb%7BR%7D%5E%7Bm%5Ctimes%20n%7D,b%5Cin%5Cmathbb%7BR%7D%5Em."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{CH} := \{ x \in \mathbb{R}^n \, | \, c^\top x = d, Ax \leq b \}, \quad c \in \mathbb{R}^n, d \in \mathbb{R}, A \in \mathbb{R}^{m \times n}, b \in \mathbb{R}^m .
-->

Constrained hyperplane are closed, convex, and degenerate sets.

In CORA, constrained hyperplanes are instantiated by

    hyp = conHyperplane(c,d,A,b);

where ``c`` is the normal vector of the hyperplane, ``d`` is the offset of the normal vector, ``A`` is the constraint matrix, and ``b`` is the constraint offset.

Example:

    % initialize constrained hyperplane
    hyp = conHyperplane([1 1],1,[0 1],1);
    % plot constrained hyperplane
    plot(hyp);

More information in Section 2.2.2.1 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help conHyperplane

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>