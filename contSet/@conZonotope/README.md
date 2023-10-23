# Constrained Zonotope

A constrained zonotope is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BCZ%7D:=%5Cbigg%5C%7Bc&plus;%5Csum_%7Bi=1%7D%5E%5Cgamma%5Cbeta_i%20G_%7B(%5Ccdot,i)%7D%5C,%5Cbigg%7C%5C,%5Csum_%7Bi=1%7D%5E%5Cgamma%5Cbeta_i%20A_%7B(%5Ccdot,i)%7D=b,%5Cbeta_i%5Cin%5B-1,1%5D%5Cbigg%5C%7D."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{CZ} := \bigg\{ c + \sum_{i=1}^\gamma \beta_i G_{(\cdot,i)} \, \bigg| \, \sum_{i=1}^\gamma \beta_i A_{(\cdot,i)} = b, \beta_i \in [-1,1] \bigg\}.
-->

Constrained zonotopes are compact, convex, and bounded sets.

In CORA, constrained zonotopes are instantiated by

    cZ = conZonotope(c,G,A,b);

where ``c`` is the zonotope center, ``G`` is the generator matrix, ``A`` is the constraint matrix, and ``b`` is the constraint offset.

Example:

    % initialize constrained zonotope
    cZ = conZonotope([0;0],[1 0 1; 1 2 -1],[-2 1 1],2);
    % plot constrained zonotope
    plot(cZ);

More information in Section 2.2.1.9 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help conZonotope

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>