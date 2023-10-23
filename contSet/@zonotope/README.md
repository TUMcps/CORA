# Zonotope

A zonotope is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BZ%7D:=%5Cbigg%5C%7Bc&plus;%5Csum_%7Bi=1%7D%5E%5Cgamma%5Cbeta_i%20G_%7B(%5Ccdot,i)%7D%5C,%5Cbigg%7C%5C,%5Cbeta_i%5Cin%5B-1,1%5D%5Cbigg%5C%7D."/>
</p>
<!--
for editor.codecogs.com:
\mathcal{Z} := \bigg\{ c + \sum_{i=1}^\gamma \beta_i G_{(\cdot,i)} \, \bigg| \, \beta_i \in [-1,1] \bigg\} .
-->

Zonotopes are compact, convex, and bounded sets.
They are point-symmetric with respect to their center.

In CORA, zonotopes are instantiated by

    Z = zonotope(c,G);

where ``c`` is the center vector and ``G`` is the generator matrix.

Example:

    % initialize zonotope
    Z = zonotope([1;2],[1 0 1; 0 1 1]);
    % plot zonotope
    plot(Z);

More information in Section 2.2.1.1 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help zonotope

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>