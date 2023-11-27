# Polynomial Zonotope

A polynomial zonotope is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BPZ%7D:=%5Cbigg%5C%7Bc&plus;%5Csum_%7Bi=1%7D%5Eh%5Cbigg(%5Cprod_%7Bk=1%7D%5Ep%5Calpha_k%5E%7BE_%7B(k,i)%7D%7D%5Cbigg)G_%7B(%5Ccdot,i)%7D&plus;%5Csum_%7Bj=1%7D%5Eq%5Cbeta_j%20G_%7BI(%5Ccdot,j)%7D%5C,%5Cbigg%7C%5C,%5Calpha_k,%5Cbeta_j%5Cin%5B-1,1%5D%5Cbigg%5C%7D." />
</p>

<!--
for editor.codecogs.com: 
\mathcal{PZ} := \bigg\{ c + \sum_{i=1}^h \bigg( \prod_{k=1}^p \alpha_k^{E_{(k,i)}} \bigg) G_{(\cdot,i)} + \sum_{j=1}^q \beta_j G_{I(\cdot,j)} \, \bigg| \, \alpha_k, \beta_j \in [-1,1] \bigg\} .
-->

Polynomial zonotopes are compact and can describe non-convex sets.

In CORA, polynomial zonotopes are instantiated by

    pZ = polyZonotope(c,G,Grest,E);

where ``c`` is the start point, ``G`` is the generator matrix with dependent factors, ``Grest`` is the  generator matrix with independent factors, and ``E`` is the exponent matrix.

Example:

    % initialize polynomial zonotope
    pZ = polyZonotope([4; 4],[2 1 2; 0 2 2],[1; 0],[1 0 3; 0 1 1]);
    % plot polynomial zonotope
    plot(pZ);

More information in Section 2.2.1.5 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help polyZonotope

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>