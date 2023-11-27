# Constrained Polynomial Zonotope

A constrained polynomial zonotope is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BCPZ%7D:=%5Cbigg%5C%7Bc&plus;%5Csum_%7Bi=1%7D%5Eh%5Cbigg(%5Cprod_%7Bk=1%7D%5Ep%5Calpha_k%5E%7BE_%7B(k,i)%7D%7D%5Cbigg)G_%7B(%5Ccdot,i)%7D&plus;%5Csum_%7Bj=1%7D%5Ed%5Cbeta_j%20G_%7BI(%5Ccdot,j)%7D%5C,%5Cbigg%7C"/>
</p>

<p style="background-color: white; padding-left:80px">
<img src="https://latex.codecogs.com/svg.image?%5Csum_%7Bi=1%7D%5Eq%5Cbigg(%5Cprod_%7Bk=1%7D%5Ep%5Calpha_k%5E%7BR_%7B(k,i)%7D%7D%5Cbigg)A_%7B(%5Ccdot,i)%7D=b,%5Calpha_k,%5Cbeta_j%5Cin%5B-1,1%5D%5Cbigg%5C%7D."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{CPZ} := \bigg\{ c + \sum_{i=1}^h \bigg( \prod_{k=1}^p \alpha_k^{E_{(k,i)}} \bigg) G_{(\cdot,i)} + \sum_{j=1}^d \beta_j G_{I(\cdot,j)} \, \bigg|
\sum_{i=1}^q \bigg( \prod_{k=1}^p \alpha_k^{R_{(k,i)}} \bigg) A_{(\cdot,i)} = b, \alpha_k, \beta_j \in [-1,1] \bigg\} .
-->

Constrained polynomial zonotopes can describe non-convex sets.

In CORA, constrained polynomial zonotopes are instantiated by

    cPZ = conPolyZono(c,G,E,A,b,R,Grest);

where ``c`` is the start point, ``G`` is the generator matrix with dependent factors, ``E`` is the exponent matrix, ``A`` is the constraint matrix, ``b`` is the constraint offset, ``R`` is the constraint exponent matrix, and ``Grest`` is the  generator matrix with independent factors.

Example:

    % initialize constrained polynomial zonotope
    cPZ = conPolyZono([0;0],[1 0 1 -1; 0 1 1 1],[1 0 1 2; 0 1 1 0; 0 0 1 1],...
        [1 -0.5 0.5],0.5,[0 1 2; 1 0 0; 0 1 0]);
    % plot polynomial zonotope
    plot(cPZ);

More information in Section 2.2.1.6 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help conPolyZono

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>