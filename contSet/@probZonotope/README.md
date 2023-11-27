# Probabilistic Zonotope

A probabilistic zonotope is a random variable (not a set, strictly speaking) defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BZ%7D_p:=%5Cbigg%5C%7Bc&plus;%5Csum_%7Bi=1%7D%5E%5Cgamma%5Cbeta_i%20G_%7B(%5Ccdot,i)%7D&plus;%5Csum_%7Bj=1%7D%5Eq%5Cmathcal%7BN%7D%5E%7B(j)%7D(0,1)%5C,P_%7B(%5Ccdot,j)%7D%5C,%5Cbigg%7C%5C,%5Cbeta_i%5Cin%5B-1,1%5D%5Cbigg%5C%7D."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{Z}_p := \bigg\{ c + \sum_{i=1}^\gamma \beta_i G_{(\cdot,i)} + \sum_{j=1}^q \mathcal{N}^{(j)}(0,1) \, P_{(\cdot,j)} \, \bigg| \, \beta_i \in [-1,1] \bigg\} .
-->

In CORA, probabilistic zonotopes are instantiated by

    probZ = probZonotope(Z,P);

where ``Z`` is the concatenated zonotope center and generator matrix, and ``P`` is the probabilistic generator matrix.

Example:

    % initialize probabilistic zonotope
    probZ = probZonotope([0 1 0; 0 0 1],[3 2; 3 -2]);
    % plot probabilistic zonotope
    plot(probZ);

More information in Section 2.2.1.10 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help probZonotope

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>