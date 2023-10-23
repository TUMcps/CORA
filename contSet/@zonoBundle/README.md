# Zonotope Bundle

A zonotope bundle is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BZB%7D:=%5Cbigcap_%7Bi=1%7D%5Es%5Cmathcal%7BZ%7D_i,"/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{ZB} := \bigcap_{i=1}^s \mathcal{Z}_i ,
-->

where each ``Z_i`` is a zonotope.

Zonotope bundles are compact, convex, and bounded sets.

In CORA, zonotope bundles are instantiated by

    zB = zonoBundle(list);

where ``list`` is a list (cell-array) of zonotope objects.

Example:

    % initialize zonotope bundle
    zB = zonoBundle({zonotope([1;1],[3 0; 0 2]),zonotope([0;0],[2 2; 2 -2])});
    % plot zonotope bundle
    plot(zB);

More information in Section 2.2.1.8 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help zonoBundle

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>