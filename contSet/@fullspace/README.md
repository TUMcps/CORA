# Fullspace

A fullspace is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BFS%7D:=%5Cmathbb%7BR%7D%5En."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{FS} := \mathbb{R}^n .
-->

Fullspaces are unbounded sets.

In CORA, fullspaces are instantiated by

    fs = fullspace(n);

where ``n`` is the dimension.

Example:

    % initialize fullspace
    fs = fullspace(2);
    % plot fullspace
    plot(fs);

Fullspaces are used in hybrid systems to model unbounded invariant sets.
More information in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help fullspace

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>