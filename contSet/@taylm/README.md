# Taylor Models

A Taylor modelis defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BT%7D(x):=%5C%7Bp(x)&plus;y%5C,%7C%5C,y%5Cin%5Cmathcal%7BI%7D%5C%7D." />
</p>

<!--
for editor.codecogs.com: 
\mathcal{T}(x) := \{ p(x) + y \, | \, y \in \mathcal{I}   \} .
-->

Taylor models are compact and can describe non-convex sets.

In CORA, Taylor models are instantiated by

    tay = taylm(int,max_order,name);

where ``int`` is an domain for the variable, ``max_order`` is the maximum polynomial degree, and ``name`` is the name of the variable.

Example:

    % initialize Taylor model
    tx = taylm(interval(-2,2),4,'x');
    tay = [tx; sin(tx)];
     
    % plot Taylor model
    plot(tay);

More information in Section 2.2.3.1 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help taylm

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>