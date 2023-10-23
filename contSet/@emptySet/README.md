# Empty Set

An empty set is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D:=%5Cemptyset."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{O} := \emptyset .
-->

In CORA, empty sets are instantiated by

    O = emptySet(n);

where ``n`` is the dimension.
This dimension is solely used for internal purposes, such as dimensionality checks.

Example:

    % initialize empty set
    O = emptySet(2);
    % plot empty set
    plot(O);

Empty sets are used in hybrid systems to model instant transitions.
More information in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help emptySet

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>