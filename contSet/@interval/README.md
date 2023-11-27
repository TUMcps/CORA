# Interval

An interval is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BI%7D:=%5C%7Bx%5Cin%5Cmathbb%7BR%7D%5En%5C,%7C%5C,%5Cunderline%7Bx%7D_i%5Cleq%20x_i%5Cleq%5Coverline%7Bx%7D_i,i=1,...,n%5C%7D."/>
</p>

Intervals are convex sets.

In CORA, intervals are instantiated by

    I = interval(lb,ub);

where ``lb`` is the lower bound and ``ub`` is the upper bound.

Example:

    % initialize interval
    I = interval([-1;1],[2;3]);
    % plot interval
    plot(I);

CORA supports unbounded intervals by setting one or more limits to ``-Inf`` or ``Inf``.

More information in Section 2.2.1.2 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help interval

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>