# Ellipsoid

An ellipsoid is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BE%7D:=%5C%7Bx%5Cin%5Cmathbb%7BR%7D%5En%5C,%7C%5C,%5Cell%5E%5Ctop%20x%5Cleq%5Cell%5E%5Ctop&plus;%5Csqrt%7B%5Cell%5E%5Ctop%20Q%5Cell%7D,%5Cforall%5Cell%5Cin%5Cmathbb%7BR%7D%5En%5C%7D,"/>
</p>

where ``Q`` is a positive semi-definite matrix.

Ellipsoids are compact and convex sets.

In CORA, ellipsoids are instantiated by

    E = ellipsoid(Q,q);

where ``Q`` is the positive semi-definite shape matrix and ``q`` is the center.

Example:

    % initialize ellipsoid
    E = ellipsoid([1 -1; -1 2],[2;1]);
    % plot ellipsoid
    plot(E);

More information in Section 2.2.1.3 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help ellipsoid

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>