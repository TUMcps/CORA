# Polytope

A polytope in halfspace representation is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BP%7D:=%5C%7Bx%5Cin%5Cmathbb%7BR%7D%5En%5C,%7C%5C,A%20x%5Cleq%20b%5C%7D,%5Cquad%20A%5Cin%5Cmathbb%7BR%7D%5E%7Bh%5Ctimes%20n%7D,b%5Cin%5Cmathbb%7BR%7D%5Eh."/>
</p>

A polytope in vertex representation is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BP%7D:=%5Cbigg%5C%7B%5Csum_%7Bi=1%7D%5Es%5Cbeta_i%20v_i%5C,%5Cbigg%7C%5C,%5Csum_%7Bi=1%7D%5Es%5Cbeta_i=1,%5Cbeta_i%5Cgeq%200%5Cbigg%5C%7D."/>
</p>

Polytopes are closed and convex sets.

In CORA, polytopes are instantiated in halfspace representation by

    P = polytope(A,b);

where ``A`` is the inequality constraint matrix, ``b`` is the inequality offset vector, and in vertex representation by

    P = polytope(V);

where ``V`` is a point cloud.

Example:

    % initialize polytope
    P = polytope([1 1; -1 1; 0 -1],[1;1;1]);
    % plot polytope
    plot(P);

More information in Section 2.2.1.4 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual.html">CORA manual</a> or type

    help polytope

<hr style="height: 1px;">
<img src="../../app/images/coraLogo_readme.svg"/>