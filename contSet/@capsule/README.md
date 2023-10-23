# Capsule

A capsule is defined as

<p style="background-color: white;">
<img src="https://latex.codecogs.com/svg.image?%5Cmathcal%7BC%7D:=%5Cmathcal%7BL%7D%5Coplus%5Cmathcal%7BS%7D,%5Cquad%5Cmathcal%7BL%7D:=%5C%7Bc&plus;g%5Calpha%5C,%7C%5C,%5Calpha%5Cin%5B-1,1%5D%5C%7D,%5Cquad%5Cmathcal%7BS%7D=%5C%7Bx%5C,%7C%5C,%7C%7Cx%7C%7C_2%5Cleq%20r%5C%7D."/>
</p>

<!--
for editor.codecogs.com: 
\mathcal{C} := \mathcal{L} \oplus \mathcal{S}, \quad \mathcal{L} := \{ c + g\alpha \, | \, \alpha \in [-1,1] \}, \quad \mathcal{S} = \{ x \, | \, ||x||_2 \leq r \}.
-->

Capsules are compact, convex, and bounded sets.

In CORA, capsules are instantiated by

    C = capsule(c,g,r);

where ``c`` is the center, ``g`` is the generator, and ``r`` is the radius.

Example:

    % initialize capsule
    C = capsule([1;0],[-1;1],0.5);
    % plot capsule
    plot(C);

More information in Section 2.2.1.7 in the <a target='_blank' href="https://tumcps.github.io/CORA/manual">CORA manual</a> or type

    help capsule

<hr style="height: 1px;">

<img src="../../app/images/coraLogo_readme.svg"/>