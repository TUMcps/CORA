function res = test_zonotope_containsPointWithRobustness

tolerance = 1e-5;

z = zonotope([0 1 0;0 0 1]);

points = [
     2,  -1,   0, 4;
     1, 0.4, 0.5, 5];
expectedRobustness = [
    -1,   0, 0.5,  -5];
expectedDirections = [
     1,  -1,   0, 0.6;
     0,   0,   1, 0.8];

robustness = 0*expectedRobustness;
directions = 0*expectedDirections;
for i=1:size(points, 2)
    [r, d] = containsPointWithRobustness(z, points(:,i));
    robustness(i) = r;
    directions(:, i) = d;
end

directionsOk = abs(directions - expectedDirections) < tolerance;
robustnessOk = abs(robustness - expectedRobustness) < tolerance;

res = all(all(directionsOk)) && all(robustnessOk);

if res
    disp('test_zonotope_containsPointWithRobustness successful');
else
    disp('test_zonotope_containsPointWithRobustness failed');
end

