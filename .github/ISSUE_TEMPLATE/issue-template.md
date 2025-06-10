---
name: Issue Template
about: General questions, bug reports, and feature requests
title: ''
labels: ''
assignees: TUMcps

---

## Description 

*Briefly describe the issue*

## Steps to Reproduce the Issue

*Ideally provide code of a minimal example. E.g.*
```
% 1. init two zonotopes
Z1 = zonotope([1;2], [1 0 1; 0 1 1]);
Z2 = zonotope([5;4], [1 0 -1; 0 1 1]);

% 2. plot the zonotopes
figure;
plot(Z1);
plot(Z2);
```

### Expected Behavior

*both zonotopes are visible in the figure*

### Current Behavior

*only the last zonotope is visible*

## Versions

- CORA version: *v2025.0.0* (type `CORAVERSION` in command window)
- Matlab version: *R2024b* (type `version` in command window)
- OS: *Linux/Windows/Mac*
