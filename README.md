## Example

1. Generate (X,Y) using a linear dynamic system: `generate_lds.m`

```
X(k+1) = A*X(k) + B*U(k) + v(k)
Y(k) = C*X(k) + D*U(k) + w(k)
Ev = 0, Evv' = Q
Ew = 0, Eww' = R
```

2. Perfrom Kalman filtering to estimate X from Y: `kalman_filt.m`

