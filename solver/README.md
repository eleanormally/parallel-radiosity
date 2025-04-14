# file format
## Spec
```
[dimension, int]
//followed by [dimension] rows of the following format:
[double x dimension, represents form factor from j to i] [vec3f emmitence] [vec3f reflection]
```
## Example
```
3
0.0 0.5 0.5 0.3 0.3 0.3 0.5 0.1 0.0
0.5 0.0 0.0 0.3 0.3 0.3 0.5 0.1 0.0
0.5 0.5 0.0 0.3 0.3 0.3 0.5 0.1 0.0
```
