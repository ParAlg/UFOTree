### Euler Tour Tree Interface

| Operation           | Description                                               |
|---------------------|-----------------------------------------------------------|
| `Link(u, v)`        | Connect two trees by adding the edge between `u` and `v`. |
| `InnerLink(u, v)`   | Connect two trees by adding the edge between `u` and `v`. |
| `OuterLink(u, v)`   | Connect two trees by adding the edge between `u` and `v`. |
| `Cut(u, v)`         | Remove the edge between `u` and `v`, splitting the tree.  |
| `IsConnected(u, v)` | Return `True` if `u` and `v` are in the same component.   |
|---------------------|-----------------------------------------------------------|

`LinkInner` splits the tour to the right of `u` and to the left of `v`, putting `v` and its whole component right after `u` in the new Euler tour.
Thus the new order of things is `u`, `uv`, `v`, the rest of `v`'s component, `vu`.

`LinkOuter` splits the tour to the left of `u` and to the left of `v`, putting `v` and its whole component right before `u` in the new Euler tour.
Thus the new order of things is `uv`, `v`, the rest of `v`'s component, `vu`, `u`.

Both of these essentially keep `v` as the root of its component by being the first occurence in its subtour, which is useful in some applications.
For example, if `v` should be a child of `u`, using `LinkInner(u,v)` may be useful.