# Template

## Fast IO
```cpp
std::ios_base::sync_with_stdio(false);
std::cin.tie(0);
```

## Disjoint Set Union
```cpp
struct DSU {
    vector<int> parent;
    vector<int> rank;
    DSU(int n) {
        parent.assign(n+5, 0);
        rank.assign(n+5, 0);
        for (int i=1;i<=n;i++) { parent[i] = i; }
    }
    int find(int u) {
        if (parent[u] == u) { return u; }
        return parent[u] = find(parent[u]);
    }
    void union_sets(int u, int v) {
        u = find(u);
        v = find(v);
        if (u != v) {
            if (rank[u] < rank[v]) { swap(u, v); }
            parent[v] = u;
            rank[u] += rank[v];
        }
    }
};
```

## Segment Tree
```cpp
// Example: max point update & max range query
struct SegTree {
    vector<int> tree; 
    SegTree(int n) {
        tree.assign(4*n+5, 0); 
    }
    int merge_segment(const int &seg1, const int &seg2) {
        return max(seg1, seg2);
    }
    void build(int a[], int v, int tl, int tr) {
        if (tl == tr) {
            tree[v] = a[tl]; 
        } else {
            int tm = (tl + tr) >> 1;
            build(a, v*2, tl, tm);
            build(a, v*2+1, tm+1, tr);
            tree[v] = merge_segment(tree[v*2], tree[v*2+1]); 
        }
    }
    void update(int v, int tl, int tr, int pos, int val) {
        if (tl == tr) {
            tree[v] = max(tree[v], val); 
        } else {
            int tm = (tl + tr) >> 1;
            if (pos <= tm) { update(v*2, tl, tm, pos, val); }
            else { update(v*2+1, tm+1, tr, pos, val); }
            tree[v] = merge_segment(tree[v*2], tree[v*2+1]); 
        }
    }
    int query(int v, int tl, int tr, int l, int r) { 
        if (tr < l || r < tl) { return 0; }
        if (l <= tl && tr <= r) {
            return tree[v];
        }
        int tm = (tl + tr) >> 1;
        return merge_segment(
            query(v*2, tl, tm, l, r),
            query(v*2+1, tm+1, tr, l, r)
        );
    }
};
```

## Lazy Propagation Segment Tree
```cpp
// Example: sum range update & range query
struct LazySegTree {
    vector<int> tree; 
    vector<int> lazy; 
    LazySegTree(int n) {
        tree.assign(4*n+5, 0); 
        lazy.assign(4*n+5, 0); 
    }
    int merge_segment(const int &seg1, const int &seg2) {
        return seg1 + seg2;
    }
    void apply(int v, int x, int tl, int tr) {
        tree[v] += x * (tr - tl + 1);
        lazy[v] += x;
    }
    void push(int v, int tl, int tr) {
        if (!lazy[v]) { return; } 
        int tm = (tl + tr) >> 1;
        apply(v*2, lazy[v], tl, tm);
        apply(v*2+1, lazy[v], tm+1, tr);
        lazy[v] = 0; 
    }
    void update(int v, int tl, int tr, int l, int r, int val) {
        if (tr < l || r < tl) { return; }
        if (l <= tl && tr <= r) {
            apply(v, val, tl, tr);
        } else {
            push(v, tl, tr);
            int tm = (tl + tr) >> 1;
            update(v*2, tl, tm, l, r, val); 
            update(v*2+1, tm+1, tr, l, r, val); 
            tree[v] = merge_segment(tree[v*2], tree[v*2+1]); 
        }
    }
    int query(int v, int tl, int tr, int l, int r) { 
        if (tr < l || r < tl) { return 0; }
        if (l <= tl && tr <= r) {
            return tree[v];
        }
        push(v, tl, tr);
        int tm = (tl + tr) >> 1;
        return merge_segment(
            query(v*2, tl, tm, l, r),
            query(v*2+1, tm+1, tr, l, r)
        );
    }
};
```

## Heavy-Light Decomposition
```cpp
vector<int> adj[MXN];
int sz[MXN], parent[MXN], depth[MXN], tin[MXN], tout[MXN], top[MXN];
int n;
void dfs1(int u, int p) {
    sz[u] = 1;
    vector<int> children;
    for (auto v : adj[u]) {
        if (v == p) { continue; }
        depth[v] = depth[u] + 1;
        dfs1(v, u);
        children.emplace_back(v);
        sz[u] += sz[v];
        parent[v] = u;
    }
    swap(adj[u], children);
}
int dfs2(int u, int timer) {
    tin[u] = timer;
    for (auto v : adj[u]) {
        if (v == adj[u][0]) {
            top[v] = top[u];
        } else {
            top[v] = v;
        }
        timer = dfs2(v, timer+1);
    }
    tout[u] = timer + 1;
    return timer;
}
void hld(int root) {
    dfs1(root, -1);
    for (int i=1;i<=n;i++) {
        sort(adj[i].begin(), adj[i].end(), [&](const int &u, const int &v) {
            return sz[u] > sz[v];
        });
    }
    top[root] = root;
    dfs2(root, 1);
}
void point_update(SegTree &segtree, int u, int val) {
    segtree.update(1, 1, n, tin[u], val);
}
int query_subtree(SegTree &segtree, int u) {
    return segtree.query(1, 1, n, tin[u], tout[u] - 1);
}
int query_path(SegTree &segtree, int u, int v) {
    int res = 0;
    while (top[u] != top[v]) {
        if (depth[top[u]] < depth[top[v]]) {
            swap(u, v);
        }
        res = segtree.merge_segment(res, segtree.query(1, 1, n, tin[top[u]], tin[u]));
        u = parent[top[u]];
    }
    res = segtree.merge_segment(res, segtree.query(1, 1, n, min(tin[u], tin[v]), max(tin[u], tin[v])));
    return res;
}
```

## Heavy-Light Decomposition with Range Update
```cpp
void range_update(LazySegTree &segtree, int u, int v, int val) {
    while (top[u] != top[v]) {
        if (depth[top[u]] < depth[top[v]]) {
            swap(u, v);
        }
        segtree.update(1, 1, n, tin[top[u]], tin[u], val);
        u = parent[top[u]];
    }
    segtree.update(1, 1, n, min(tin[u], tin[v]), max(tin[u], tin[v]), val);
}
```

## Mo's algorithm
```cpp
#define BLOCK_SIZE 450
void add(int idx) {
}
void remove(int idx) {
}
struct Query {
    int l, r, idx;
    Query(int l, int r, int idx) {
        this->l = l;
        this->r = r;
        this->idx = idx;
    }
    bool operator<(const Query &other) const {
        if (l / BLOCK_SIZE != other.l / BLOCK_SIZE) {
            return make_pair(l / BLOCK_SIZE, r) <
                   make_pair(other.l / BLOCK_SIZE, other.r);
        }
        return (l / BLOCK_SIZE % 2) ? (r < other.r) : (r > other.r);
    }
};
vector<int> mo(int q, vector<Query> &queries) {
    vector<int> ans(q);
    sort(queries.begin(), queries.end());
    int l = 1, r = 0;
    for (Query q : queries) {
        while (r < q.r) { add(++r); }
        while (r > q.r) { remove(r--); }
        while (l > q.l) { add(--l); }
        while (l < q.l) { remove(l++); }
        // ans[q.idx] = ...
    }
    return ans;
}
```

## Strongly Connected Component
```cpp
void dfs(int u, vector<int> adj[],  vector<bool> &vis, vector<int> &output) {
    vis[u] = true;
    for (auto v : adj[u]) {
        if (vis[v]) { continue; }
        dfs(v, adj, vis, output);
    }
    output.emplace_back(u);
}

vector<int> find_scc(int n, vector<int> adj[]) {
    vector<int> scc(n+5);

    vector<int> order;
    vector<bool> vis(n+5);
    for (int v=1;v<=n;v++) {
        if (vis[v]) { continue; }
        dfs(v, adj, vis, order);
    }

    vector<int> adj_T[n+5];
    for (int u=1;u<=n;u++) { adj[u].clear(); }
    for (int u=1;u<=n;u++) {
        for (auto v: adj[u]) {
            adj_T[v].emplace_back(u);
        }
    }

    vis.assign(n+5, false);
    reverse(order.begin(), order.end());
    for (auto u : order) {
        if (vis[u]) { continue; }
        vector<int> component;
        dfs(u, adj_T,  vis, component);
        int root = *min_element(component.begin(), component.end());
        for (auto v : adj[u]) {
            scc[v] = u;
        }
    }

    return scc;
}
```
