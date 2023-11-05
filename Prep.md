## Approaches:

- Re-create the stage 1 subnetworks

```
data[i] = {"edgeCount": edge_count, "subNet": subnetwork}
```

```
df = pd.DataFrame.from_dict(data, orient='index')
```

- 