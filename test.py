import networkx as nx
G_loaded = nx.read_graphml('graph.graphml')
is_connected = nx.is_connected(G_loaded)
print("图是否连通:", is_connected)
components = list(nx.connected_components(G_loaded))  # 对于无向图
print("连通分量:", components)