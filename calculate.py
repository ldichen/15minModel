import networkx as nx
import geopandas as gpd
import numpy as np
from scipy.spatial import KDTree
from shapely.geometry import Point
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter



def calculate_distance(start_node):
    # 构建 KDTree
    graph_nodes = np.array(list(node_positions.values()))
    kdtree = KDTree(graph_nodes)

    # 查找给定坐标的最近节点
    dist, idx = kdtree.query(start_node)
    start_nearest_node = list(node_positions.keys())[idx]

    # 使用Dijkstra算法计算从起始节点到其他节点的最短路径
    lengths, paths = nx.single_source_dijkstra(G_loaded, source=start_nearest_node, cutoff=1000,
                                               weight='weight')
    return lengths, paths

def calculate_node(lengths, paths):
    poi_true_counts = {}
    # 过滤出在距离阈值内的节点
    reachable_nodes = [node for node, length in lengths.items() if length <= 1000 if node in poi_positions]

    # 获取每个 POI 的路径
    poi_paths = {node: paths[node] for node in reachable_nodes}
    # 遍历每个POI路径
    for poi, path in poi_paths.items():
        true_count = 0
        for node in path:
            if G_loaded.nodes[node].get('node') == 'true':
                true_count += 1
        poi_true_counts[poi] = max(true_count - 3, 0)

    # 遍历 poi_true_counts 并更新 lengths 中的对应节点长度
    for node, true_count in poi_true_counts.items():
        # 计算加去的长度
        plus = true_count * 10
        # 更新 lengths 中对应节点的值
        new_length = lengths[node] + plus
        # 确保长度不小于 5
        lengths[node] = max(new_length, 5)

    reachable_cross_nodes = [node for node, length in lengths.items() if length <= 1000 if node in poi_positions]

    # 创建GeoDataFrame来存储这些节点的位置
    geometry = [Point(poi_positions[node]) for node in reachable_cross_nodes]

    start_geometry = [Point(start_node)]

    # 创建图形和轴
    fig, ax = plt.subplots()

    roads.plot(ax=ax, color='blue', linewidth=0.5, label='Road Network')

    gdf = gpd.GeoDataFrame({'geometry': geometry})
    #
    start_gdf = gpd.GeoDataFrame({'geometry': start_geometry})
    start_gdf.plot(ax=ax, marker='x', color='green', markersize=200, label='Start Node')
    gdf.set_crs(epsg=2473, inplace=True)
    # 绘制数据
    gdf.plot(ax=ax, marker='o', color='red', markersize=5, label='node')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    G_loaded = nx.read_graphml('graph.graphml')

    roads = gpd.read_file('C:/Users/admin/Desktop/hp_roaln.geojson')
    # 初始化存储POI路径中true节点数量的字典
    # 获取节点坐标
    node_positions = {node: (data['x'], data['y']) for node, data in G_loaded.nodes(data=True)}
    poi_positions = {node: (data['x'], data['y']) for node, data in G_loaded.nodes(data=True) if
                     data.get('poi') == 'true'}
    # poi_positions = {node: pos for node, pos in node_positions.items() if pos is not None}
    start_node = (121.4799148, 31.2285818)
    lengths, paths = calculate_distance(start_node)
    calculate_node(lengths, paths)

# # 初始化绘图
# plt.figure(figsize=(10, 8))
#
# # 遍历并绘制每条路径
# for path in paths.values():
#     path_edges = list(zip(path[:-1], path[1:]))  # 将路径转化为边
#     nx.draw_networkx_edges(G_loaded, pos=node_positions, edgelist=path_edges, edge_color='red', width=2)
#
# # 设置图形的标题
# plt.title("Shortest Paths from Source Node")
#
# # 显示图像
# plt.show()