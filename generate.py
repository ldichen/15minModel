import geopandas as gpd
import networkx as nx
from networkx.algorithms.shortest_paths.weighted import dijkstra_path_length

from math import sqrt
from scipy.spatial import KDTree
from osgeo import ogr, osr
import numpy as np

# 读取 GeoJSON 文件
roads = gpd.read_file('C:/Users/admin/Desktop/hp_roaln.geojson')
pois = gpd.read_file('C:/Users/admin/Desktop/hp_poi.geojson')

# 定义源坐标系
source_ref = osr.SpatialReference()
source_ref.ImportFromEPSG(4326)

# 定义目标坐标系
target_ref = osr.SpatialReference()
target_ref.ImportFromEPSG(2437)


def calculate_distance_in_meters(Point1, Point2, Source_ref=source_ref, Target_ref=target_ref):
    """
    将两个地理坐标 (经纬度) 转换为投影坐标，并计算它们之间的距离
    """

    # 创建坐标转换对象
    transform = osr.CoordinateTransformation(Source_ref, Target_ref)

    # 创建并转换第一个点
    point1 = ogr.Geometry(ogr.wkbPoint)
    point1.AddPoint(Point1[1], Point1[0])
    point1.Transform(transform)

    # 创建并转换第二个点
    point2 = ogr.Geometry(ogr.wkbPoint)
    point2.AddPoint(Point2[1], Point2[0])
    point2.Transform(transform)

    # 计算距离
    distance = point1.Distance(point2)

    return distance


num = 0


def construct_G():
    global num
    # 构建道路网络图
    G = nx.Graph()
    node_index = {}
    for _, row in roads.iterrows():
        geometry = row.geometry
        if geometry.geom_type == 'MultiLineString':
            for linestring in geometry.geoms:
                coords = list(linestring.coords)
                for i in range(len(coords) - 1):
                    if coords[i] in node_index:  # 如果线的第一个点存在相同的节点
                        if coords[i + 1] in node_index:  # 如果线的第二个点存在相同的节点
                            # 节点均存在，沟通线
                            distance = calculate_distance_in_meters(coords[i], coords[i + 1])
                            G.add_edge(node_index[coords[i]], node_index[coords[i + 1]], weight=distance)
                        else:  # 如果线的第二个点不存在相同的节点
                            if i == len(coords) - 2:
                                G.add_node(num, x=coords[i + 1][0], y=coords[i + 1][1],poi='false',node='true')
                                node_index[coords[i + 1]] = num
                            else:
                                G.add_node(num, x=coords[i + 1][0], y=coords[i + 1][1],poi='false',node='false')
                                node_index[coords[i + 1]] = num
                            # 计算边的长度作为权重
                            distance = calculate_distance_in_meters(coords[i], coords[i + 1])
                            G.add_edge(node_index[coords[i]], num, weight=distance)
                            num = num + 1
                    else:  # 如果线的第一个点不存在相同的节点
                        if coords[i + 1] in node_index:  # 如果线的第二个点存在相同的节点
                            if i == 0:
                                G.add_node(num, x=coords[i][0], y=coords[i][1],poi='false',node='true')
                                node_index[coords[i]] = num
                            # 计算边的长度作为权重
                            distance = calculate_distance_in_meters(coords[i], coords[i + 1])
                            G.add_edge(num, node_index[coords[i + 1]], weight=distance)
                            num = num + 1
                        else:  # 如果线的第二个点不存在相同的节点
                            if i == 0:
                                G.add_node(num, x=coords[i][0], y=coords[i][1],poi='false',node='true')
                                node_index[coords[i]] = num
                            if i == len(coords) - 2:
                                G.add_node(num + 1, x=coords[i + 1][0], y=coords[i + 1][1],poi='false',node='true')
                                node_index[coords[i + 1]] = num + 1
                            else:
                                G.add_node(num + 1, x=coords[i + 1][0], y=coords[i + 1][1],poi='false',node='false')
                                node_index[coords[i + 1]] = num + 1
                            # 计算边的长度作为权重
                            distance = calculate_distance_in_meters(coords[i], coords[i + 1])
                            G.add_edge(num, num + 1, weight=distance)
                            num = num + 2
        elif geometry.geom_type == 'LineString':
            coords = list(geometry.coords)
            for i in range(len(coords) - 1):
                # 计算边的长度作为权重
                G.add_node(num, x=coords[i][0], y=coords[i][1],poi='false')
                # 计算边的长度作为权重
                distance = calculate_distance_in_meters(coords[i], coords[i + 1])
                G.add_edge(num, num + 1, weight=distance)
                num = num + 1
    # 构建路网节点的KD树
    graph_nodes = np.array([[data['x'], data['y']] for _, data in G.nodes(data=True) if data.get('node') == 'true'])
    kdtree = KDTree(graph_nodes)

    for _, row in pois.iterrows():
        poi_coords = row.geometry.coords[0]  # 假设 POI 是一个点
        # 查找最近的节点
        dist, idx = kdtree.query(poi_coords)
        nearest_node = tuple(graph_nodes[idx])
        if poi_coords not in graph_nodes:
            G.add_node(num, x=poi_coords[0], y=poi_coords[1],poi='true',node='true')
            # 计算边的长度作为权重
            distance = calculate_distance_in_meters(poi_coords, nearest_node)
            G.add_edge(num, node_index[nearest_node], weight=distance)
            num = num + 1
    # 保存为 GraphML 文件
    nx.write_graphml(G, 'graph.graphml')


if __name__ == '__main__':
    construct_G()
