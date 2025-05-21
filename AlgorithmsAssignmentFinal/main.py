
try:
    import heapq
    import math
    from typing import List
    import networkx as nx
    import matplotlib.pyplot as plt
    from itertools import combinations
except ImportError:
    print("Imports could not be found. Some functionalities may not work.")

import functools

@functools.total_ordering # Decorator to automatically generate other comparison methods
class Edge:
    '''
    We decided to create the edge class in the main.py file so that this code can be copied
    and pasted into jupyter notebook and run with ease if no imports are available
    to be downloaded.
    We used the Edge class to represent all edges between nodes and the source node to
    implement with them Prim's Algorithm
    '''

    def __init__(self, stnode: int, endnode: int, weight: int):

        self.stnode = stnode
        self.endnode = endnode
        self.weight = weight

    def __lt__(self, other: 'Edge') -> bool:    # for putting the edges in the priority list sorted by smallest cost

        return self.weight < other.weight

    def __eq__(self, other: object) -> bool:

        if not isinstance(other, Edge):
            return NotImplemented
        return (self.stnode == other.stnode and
                self.endnode == other.endnode and
                self.weight == other.weight)


def main():
    global path # to fix the error of path variable shown below
    try:
        path = input("Write the path of input file: ")

        with open(path, 'r') as filereader:

            first_line = filereader.readline().strip().split()  #creating an array of the values in first line of file
            if len(first_line) < 4:
                print("Error: The first line of the input file should contain dynos, bonds, bucketcost, and bondcost.")
                return

            dynos = int(first_line[0])  #number of dynos from file
            bonds = int(first_line[1])  #number of bonds from file
            bucketcost = int(first_line[2]) #cost of bucket from file
            bondcost = int(first_line[3])   #cost of bond from file

            num_nodes = dynos + 1   #number of dynos plus the virtual bucket source node

            graph: List[List[Edge]] = [[] for _ in range(num_nodes)] #creating the adjacency list for graph representation

            for _ in range(bonds):
                line = filereader.readline().strip().split()    #reads each line below first line to create bidirectional edges with bondcost weight and put them on the adjacency list
                if len(line) < 2:
                    print(f"Error: Bond definition line is malformed: {' '.join(line)}")
                    return
                u, v = int(line[0]), int(line[1])
                graph[u].append(Edge(u, v, bondcost))
                graph[v].append(Edge(v, u, bondcost))

            for i in range(1, dynos + 1):
                '''
                Establishing a bidirectional edge to all the dynos with the bucket source node 0, for fixing the
                the disconnected graph situation in order to apply Prim's algorithm and finding the MST cost. These edges
                always weigh as much as the bucketcost
                '''
                graph[0].append(Edge(0, i, bucketcost))
                graph[i].append(Edge(i, 0, bucketcost))

    except FileNotFoundError:
        print(f"Error: The file '{path}' was not found.")

    #This is where Prim's Algorithm starts
    cost: List[float] = [math.inf] * num_nodes  #a list of minimum costs for each node
    prev: List[int] = [-1] * num_nodes  # a list of parents for each node
    visited: List[bool] = [False] * num_nodes   # initializing all values to not visited (false)

    start_node = 0  # always starting from the virtual source
    cost[start_node] = 0    # cost of start node put to 0 since it starts from there ("it's own parent")

    mst_queue: List[Edge] = []  #the priority queue created and empty to store the edges from the ones with smallest cost to largest

    for edge in graph[start_node]:
        heapq.heappush(mst_queue, edge)    #pushing all edges connected to node 0, which all have the same weight (bucketcost)

    visited[start_node] = True  #the starting node is visitet, it doesn't have to be visited again,

    mst_edges: List[Edge] = []
    mst_cost: int = 0   #the overall minimum cost is initialized to 0
    mst_bond_count: int = 0  # the number of bonds should always be less or equal to number of dynos, in evey scenario with every price of bucket hosting and bonds

    while mst_queue and mst_bond_count < num_nodes - 1:
        current_edge = heapq.heappop(mst_queue) #it iteratively takes each edge in priorty with minimum cost

        v = current_edge.stnode
        z = current_edge.endnode
        weight = current_edge.weight

        if visited[z]:  #if the dyno is visited which means we have found its minimum edge, then skip
            continue

        visited[z] = True #if not visited, mark as visited

        cost[z] = weight    # add its cost to the minimum cost list for each dyno
        prev[z] = v         # add its parent for the MST graph

        mst_edges.append(current_edge)  # add the edge to the list of edges of MST
        mst_cost += weight  #add the weight to the mst cost
        mst_bond_count += 1 #added 1 more bond to the graph

        for neighbor_edge in graph[z]:
            '''
            Now we push in the priority queue all the other edges connected to the node z 
            which are not connected to already visited dynos. For eg. if dyno z = 1 has 
            an edge to dyno 2, its cost is smaller than 0 to 2 so it goes in the beginning 
            of queue based on priority, ensuring finding the mst for all the nodes.
            In case of a disconnected graph, we know each component will start from virtual node 0,
            so we will always have at least a bucket hosted for each component
            '''
            if not visited[neighbor_edge.endnode]:
                heapq.heappush(mst_queue, neighbor_edge)

    if mst_bond_count < num_nodes - 1 and num_nodes > 1:
        pass

    print(f"Minimum Spanning Tree Cost: {mst_cost}")    # Printing the value of MST cost

    try:
        with open(path) as f:
            n, m, bc, bdc = map(int, f.readline().split())
            G = nx.Graph()
            for _ in range(m):
                u, v = map(int, f.readline().split())
                G.add_edge(u, v, weight=bdc)
            for i in range(1, n + 1):
                G.add_edge(0, i, weight=bc)

        T = nx.minimum_spanning_tree(G)
        mst_weight = sum(d['weight'] for _, _, d in T.edges(data=True))

        def count_all_msts(G, weight):
            count = 0
            for subset in combinations(G.edges(data=True), len(G.nodes) - 1):
                if sum(d['weight'] for _, _, d in subset) == weight:
                    H = nx.Graph()
                    H.add_nodes_from(G.nodes)
                    H.add_edges_from((u, v) for u, v, _ in subset)
                    if nx.is_connected(H): count += 1
            return count

        total_msts = count_all_msts(G, mst_weight)
        print("Total MSTs with same minimum cost:", total_msts)

        real_nodes = [node for node in G.nodes if node != 0]
        pos = nx.spring_layout(G.subgraph(real_nodes), seed=42)

        bucket_nodes = {v if u == 0 else u for u, v in T.edges() if 0 in (u, v)}

        node_colors = ['red' if node in bucket_nodes else 'skyblue' for node in real_nodes]

        nx.draw_networkx_nodes(G, pos, nodelist=[node for node, color in zip(real_nodes, node_colors) if color == 'red'], node_color='red', label="Dyno linked to a Bucket via a Bond")
        nx.draw_networkx_nodes(G, pos, nodelist=[node for node, color in zip(real_nodes, node_colors) if color == 'skyblue'], node_color='skyblue', label="Dyno Hosting a Bucket")
        nx.draw_networkx_labels(G, pos, labels={n: n for n in real_nodes})

        mst_edges = []
        for u, v in T.edges():
            if 0 not in (u, v):
                mst_edges.append((u, v))
        nx.draw_networkx_edges(T, pos, edgelist=mst_edges, edge_color='green', width=2, label="MST Bond")

        plt.title(f"One MST Visualization\nTotal distinct MSTs with same cost: {total_msts}")
        plt.legend()
        plt.show()

    except Exception as e:
        print(f" Imports not found to complete the graph visualization and count the number of possible MSTs. Try importing networkx and matplotlib.pyplot")

if __name__ == "__main__":
    main()
