import argparse

def get_number_edges(query_string: str):
    count = 0
    for char in query_string:
        if char == '1':
            count += 1
    return count

def get_edges(query_size: int, query_string:str):
    edges = []
    for i in range(query_size):
        for j in range(query_size):
            idx = i * query_size + j
            if query_string[idx] == '1':
                edges.append((i, j))
    return get_number_edges(query_string), edges

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script that takes a query size the adjacent matrix and outputs the corresponding format readable by Dryadic.")
    parser.add_argument("query_size", type=int, help="An integer input.")
    parser.add_argument("query_str", type=str, help="A string input.")
    args = parser.parse_args()

    num_edges, edges = get_edges(args.query_size, args.query_str)
    print(1) # number of pattern
    print(args.query_size, num_edges) # number of vertices, number of edges
    print("0 " * args.query_size) # labels
    for edge in edges:
        u, v = edge
        print(u, v)    