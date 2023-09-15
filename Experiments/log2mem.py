import os, csv, dataclasses, argparse

tested_systems = ["GraphPi", "Dryadic", "GraphMini", "Base"]
tested_graphs = ["wiki", "patents", "youtube", "lj", "orkut", "friendster"]
# tested_graphs = ["wiki"]
tested_queries = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]
tested_query_types = ["VertexInduced", "EdgeInduced", "EdgeInducedIEP"]

@dataclasses.dataclass
class Memlog:
    is_available: bool = False
    is_oot: bool = False
    vertex_set_alloc: str = "None"
    aux_graph_alloc: str = "None"
    num_threads: int = 64
    graph: str = "None"
    query: str = "None"
    qtype: str = "None"
    number_runs: int = 0
    result: int = 0
    
def read_mem(system, graph, query, qtype):
    log = Memlog()
    log.system = system
    log.graph = graph
    log.query = query
    log.qtype = qtype
    
    dir_path = os.path.dirname(__file__)
    file_path = os.path.join(dir_path, system, graph, query, qtype, "run.log")
    
    try:
        with open(file_path, "r") as file:
            for line in file:
                if "VertexSetAllocated=" in line:
                    log.vertex_set_alloc = line.split("VertexSetAllocated=")[1].strip()
                if "MiniGraphAllocated=" in line:
                    log.aux_graph_alloc = line.split("MiniGraphAllocated=")[1].strip()
                if "Threads=" in line:
                    log.num_threads = int(line.split("Threads=")[1].strip())
                if "RESULT=" in line:
                    log.result = int(line.split("RESULT=")[1].strip())
                    log.number_runs += 1
            
            if log.number_runs >= 1:
                log.is_available = True
                log.is_oot = False
            else:
                log.is_available = True
                log.is_oot = True

    except FileNotFoundError:
        # the file does not exist -> system does not support it
        log.is_available = False
    
    return log

def find_log(logs, system, graph, query, qtype) -> Memlog:
    for log in logs:
        if log.system == system and log.graph == graph and log.query == query and log.qtype == qtype:
            return log
    return Memlog() # not available
                    
def log2mem(baseline:str, qtype: str):
    
    def read_logs(system: str):
        logs = []
        for graph in tested_graphs:
            for query in tested_queries:
                log = read_mem(system, graph, query, qtype)
                logs.append(log)
        return logs
    dir_path = os.path.dirname(__file__)
    file_path = os.path.join(dir_path, f"{baseline}_{qtype}_memory.csv")
    headers = ["Graph", "Query", "QueryType", "VertexSetAllocated", "AuxGraphAllocated", "AuxGraphPerThread", "Threads"]
    
    logs = read_logs(baseline)
    lines = []
    for graph in tested_graphs:
        for query in tested_queries:
            log = find_log(logs, baseline, graph, query, qtype)
            cur_line = [graph, query, qtype]
            if log.is_available and not log.is_oot:
                auxgraph = float(log.aux_graph_alloc.split()[0])
                perthread_auxgraph = round(auxgraph / log.num_threads, 3)
                cur_line.append(str(log.vertex_set_alloc))
                cur_line.append(str(log.aux_graph_alloc))
                cur_line.append(str(perthread_auxgraph) + " MB")
                cur_line.append(str(log.num_threads))
                
            elif log.is_available and log.is_oot:
                cur_line.append("out of time")
                cur_line.append("out of time")
                cur_line.append("out of time")
                cur_line.append("out of time")
            else:
                cur_line.append("not available")
                cur_line.append("not available")
                cur_line.append("not available")
                cur_line.append("not available")
            lines.append(cur_line)
        
    with open(file_path, "w") as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        writer.writerows(lines)
    
    print(f"Result is writtern to {file_path}")
    
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='convert log file to speed up')
    parser.add_argument('--adjtype', default=str, type=str, help='Query Type', choices=["VertexInduced", "EdgeInduced", "EdgeInducedIEP"])
    args = parser.parse_args()
    
    log2mem("GraphMini", args.adjtype)    