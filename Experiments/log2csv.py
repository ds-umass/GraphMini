import os, csv, dataclasses

tested_systems = ["GraphPi", "Dryadic", "MiniGraph"]
tested_graphs = ["wiki", "patents", "youtube", "lj", "orkut", "friendster"]
tested_queries = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]
# tested_queries = ["P1"]
tested_query_types = ["VertexInduced", "EdgeInduced", "EdgeInducedIEP"]

@dataclasses.dataclass
class Runlog:
    is_available: bool = False
    is_oot: bool = False
    code_generation_time: float = 0.0
    code_execution_time: float = 0.0
    system: str = "None"
    graph: str = "None"
    query: str = "None"
    qtype: str = "None"
    number_runs: int = 0
    result: int = 0
    
def read_log(system, graph, query, qtype):
    log = Runlog()
    log.system = system
    log.graph = graph
    log.query = query
    log.qtype = qtype
    
    dir_path = os.path.dirname(__file__)
    file_path = os.path.join(dir_path, system, graph, query, qtype, "run.log")
    
    try:
        with open(file_path, "r") as file:
            for line in file:
                if "CODE_GENERATION_TIME(s)=" in line:
                    log.code_generation_time += float(line.split("CODE_GENERATION_TIME(s)=")[1].strip())
                if "COMPILATION_TIME(s)=" in line:
                    log.code_generation_time += float(line.split("COMPILATION_TIME(s)=")[1].strip())
                if "CODE_EXECUTION_TIME(s)=" in line:
                    log.code_execution_time += float(line.split("CODE_EXECUTION_TIME(s)=")[1].strip())
                if "RESULT=" in line:
                    log.result = int(line.split("RESULT=")[1].strip())
                    log.number_runs += 1
            
            if log.number_runs >= 1:
                log.code_execution_time /= log.number_runs
                log.code_generation_time /= log.number_runs
                log.is_available = True
                log.is_oot = False
            else:
                log.is_available = True
                log.is_oot = True

    except FileNotFoundError:
        # the file does not exist -> system does not support it
        log.is_available = False
    
    return log
        
def log2csv(logs: list[Runlog]):
    def find_log(system, graph, query, qtype):
        for log in logs:
            if log.system == system and log.graph == graph and log.query == query and log.qtype == qtype:
                return log
        return Runlog() # not available
    
    dir_path = os.path.dirname(__file__)
    file_path = os.path.join(dir_path, "log.csv")
    headers = ["Graph", "Query", "QueryType", \
        "GraphPi-Compilation(s)", "GraphPi-Execution(s)", "GraphPi-Runs", \
        "Dryadic-Compilation(s)", "Dryadic-Execution(s)", "Dryadic-Runs", \
        "MiniGraph-Compilation(s)", "MiniGraph-Code-Execution(s)", "MiniGraph-Runs"]
    
    with open(file_path, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        
        for graph in tested_graphs:
            for query in tested_queries:
                for qtype in tested_query_types:
                    cur_line = [graph, query, qtype]
                    
                    for system in tested_systems:
                        log = find_log(system, graph, query, qtype)
                        if log.is_available and not log.is_oot:
                            cur_line.append(str(log.code_generation_time))
                            cur_line.append(str(log.code_execution_time))
                            cur_line.append(str(log.result))
                            
                        elif log.is_available and log.is_oot:
                            cur_line.append("out of time")
                            cur_line.append("out of time")
                            cur_line.append("out of time")
                        else:
                            cur_line.append("not available")
                            cur_line.append("not available")
                            cur_line.append("not available")
                    
                    writer.writerow(cur_line)
                    
                    
if __name__ == "__main__":
    
    logs = []
    for system in tested_systems:
        for graph in tested_graphs:
            for query in tested_queries:
                for qtype in tested_query_types:
                    log = read_log(system, graph, query, qtype)
                    logs.append(log)
    
    log2csv(logs)
    