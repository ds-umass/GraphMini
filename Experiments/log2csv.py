import os, csv, dataclasses, argparse

tested_systems = ["GraphPi", "Dryadic", "GraphMini", "Base"]
tested_graphs = ["wiki", "patents", "youtube", "lj", "orkut", "friendster"]
# tested_graphs = ["wiki"]
tested_queries = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]
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

def find_log(logs, system, graph, query, qtype) -> Runlog:
    for log in logs:
        if log.system == system and log.graph == graph and log.query == query and log.qtype == qtype:
            return log
    return Runlog() # not available
                    
def log2csv(baseline: str, target: str, qtype: str):
    
    def read_logs(system: str):
        logs = []
        for graph in tested_graphs:
            for query in tested_queries:
                log = read_log(system, graph, query, qtype)
                logs.append(log)
        return logs
    dir_path = os.path.dirname(__file__)
    file_path = os.path.join(dir_path, f"{baseline}_vs_{target}_{qtype}.csv")
    headers = ["Graph", "Query", "QueryType", \
        f"{baseline}-Compilation(s)", f"{baseline}-Execution(s)", f"{baseline}-Results", \
        f"{target}-Compilation(s)", f"{target}-Execution(s)", f"{target}-Results", \
        "Speed Up"]
    
    baseline_logs = read_logs(baseline)
    target_logs = read_logs(target)
    lines = []
    for graph in tested_graphs:
        for query in tested_queries:
            base_log = find_log(baseline_logs, baseline, graph, query, qtype)
            target_log = find_log(target_logs, target, graph, query, qtype)
            cur_line = [graph, query, qtype]
            for log in [base_log, target_log]:
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
            if base_log.is_available and target_log.is_available:
                cur_line.append(round(base_log.code_execution_time / target_log.code_execution_time, 3))
            elif not base_log.is_available and target_log.is_available:
                cur_line.append("INF")
            elif base_log.is_available and not target_log.is_available:
                cur_line.append("0")
                print('targe_log is not available')
            else:
                cur_line.append("out of time")
            lines.append(cur_line)
        
    with open(file_path, "w") as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        writer.writerows(lines)
    
    print(f"Result is writtern to {file_path}")
    
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='convert log file to speed up')
    parser.add_argument('--baseline', type=str, help='baseline system', choices=['GraphPi', 'GraphMini', "Dryadic", "Base"])
    parser.add_argument('--target', type=str, help='target system', choices=['GraphPi', 'GraphMini', "Dryadic", "Base"])
    parser.add_argument('--adjtype', default=str, type=str, help='Query Type', choices=["VertexInduced", "EdgeInduced", "EdgeInducedIEP"])
    args = parser.parse_args()
    
    log2csv(args.baseline, args.target, args.adjtype)    