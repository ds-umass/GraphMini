# GraphMini
Artifacts for "GraphMini: Accelerating Graph Pattern Matching Using Auxiliary Graphs"

Github link: [https://github.com/ds-umass/GraphMini/](https://github.com/ds-umass/GraphMini/)

# Hardware Requirements
1. 128GB of free RAM to preprocess the graph Friendster correctly.
2. 180GB of free disk space to store preprocessed graphs.

# Software Requirements
### Operating System
1. Ubuntu 22.04 (tested)

### Libraries
1. CMake (Version >= 3.20)
2. Make
3. MPI
4. GCC compiler (>= 7)
5. Python (>= 3.8)
6. bc
7. clang-format (optional)

### Install libraries
```bash
sudo apt install bc cmake mpich clang-format -y
```

# Tested Graph Data
1. Wiki
2. YouTube
3. Patents
4. LiveJournal
5. Orkut
6. Friendster


All graphs are treated as unlabelled and undirected graphs.
We also provide scripts to download and preprocess these data automatically. See details below.

# Project Structure (including two baseline systems: GraphPi and Dryadic)
```markdown
- Datasets
    Subfolders:
    - Dryadic (preprocessed graph data for Dryadic)
    - GraphPi (preprocessed graph data for GraphPi)
    - GraphMini (preprocessed graph data for GraphMini)
    - TXT (downloaded graph data from SNAP)
    Scripts: the following scripts are used to download and preprocess datasets used in the experiments
    - download.sh (downloads data from SNAP)
    - dryadic_prep.sh scripts for preprocessing datasets into format can be handled by Dryadic
    - graphpi.sh scripts for preprocessing datasets into format can be handled by GraphPi
    - graphmini_prep.sh scripts for preprocessing datasets into format can be handled by GraphMini
    - prep.sh run all the scripts above
- Dryadic (source code/binaries of Dryadic)
- GraphPi (source code/binaries of GraphPi)
- GraphMini (source code/binaries of GraphMini)
- Experiments 
    Subfolders:
    - Dryadic (experiment results for Dryadic)
    - GraphMini (experiment results for GraphMini)
    - GraphPi (experiment results for GraphPi)
    Scripts: the following scripts are used to automate running the benchmark experiments that appeared in the paper
    - queries.sh (defines what queries and graphs to run in the experiment)
    - run_dryadic_edge.sh (test dryadic on edge-induced pattern defined in queries.sh)
    - run_dryadic_vertex.sh (test dryadic on vertex-induced pattern defined in queries.sh)
    - run_graphpi_edge_iep.sh (test graphpi on edge-induced pattern defined in queries.sh, with inclusion-exclusion optimization)
    - run_graphpi_edge.sh (test graphpi on edge-induced pattern defined in queries.sh, without inclusion-exclusion optimization)
    - run_graphmini_vertex.sh (test GraphMini on vertex-induced pattern defined in queries.sh, without inclusion-exclusion optimization)
    - run_graphmini_edge.sh (test GraphMini on edge-induced pattern defined in queries.sh, without inclusion-exclusion optimization)
    - run_graphmini_edge_iep.sh (test GraphMini on edge-induced pattern defined in queries.sh, with inclusion-exclusion optimization)
    - run_base_vertex.sh (test Base on vertex-induced pattern defined in queries.sh, without inclusion-exclusion optimization)
    - run_base_edge.sh (test Base on edge-induced pattern defined in queries.sh, without inclusion-exclusion optimization)
    - run_base_edge_iep.sh (test Base on edge-induced pattern defined in queries.sh, with inclusion-exclusion optimization)
```

# How to reproduce the experiment results
1. Build the source code for Dryadic, GraphPi and GraphMini:

```bash
bash ./build.sh
```

2. Download the dataset from SNAP and preprocess the datasets for Dryadic, GraphPi, and GraphMini:

```bash
bash ./Datasets/prep.sh
```

3. Define the queries to run and graphs to test by modifying the following file 

```bash
./Experiments/queries.sh
```

The `Queries` variable defines the queries to test. It takes inputs in (undirected) adjacency matrix format.

The `QuerySizes` variable defines the number of vertex in each of the query patterns. 

The `GraphNames` defines the data graphs to run the experiments. Examples are given in the script. Notice that the default setting will run experiments on all the graphs, which can take a significant amount of time, you can modify the `GraphNames` to run experiments on graphs of interest only.

The `TIMEOUT` defines the maximum amount of time before terminating a query. Notice that some queries can run for a very long time (>24h). You can modify the "TIMEOUT" variable to change the upper bound of query execution time. GraphMini can finish most queries in less than 12 hours on a 32-core system. 

4. Reproduce the benchmark experiments:

```bash
nohup bash ./Experiments/run.sh &
```

This script will run all the tested systems (GraphPi, Dryadic, GraphMini, Base) by invoking all the tested scripts in the directory (e.g. run_graphmini_vertex.sh). 
Each script can take more than 3 days to run. So if you want to get results quickly, you can use multiple machines and run them separately by invoking each script one by one.

The script will generate a log file for each tested query. The path to the log file is named as:
```
./Experiments/{System}_{Graph}_{Query}_{QueryType}/run.log
```

5. Generating Tables

To generate CSV tables corresponding to "Fig 5: GraphMini vs State-of-The-Art"

Fig 5a:

```bash
cd Experiments && python3 log2csv.py --baseline=Dryadic --target=GraphMini --adjtype=VertexInduced
```

Fig 5b:

```bash
cd Experiments && python3 log2csv.py --baseline=Dryadic --target=GraphMini --adjtype=EdgeInduced
```

Fig 5c:

```bash
cd Experiments && python3 log2csv.py --baseline=GraphPi --target=GraphMini --adjtype=EdgeInduced
```

Fig 5d:

```bash
cd Experiments && python3 log2csv.py --baseline=GraphPi --target=GraphMini --adjtype=EdgeInducedIEP
```

To generate tables corresponding to "Fig 6: GraphMini vs Base"
Fig 6a:

```bash
cd Experiments && python3 log2csv.py --baseline=Base --target=GraphMini --adjtype=EdgeInduced
```

Fig 6b:

```bash
cd Experiments && python3 log2csv.py --baseline=Base --target=GraphMini --adjtype=VertexInduced
```

# To learn more about GraphMini
See more detail in GraphMini/README.md
