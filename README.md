# GraphMini-Artifacts
Artifacts for GraphMini (MiniGraph was its initial name)

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

### Install libraries
```bash
sudo apt install cmake mpich -y
```

# Tested Graph Data
1. Wiki
2. YouTube
3. Patents
4. LiveJournal
5. Orkut
6. Friendster


All graphs are treated as unlabelled and undirected graph.
We also provide scripts to download and preprocess these data automatically. See details below.

# Project Structure (including two baseline systems: GraphPi and Dryadic)
```markdown
- Datasets
    Subfolders:
    - Dryadic (preprocessd graph data for Dryadic)
    - GraphPi (preprocessd graph data for GraphPi)
    - MiniGraph (preprocessd graph data for MiniGraph)
    - TXT (downloaded graph data from SNAP)
    Scripts: the following scripts are used to download and preprocess datasets used in the experiments
    - download.sh (downloads data from SNAP)
    - dryadic_prep.sh scripts for preprocessing datasets into format can be handled by Dryadic
    - graphpi.sh scripts for preprocessing datasets into format can be handled by GraphPi
    - minigraph_prep.sh scripts for preprocessing datasets into format can be handled by MiniGraph
    - prep.sh run all the scripts above
- Dryadic (source code / binaries of Dryadic)
- GraphPi (source code / binaries of GraphPi)
- MiniGraph (source code / binaries of MiniGraph)
- Experiments 
    Subfolders:
    - Dryadic (experiment results for Dryadic)
    - MiniGraph (experiment results for MiniGraph)
    - GraphPi (experiment results for GraphPi)
    Scripts: the following scripts are used to automate running the benchmark experiments appeared in the paper
    - queries.sh (defines what queries and graphs to run in the experiment)
    - run_dryadic_edge.sh (test dryadic on edge-induced pattern defined in queries.sh)
    - run_dryadic_vertex.sh (test dryadic on vertex-induced pattern defined in queries.sh)
    - run_graphpi_edge_iep.sh (test graphpi on edge-induced pattern defined in queries.sh, with inclusion-exclusion optimization)
    - run_graphpi_edge.sh (test graphpi on edge-induced pattern defined in queries.sh, without inclusion-exclusion optimization)
    - run_minigraph_vertex.sh (test minigraph on vertex-induced pattern defined in queries.sh, without inclusion-exclusion optimization)
    - run_minigraph_edge.sh (test minigraph on edge-induced pattern defined in queries.sh, without inclusion-exclusion optimization)
    - run_minigraph_edge_iep.sh (test minigraph on edge-induced pattern defined in queries.sh, with inclusion-exclusion optimization)
```

# How to reproduce the experiment results
1. Build the source code for Dryadic, GraphPi and GraphMini:

```bash
bash ./build.sh
```

2. Download dataset from snap and preprocess the datasets for Dryadic, GraphPi and GraphMini:

```bash
bash ./Datasets/prep.sh
```

3. Define the queries to run and graphs to test by modifying the following file 

```bash
./Experiments/queries.sh
```

The "Queries" variable defines the queries to test. It takes inputs in (undirected) adjacency matrix format
The "QuerySizes" variable defines the number of vertex in the each of the query pattern. 
The "GraphNames" defines the data graphs to run the experiments. Examples are given in the script.

4. Reproduce the benchmark experiments (ex. vertex-induced) of GraphMini:

```bash
nohup bash ./Experiments/run_minigraph_vertex.sh &
```
Similarly, you can use other scripts in the `Experiment` folder to reproduce the results for different baseline systems and settings.

# Convert log files into csv output.

Each query will produce a run.log file to record the query compilation time and query execution time.
You can easily convert the results using the script `./Experiments/log2csv.py` to convert them into a csv file called `log.csv`.

# How to run a single query with GraphMini
1. You must first preprocess the dataset and build the binaries of GraphMini as described above.

The binary takes 5 inputs to execute:
```bash
./MiniGraph/build/bin/run [graph_nickname] [path_to_graph] [query_nickname] [query_adjmat] [query_type]
```
- graph_nickname: nickname for tested graph (ex. wiki)
- path_to_graph: path to the directory that contains the preprocessed graph. (ex ./Datasets/MiniGraph/wiki)
- query_nickname: nickname for tested query (ex. P1)
- query_adjmat: adjacency matrix of the tested query (ex. "0111101111011110" 4-clique)
- query_type: 
    - 0: vertex-induced, 
    - 1: edge-induced, 
    - 2: edge-induced with IEP optimization

An example query:

```bash
./MiniGraph/build/bin/run wiki ./Datasets/MiniGraph/wiki P1 0111101111011110 0
```

# How to run GraphMini with different configuration
The default setting uses a cost model to determine which vertices to prune and uses nested parallelism to speed up execution.
In `MiniGraph/src/run.cpp`, you can modify the following two lines to change the configuration. You also need to rebuild the binaries to ensure the changes are made.

```cpp
line 210: conf.pruningType = PruningType::CostModel;
line 211: conf.parType     = ParallelType::NestedRt;
```
Acceptable configuration for `conf.pruningType` includes:
```cpp
enum class PruningType {
    None = 0, // not using minigraph
    Static = 1, // use if and only if not introducing any redundant set operation
    Eager = 2, // use for all potential cases
    Online = 3, // use iff pruned adj will be used once
    CostModel = 4, // lazy + cost model
};

enum class ParallelType {
    OpenMP,
    TbbTop, // only parallel at top level like OpenMP
    Nested, // aggressively nested
    NestedRt, // decide whether to parallel based on cost model
};
```
    
