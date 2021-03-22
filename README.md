# RapidMatch
## Introduction
A subgraph query searches for all embeddings in a data
graph that are identical to a query graph. Two kinds of
algorithms, either graph exploration based or join based,
have been developed for processing subgraph queries. Due
to algorithmic and implementational differences, join-based
systems can handle query graphs of a few vertices efficiently
whereas exploration-based approaches typically process up
to several tens of vertices in the query graph. In this paper,
we first compare these two kinds of methods and prove
that the complexity of result enumeration in state-of-the-art
exploration-based methods matches that of the worst-case optimal join.
Furthermore, we propose RapidMatch, a holistic subgraph query processing
framework integrating the two approaches. Specifically, RapidMatch not
only runs relational operators such as selections and joins, but also
utilizes graph structural information, as in graph exploration,
for filtering and join plan generation. Consequently, it outperforms the
state of the art in both approaches on a wide range of query workloads.

For the details, please refer to our VLDB'2021 paper
"RapidMatch: A Holistic Approach to Subgraph Query Processing [[Preview](https://www.comp.nus.edu.sg/~hebs/pub/rapidmatch-vldb21.pdf)]"
by [Shixuan Sun](https://github.com/shixuansun), [Xibo Sun](https://github.com/xibosun), [Yulin Che](https://github.com/CheYulin),
[Prof. Qiong Luo](http://www.cse.ust.hk/~luo/), and [Prof. Bingsheng He](https://www.comp.nus.edu.sg/~hebs/).
If you have any further questions, please feel free to contact us.

Please cite our paper, if you use our source code.

* "Shixuan Sun, Xibo Sun, Yulin Che, Qiong Luo, and Bingsheng He. RapidMatch: A Holistic Approach to Subgraph
Query Processing. VLDB 2021."


## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

## Test
Execute the following commands to test the correctness of the binary file.

```zsh
cd test
python test.py ../build/matching/RapidMatch.out
```

## Execute
After compiling the source code, you can find the binary file 'RapidMatch.out'
under the 'build/matching' directory.  Execute the binary with the following
command './RapidMatch.out -d data_graphs -q query_graphs -order nd
-preprocess true -num number_of_embeddings -time_limit time_in_seconds',
in which '-d' specifies the input of the data graphs and '-q' specifies the
input of the query graphs. The '-order' parameter gives the ordering method, which
is 'nd'. 'nd' denotes the join plan generation method based on the nucleus
decomposition. Set '-preprocess' as 'true' to enable the filtering method based
on semi-join operations. The '-num' parameter sets the maximum number of
embeddings that you would like to find. If the number of embeddings enumerated
reaches the limit or all the results have been found, then the program will terminate.
Set '-num' as 'MAX' to find all results. The '-time_limit' parameter configures the
time budget for the query. If the query cannot be completed within the time limit,
then the program will terminate the query and return the number of results found.

Example (Execute the query with the filtering method enabled, and find all results.
The time limit is 60 seconds.):

```zsh
./RapidMatch.out -d ../../dataset/simple_dataset/test_case_1.graph -q ../../dataset/simple_dataset/query1_positive.graph -order nd -preprocess true -num MAX -time_limit 60
```

## Input
Both the input query graph and data graph are vertex-labeled.
Each graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted
as 'v VertexID LabelId Degree' and 'e VertexId VertexId' respectively. Note that we require that the vertex
id is started from 0 and the range is [0,N - 1] where V is the vertex set. The following
is an input sample. You can also find sample data sets and query sets under the test folder.

Example:

```zsh
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

## Configuration

You can configure the data layout (Encoded Trie, Hash Table, Trie),
set intersection algorithms (Merge, Hybrid, Merge+AVX2, Hybrid+AVX2, QFilter),
optimization techniques (Intersection Caching, Failing Set Pruning)
and result types (Homomorphism, Isomorphism) by defining macros in
'configuration/config.h'.

| Macro | Description |
| :-----------------------------------------------: | :-------------: |
|HYBRID 0| a hybrid method handling the cardinality skew by integrating the merge-based method with the galloping-based method  |
|HYBRID 1| the merge-based set intersection |
|SI 0 | Accelerate the set intersection with AVX2 |
|SI 1 | Accelerate the set intersection with AVX512 |
|SI 2 | Scalar set intersection |
|RELATION_STRUCTURE 0 | Encoded Trie | 
|RELATION_STRUCTURE 1 | Hash Table |
|RELATION_STRUCTURE 2 | Trie |
|SPARSE_BITMAP| Enable the QFilter set intersection method|
|FAILING_SET_PRUNING| Enable the failing set pruning method|
|INTERSECTION_CACHE| Enable the intersection caching method|
|HOMOMORPHISM| Find the subgraph homomorphisms|

In our paper, we execute the large queries with the following configuration, which is the default setting. We set
the time limit as 300 seconds (5 minutes) and the number of embeddings as 100000.

| Macro | Description |
| :-----------------------------------------------: |  :-------------: |
|HYBRID 0| Hybrid|
|SI 0 | AVX2 |
|RELATION_STRUCTURE 0 | Encoded Trie | 
|FAILING_SET_PRUNING| Enable the failing set pruning method|

We execute the small queries with the following configuration. We set the time limit as 86400 seconds (24 hours) and the number
of embeddings as MAX. **Note that when finding homomorphisms, you need to disable the failing set pruning technique by commenting out "FAILING_SET_PRUNING" in config.h because this optimization is based on the definition of isomorphism.**


| Macro | Description |
| :-----------------------------------------------: |  :-------------: |
|HYBRID 0| Hybrid|
|SI 0 | AVX2 |
|RELATION_STRUCTURE 0 | Encoded Trie |
|SPARSE_BITMAP| Enable the QFilter set intersection method|
|INTERSECTION_CACHE| Enable the intersection caching method|
|HOMOMORPHISM| Find the subgraph homomorphisms|

## Experiment Datasets
The real world datasets and the corresponding query sets used in our paper can be downloaded [here](https://hkustconnect-my.sharepoint.com/:u:/g/personal/ssunah_connect_ust_hk/EWwS7ixh4NBHriiPHNpUMAkBu8vbH1f37Ug8CPWQdUXj4w?e=0GXEMg).

## Outside Code
Our project utilizes some outside source code, which is listed in the following.

| Description | GitHub Link |
| :-----------------------------------------------: |  :-------------: |
|Hash Join | https://github.com/wagjamin/HashJoins|
|QFilter | https://github.com/pkumod/GraphSetIntersection |
|Sparsepp | https://github.com/greg7mdp/sparsepp |
|Nucleus Decomposition | https://github.com/sariyuce/nucleus|
