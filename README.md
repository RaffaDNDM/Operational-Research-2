# Operational-Research-2
### Algorithms of integer linear programming applied to TSP problem
-CPLEX
  - Exact solvers.
    - Compact solvers.
      1. [Gavish Graves solver.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/gg_solver.h)
      2. [MTZ solver.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/mtz_solver.h)
    - Non compact solvers.
      1. [Loop solver.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/loop_solver.h)
      2. [Lazy Constraint callback.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
      3. [Lazy Constraint callback + Heuristic callback (Patching).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
      4. [Generic callback.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
      5. [Generic callback (Patching).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
  - Math-heuristic solvers.
    1. [Hard Fixing.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
    2. [Soft fixing.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
- Heuristic solvers.
  - Construction algorithm.
    1. [Nearest Neighborhood.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)
    2. [Insertion.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)
    3. [GASP variants of previous algorithms.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)
  - Refinement algorithm.
    1. [Greedy Refinement (2-opt).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)
  - Meta-heuristich algorithm.
    1. [Multi-start variations (using multithreading).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)
    2. [Hybrid Variable Neighborhood Search (VNS).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)
    3. [Reactive Tabu Search.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)
    4. [Simulating anealing.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)
    5. [Genetic algorithms.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic.h)

### Report about used techniques

All the previous implemented algorithms and the used CPLEX functions are explained in the [report](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/Report/Report.pdf).
In the same file you can find also the explanation of tsp instances format in TSPlib dataset, used in training and test phases. In the report, there is also the explanation of Gnuplot tool, used to show results of algorithm and python programs used to keep input instances and generate performance profile of several algorithms.
