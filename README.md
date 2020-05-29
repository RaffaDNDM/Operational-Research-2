# Operational-Research-2
### Algorithms of integer linear programming applied to TSP problem
- [CPLEX exact solvers.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/cplex_solver.h)
  - [Branch&Cut.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
  - Lazy callbacks.
    1. [Subtour Elimination Constraint (SEC).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
    2. [Patching (SEC).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
  - General lazy callbacks.
    1. [Subtour Elimination Constraint (SEC).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
    2. [Patching (SEC).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
- Compact solvers.
  - [Gavish Graves solver.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/gg_solver.h)
  - [MTZ solver.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/mtz_solver.h)
- Math-heuristic approaches
  - [Loop solver.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/loop_solver.h)
  - [Hard Fixing.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
  - [Soft fixing.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/bc_solver.h)
- [Heuristic solvers.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
  - Construction algorithm.
    1. [Nearest Neighborhood.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
    2. [Insertion.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
    3. [GASP variants of previous algorithms.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
  - Refinement algorithm.
    1. [Greedy Refinement.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
  - Meta-heuristich algorithm.
    1. [Variable Neighborhood Search (VNS).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
    2. [Tabu Search.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
    3. [Genetic algorithms.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
    4. [Simulating anealing.](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
    5. [Multi-start variations (using multithreading).](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/TSP/TSP/heuristic_solver.h)
### Report about used techniques
    All the previous implemented algorithms and the used CPLEX functions are explained in the [report](https://github.com/RaffaDNDM/Operational-Research-2/blob/master/Report/Report.pdf). In the same file you can find also the explanation of tsp instances format in tsplib dataset, used in training and test phases. In the report, there is also the explanation of Gnuplot tool, used to show result of algorithm and python programs used to keep input instances and generate performance profile of several algorithms.
