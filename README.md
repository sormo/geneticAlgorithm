### Genetic Algorithm
Implementation of genetic algorithm as a generic function that is applied on different examples. Json parsing is done with [this](https://github.com/nlohmann/json) parser. Visualization with [this](https://github.com/adishavit/simple-svg) header for svg drawing.   

#### 0/1 Knapsack
Try to get optimal value for Knapsack problem. Problem is defined as an array of pairs where pair represents item's weight and value.  
Gene is single bit. Position of bit (gene) in bit string is index of particular item. If bit is set, item is taken into backapack. Fitness of chromosome is sum of values of all items until maximum weight is reached. Mutation is done by toggling genes.  

#### Traveling Salesman
Try to get optimal vaue for Traveling Salesman problem. Problem is defined as an array of points (x, y value) which represents cities that must be visited.  
In this case gene is unsigned char. Chromosome is ordered array of indices into array of cities. Fitness is computed by summing distances between cities in order in which they appear in chromosome. Mutation is done by swapping genes.

#### Digits and Operators
Example taken from [here](http://www.ai-junkie.com/ga/intro/gat3.html).
> Given the digits 0 through 9 and the operators +, -, * and /,  find a sequence that will represent a given target number. The operators will be applied sequentially from left to right as you read.  

Chromosome is represented as bit string. Each 4 bits is mapped into one of allowed characters: numbers 0-9 and operators: +-*/ . Fitness is computed by converting chromosome into expression (while ingoring invalid characters). Value of fitness is inversed distance between target number and value of expression. Mutation is done by toggling genes.

#### Circles
This example is also taken from [here](http://www.ai-junkie.com/ga/intro/gat3.html).
Problem is defined as an array of circles. Find circle with maximum radius that does not cross another circle or frame.  
Gene in this case is double. Chromosome represents single circle as 3 doubles (x, y coordinates and radius). Fitness is either radius of circle if circle is valid or zero if it's not. Mutation is done by adding random value (in some range) to current value of gene.

#### Rectangles
Try to find path between randomly generated rectangles from bottom up. Gene is double and chromosome represents points of path (pairs of doubles represents single point). Fitness is height of last point of path.  
This example usually gets stuck in local maximum and is not able to find path even there is route.
