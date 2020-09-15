# Bachelor Thesis of Samuel Holderbach

## What this is about
A simple graph puzzle called <i>Minorensucher</i> or <i>minorfinder</i>.
The game begins by showing the player a bigger graph (the original graph) and a smaller one (the minor). To solve the puzzle the player has to contract eges untill the original graph is isomorphic to the minor.

## Build the game
### Generate Unix Makefile
- Install any gtk version (2 or 3)
- Run the perl script `mkfiles.pl`
- Run the batch scipt `mkauto.sh`
- Run `./configure`
### Actual build and subsequent builds
- Finally run `make`

## Run the game
- Run `./minorfinder`

## Play the game
### Possible moves
There are three possible moves to solve the puzzle:
1. Contract edges (Merge the incident vertices of an edge and preserve their incident edges)
2. Delete edges
3. Delete vertices with no incident edges
### Control
1. Contract edges: **LEFT MOUSE BUTTON**,
  click on the edge that you want to contract
2. Delete edges: **RIGHT MOUSE BUTTON**,
  click on the edge that you want to delete
3. Delete vertices with no incident edges: **RIGHT MOUSE BUTTON**,
  click on the vertex that you want to delete
### Goal
Repeatedly apply one of the possible moves to the bigger graph until it appears to be isomorphic to the minor.
You will notice a quick flash when you've found a solution. If you can't find a solution you can use the solve
function which will try to find a solution for you.
