# Bachelor Thesis of Samuel Holderbach

# What this is about
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