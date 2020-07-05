# -*- makefile -*-

# The `nullgame' source file is a largely blank one, which contains
# all the correct function definitions to compile and link, but
# which defines the null game in which nothing is ever drawn and
# there are no valid moves. Its main purpose is to act as a
# template for writing new game definition source files. I include
# it in the Makefile because it will be worse than useless if it
# ever fails to compile, so it's important that it should actually
# be built on a regular basis.
minorfinder : [X] GTK COMMON minorfinder minorfinder-icon|no-icon
minorfinder : [G] WINDOWS COMMON minorfinder minorfinder.res|noicon.res