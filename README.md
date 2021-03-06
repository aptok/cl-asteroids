
# Motivation

When I wanted to calculate the heliocentric position of asteroids from
the elements provided by the JPL, I startet writing these Lisp
functions.


# Functions

One can use `download-element-files-from-jpl` to download the orbital
elements from the Jet Propulsion Laboratory. They come in two files,
ELEMENTS.UNNUM and ELEMENTS.UNNUM, which contain more than one million
elements of asteroids.

Because of small differences in the two files and me to lazy to
abstract from them, there are two functions to read them into lists of
asteroid objects, `read-elements-numbered` and
`read-elements-unnumbered`.

Asteroid object are defined by `defclass asteroid`. The class is only
used to store the information contained in the element files. The `id`
of unnumbered asteroids is the same as there designation.

To calculate the heliocentric ecliptic position of an asteroid one can
use the function `asteroid-heliocentric-position-at-date` with an
asteroid object and a modified julian date.

