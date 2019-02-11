# AutO
__An automatic orienteering event creator.__

Given a set of potential control locations, AutO selects a set of controls and creates a score orienteering course.  It produces files for upload to the MapRun server, which can be used in the MapRun app.

## Setup
Requires:
* A set of coordinates of potential controls (e.g. the sets of street light locations that are available from councils).
* A parameter file: a KML with a defined __limit__ polygon for the course, and a __start__ location.
* An OSRM server to calculate routes between points.
* An OOMap server to draw the map.
