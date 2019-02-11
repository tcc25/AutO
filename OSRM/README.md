AutO requires an OSRM server to be set up.  

To set up, use the [osrm-backend](https://hub.docker.com/r/osrm/osrm-backend/) docker container (v5.20.0 works), with modified_foot.lua

This makes the modifications to regular foot.lua:
* Allow travel along 'trunk' roads.
* Allow travel along cycle paths.
