#!/bin/sh

# The first command compiles the model, so is not included in timing
libbi sample --model-file sir.bi --target prediction --end-time 100000 --nparticles 10
/usr/bin/time libbi sample --model-file sir.bi --target prediction --end-time 100000 --nparticles 10
