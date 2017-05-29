#!/bin/bash
mencoder "mf://*.png" -mf fps=4 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4
