#!/bin/bash
mencoder "mf://*.png" -mf fps=2 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4
