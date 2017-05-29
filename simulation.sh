#!/bin/bash
rm img/*.png
octave printInputs.m
./m2slv01
octave printOutImg.m
