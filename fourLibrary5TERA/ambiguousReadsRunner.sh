#!/bin/bash

# I want to have this script run the ambiguousReads.py script on a short list of libraries:
for LIBRARY in newerN2 newerS5 newerS6 oldN2 oldS6 thirdN2 thirdS6 thirdS5 temp25cN2 temp25cS5 temp25cS6 temp25cS7
do
    python ambiguousReads.py $LIBRARY input
done
