#!/bin/bash

xargs -L 1 curl -O -J -L < ENCODE_files.txt
