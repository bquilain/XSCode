#!/bin/bash

for f in `ls *.so`; do ln -s $f lib$f; done
