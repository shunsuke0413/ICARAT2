#!/bin/bash

find ../ -type f -name "*.vtk" -exec rm -f {} \;
find ../ -type f -name "*.vtu" -exec rm -f {} \;
find ../ -type f -name "*.pvtu" -exec rm -f {} \;
find ../ -type f -name "*.out" -exec rm -f {} \;
find ../ -type d -name "objs" -exec rm -rf {} \;