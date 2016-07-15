#!/usr/bin/env bash
git submodule update --init

cd cgen
git am ../patches/0001-ctye-to-dtype-handle-void.patch
cd ..

python setup.py install --user
