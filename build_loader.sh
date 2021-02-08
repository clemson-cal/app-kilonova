#!/bin/bash

cd knc_loader
maturin build --release
pip3 install --force-reinstall target/wheels/*.whl
