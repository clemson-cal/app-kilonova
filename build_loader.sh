#!/bin/bash

pathadd() {
    if [ -d "$1" ] && [[ ":$PYTHONPATH:" != *":$1:"* ]]; then
        export PYTHONPATH="${PYTHONPATH:+"$PYTHONPATH:"}$1"
    fi
}

target_dir=knc_loader/target/release
cd knc_loader; cargo build --release; cd ..

rslib=$(ls $target_dir | egrep '.so|.dylib')
pylib=$(pwd)/lib/knc_loader.so
mkdir -p lib
rm -f $pylib
ln -s $(pwd)/$target_dir/$lib $pylib
pathadd $(pwd)/lib

echo "adding knc_loader to your Python path:"
echo "PYTHONPATH=$PYTHONPATH"
