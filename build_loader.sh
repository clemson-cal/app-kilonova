#!/bin/bash
target_dir=knc_loader/target/release
cd knc_loader; cargo build --release; cd ..
rslib=$(ls $target_dir | egrep '.so|.dylib')
pylib=$(pwd)/lib/knc_loader.so
mkdir -p lib
rm -f $pylib
ln -s $(pwd)/$target_dir/$rslib $pylib
pypath=$(pwd)/lib
if [ -d "$pypath" ] && [[ ":$PYTHONPATH:" != *":$pypath:"* ]]; then
    export PYTHONPATH="${PYTHONPATH:+"$PYTHONPATH:"}$pypath"
fi
echo "link $(pwd)/$target_dir/$rslib -> $pylib"
echo "add knc_loader to your Python path: PYTHONPATH=$PYTHONPATH"
unset target_dir
unset rslib
unset pylib
unset pypath