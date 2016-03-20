#!/bin/sh

code_formatter="clang-format"

if [[ $# -eq 0 ]] ; then
    echo "Need to provide a source directory"
    exit 1
fi

if [ ! -d "$1" ]; then
    echo "No such directory:  $1"
    exit 1
fi

if ! [ -x "$(command -v $code_formatter)" ]; then
    echo "$code_formatter not found; `basename "$0"` requires $code_formatter to be installed!" >&2
    exit 1
fi

find $1 -name '*.c' -o -name '*.cpp' -o -name '*.h' -o -name '*.hpp' -o -name '*.cc' -o -name '*.cu' | while read -r source_file; do
    echo "running:  $code_formatter -i $source_file";
    # $code_formatter -i $source_file;
done