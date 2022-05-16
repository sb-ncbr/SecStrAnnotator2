#!/bin/bash

# Usage:
#     bash build.sh
#     bash build.sh --release

set -e 

DIR=$(realpath $(dirname $0))

echo $DIR


CONFIG="Debug";

for ARG in $@; do
    if [ "$ARG" = "--release" ]; then
        CONFIG="Release";
    fi;
done;

date --utc +'%F %T' > $DIR/buildtime
dotnet build $DIR/SecStrAnnotator.csproj -c $CONFIG
cp $DIR/SecStrAnnotator_config.json $DIR/bin/$CONFIG/net6.0/SecStrAnnotator_config.json;
