#!/usr/bin/env bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#echo "Script dir: $script_dir"
PROTO_DIR=$script_dir/../dtcc_core/proto

PYTHON_DIR=$script_dir/../dtcc_core/model
#
#
echo "Building Python classes..."
protoc --python_out=$PYTHON_DIR --proto_path=$PROTO_DIR dtcc.proto
#
#CPP_DIR=./src/cpp/protobuf
#echo "Building C++ classes..."
#protoc --cpp_out=$CPP_DIR --proto_path=$PROTO_DIR dtcc.proto
