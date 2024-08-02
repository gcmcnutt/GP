#! /bin/sh

# ok we come in with a port -- let's adjust

PROGDIR=../../crsim/crrcsim-0.9.13
PROG=build/crrcsim
PORT=$1

cd $PROGDIR
$PROG -p autoc-config.xml -p "$PORT" -i AUTOC