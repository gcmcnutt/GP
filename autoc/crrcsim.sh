#! /bin/sh

# ok we come in with a port -- let's adjust

PROGDIR=../../crsim/crrcsim-0.9.13
PROG=build/crrcsim
PORT=$1
CRRCSIM_LOG=autoc_crrcsim.$$.log

# point to an instance of Xvfb
export DISPLAY=:2

cd $PROGDIR
$PROG -g autoc_config.xml -p "$PORT" -i AUTOC > $CRRCSIM_LOG 2>&1