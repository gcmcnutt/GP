#! /bin/bash

# ok we come in with an instanceID and port -- let's adjust
INSTANCE=$1
PORT=$2

PROGDIR=${HOME}/crsim/crrcsim-0.9.13
PROG=build/crrcsim
CRRCSIM_LOG=/tmp/autoc_crrcsim.$$.log

# point to an instance of Xvfb for all but the first
# Xvfb :2 -screen 0 1024x768x8 &

case "$INSTANCE" in
  # "0")
  #   DISPLAY=:0
  #   ;;
  *)
    DISPLAY=:2
    ;;
esac

cd $PROGDIR
DISPLAY=$DISPLAY $PROG -g autoc_config.xml -p "$PORT" -i AUTOC > $CRRCSIM_LOG 2>&1