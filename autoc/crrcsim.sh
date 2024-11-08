#! /bin/bash

# ok we come in with an instanceID and port -- let's adjust
BASE=$1
INSTANCE=$2
PORT=$3

PROGDIR=${HOME}/crsim/crrcsim-0.9.13
PROG=build/crrcsim
CRRCSIM_LOGDIR=/tmp/crrcsim
CRRCSIM_LOG=${CRRCSIM_LOGDIR}/autoc_crrcsim-$1.$$.log

# point to an instance of Xvfb for all but the first
# Xvfb :2 -screen 0 1024x768x8 &
mkdir -p "${CRRCSIM_LOGDIR}"

case "$INSTANCE" in
  # "0")
  #   DISPLAY=:0
  #   ;;
  *)
    DISPLAY=:0
    ;;
esac

cd $PROGDIR
DISPLAY=$DISPLAY stdbuf -o0 -e0 $PROG -g autoc_config.xml -p "$PORT" -i AUTOC > $CRRCSIM_LOG 2>&1