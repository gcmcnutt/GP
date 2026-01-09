#! /bin/bash

# ok we come in with an instanceID and port -- let's adjust
BASE=$1
INSTANCE=$2
PORT=$3

PROGDIR=${HOME}/crsim/crrcsim-0.9.13
PROG=build/crrcsim
CRRCSIM_LOGDIR=/tmp/crrcsim
# Include worker index and pid in the log name for easier correlation
CRRCSIM_LOG=${CRRCSIM_LOGDIR}/autoc_crrcsim-${BASE}.${INSTANCE}.$$.log
export AUTOC_WORKER_INDEX=${INSTANCE}

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


DISPLAY=$DISPLAY stdbuf -o0 -e0 $VALGRIND_CMD $PROG -g autoc_config.xml -p "$PORT" -i AUTOC > $CRRCSIM_LOG 2>&1