#! /bin/sh

# ok we come in with a port -- let's adjust

PROGDIR=../../crsim/crrcsim-0.9.13
PROG=build/crrcsim
INSTANCE=$1
PORT=$2
CRRCSIM_LOG=autoc_crrcsim.$$.log

# point to an instance of Xvfb for all but the first

INSTANCE=99

case "$INSTANCE" in
  "0")
    DISPLAY=:0
    ;;
  *)
    DISPLAY=:2
    ;;
esac

cd $PROGDIR
DISPLAY=$DISPLAY $PROG -g autoc_config.xml -p "$PORT" -i AUTOC > $CRRCSIM_LOG 2>&1