#!/bin/bash
# Script to analyze DTEST logs from crrcsim workers to find divergence points

LOGDIR=/tmp/crrcsim

echo "=== DTEST Log Analysis ==="
echo ""

# Take base process ID as argument, or find most recent
if [ -n "$1" ]; then
    LATEST_BASE=$1
    echo "Using specified base process ID: $LATEST_BASE"
else
    # Find most recent autoc run by looking at log timestamps
    LATEST_BASE=$(ls -t ${LOGDIR}/autoc_crrcsim-*.log 2>/dev/null | head -1 | sed 's/.*autoc_crrcsim-\([0-9]*\)\..*/\1/')
    if [ -z "$LATEST_BASE" ]; then
        echo "No crrcsim logs found in $LOGDIR"
        echo "Usage: $0 [base_process_id]"
        exit 1
    fi
    echo "Auto-detected most recent base process ID: $LATEST_BASE"
fi

echo "Analyzing logs from run: $LATEST_BASE"
echo ""

# Count how many worker logs we have
WORKER_LOGS=$(ls ${LOGDIR}/autoc_crrcsim-${LATEST_BASE}.*.log 2>/dev/null)
NUM_WORKERS=$(echo "$WORKER_LOGS" | wc -l)
echo "Found $NUM_WORKERS worker logs"
echo ""

# Extract DTEST_RESET entries from each worker
echo "=== DTEST_RESET Comparison (First evaluation) ==="
for log in $WORKER_LOGS; do
    PID=$(basename $log | sed 's/.*\.\([0-9]*\)\.log/\1/')
    FIRST_RESET=$(grep "DTEST_RESET" $log 2>/dev/null | head -1)
    if [ -n "$FIRST_RESET" ]; then
        echo "Worker PID=$PID: $FIRST_RESET"
    fi
done
echo ""

# Extract DTEST_FIRST entries from each worker
echo "=== DTEST_FIRST Comparison (First timestep after reset) ==="
for log in $WORKER_LOGS; do
    PID=$(basename $log | sed 's/.*\.\([0-9]*\)\.log/\1/')
    FIRST_STEP=$(grep "DTEST_FIRST" $log 2>/dev/null | head -1)
    if [ -n "$FIRST_STEP" ]; then
        echo "Worker PID=$PID: $FIRST_STEP"
    fi
done
echo ""

# Check if there are any differences in the position/quaternion values
echo "=== Checking for differences in initial state ==="
echo "Extracting position coordinates from DTEST_RESET..."
grep "DTEST_RESET" ${LOGDIR}/autoc_crrcsim-${LATEST_BASE}.*.log 2>/dev/null | \
    head -${NUM_WORKERS} | \
    sed 's/.*pos\[\([^]]*\)\].*/\1/' | \
    sort | uniq -c
echo ""

echo "Extracting quaternion from DTEST_RESET..."
grep "DTEST_RESET" ${LOGDIR}/autoc_crrcsim-${LATEST_BASE}.*.log 2>/dev/null | \
    head -${NUM_WORKERS} | \
    sed 's/.*quat\[\([^]]*\)\].*/\1/' | \
    sort | uniq -c
echo ""

echo "Extracting velocity from DTEST_RESET..."
grep "DTEST_RESET" ${LOGDIR}/autoc_crrcsim-${LATEST_BASE}.*.log 2>/dev/null | \
    head -${NUM_WORKERS} | \
    sed 's/.*vel\[\([^]]*\)\].*/\1/' | \
    sort | uniq -c
echo ""

# Check for DTEST_FIRST entries
echo "=== DTEST_FIRST Analysis ==="
FIRST_COUNT=$(grep -h "DTEST_FIRST" ${LOGDIR}/autoc_crrcsim-${LATEST_BASE}.*.log 2>/dev/null | wc -l)
if [ "$FIRST_COUNT" -gt 0 ]; then
    echo "Found $FIRST_COUNT DTEST_FIRST entries"
    echo ""
    echo "First timestep position comparison:"
    grep "DTEST_FIRST" ${LOGDIR}/autoc_crrcsim-${LATEST_BASE}.*.log 2>/dev/null | \
        head -${NUM_WORKERS} | \
        sed 's/.*pos\[\([^]]*\)\].*/\1/' | \
        sort | uniq -c
    echo ""
    echo "First timestep quaternion comparison:"
    grep "DTEST_FIRST" ${LOGDIR}/autoc_crrcsim-${LATEST_BASE}.*.log 2>/dev/null | \
        head -${NUM_WORKERS} | \
        sed 's/.*quat\[\([^]]*\)\].*/\1/' | \
        sort | uniq -c
    echo ""
    echo "First timestep time comparison:"
    grep "DTEST_FIRST" ${LOGDIR}/autoc_crrcsim-${LATEST_BASE}.*.log 2>/dev/null | \
        head -${NUM_WORKERS} | \
        sed 's/.*time=\([0-9]*\).*/\1/' | \
        sort | uniq -c
else
    echo "WARNING: No DTEST_FIRST entries found!"
    echo "This may indicate the first timestep logging is not working."
fi
echo ""

echo "=== Analysis complete ==="
