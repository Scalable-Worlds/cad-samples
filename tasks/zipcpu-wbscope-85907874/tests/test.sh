#!/usr/bin/env bash
# Tests run as root. Programs they launch run as UID 65534: they can read the
# inputs and their own submission, but not the tests or the reward.
set -euo pipefail
PYTHON=/usr/local/bin/python3
OUT=/app/out TESTS=/tests LOGS=/logs/verifier REPORT=/logs/verifier/ctrf.json

[ "$EUID" -eq 0 ]
rm -f /logs/verifier/reward.txt /logs/verifier/reward.json "$REPORT"
[ ! -L "$OUT" ]
mkdir -p /app/inputs "$OUT" "$LOGS"
[ "$(stat -c %u "$LOGS")" != 65534 ]

chown 0:0 /app; chmod go-w /app
chown -hR 0:0 /app/inputs "$TESTS"; chmod -R a+rX,go-w /app/inputs
# Drop links from the submission: root would read a protected file through one.
find "$OUT" -type l -delete
chown -hR 65534:65534 "$OUT"
chmod 700 "$TESTS" "$LOGS"

cd "$TESTS"
"$PYTHON" -P -c 'import pytest, ctrf'
status=0
PYTHONDONTWRITEBYTECODE=1 "$PYTHON" -P -m pytest -p no:cacheprovider -rA \
    --ctrf "$REPORT" "$TESTS/test_outputs.py" || status=$?

rm -f /logs/verifier/reward.txt /logs/verifier/reward.json
case "$status" in
    0) echo 1 > /logs/verifier/reward.txt ;;
    2|3|4|5) echo "pytest exit $status: the tests did not run, no reward" ;;
    *) echo 0 > /logs/verifier/reward.txt ;;
esac
exit "$status"
