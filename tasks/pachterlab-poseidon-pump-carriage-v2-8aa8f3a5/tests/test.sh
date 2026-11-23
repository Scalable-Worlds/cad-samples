#!/usr/bin/env bash
# Tests run as root. Tests must drop to UID 65534 before executing submitted code.
set -euo pipefail

# Set up the verifier.
PYTHON=/usr/local/bin/python3
OUT=/app/out TESTS=/tests LOGS=/logs/verifier
REWARD=/logs/verifier/reward.txt REPORT=/logs/verifier/ctrf.json
[ "$EUID" -eq 0 ]
rm -f "$REWARD" "$LOGS/reward.json" "$REPORT"
[ ! -L "$OUT" ]
mkdir -p /app/inputs "$OUT" "$LOGS"
[ "$(stat -c %u "$LOGS")" != 65534 ]

# Lock down the verifier.
chown 0:0 /app; chmod go-w /app
chown -hR 0:0 /app/inputs "$TESTS"
chmod -R a+rX,go-w /app/inputs
# Drop links and special files before root reads the submission.
find "$OUT" ! -type d ! -type f -delete
chown -hR 65534:65534 "$OUT"
chmod 700 "$TESTS" "$LOGS"

# Run pytest.
cd "$TESTS"
"$PYTHON" -P -c 'import pytest, ctrf'
status=0
PYTHONDONTWRITEBYTECODE=1 "$PYTHON" -P -m pytest -p no:cacheprovider -rA \
    --ctrf "$REPORT" "$TESTS/test_outputs.py" || status=$?

# Write the reward.
rm -f "$REWARD" "$LOGS/reward.json"
case "$status" in
    0) echo 1 > "$REWARD" ;;
    2|3|4|5) echo "pytest exit $status: the tests did not run, no reward" ;;
    *) echo 0 > "$REWARD" ;;
esac
exit "$status"
