#!/usr/bin/env bash
# Post the live dataset check counts to Slack. Used by the live-dataset-checks
# workflow; reads NOTIFICATION_TITLE, SLACK_WEBHOOK_URL and the GitHub run
# variables from the environment.
set -euo pipefail

if [[ -z "${SLACK_WEBHOOK_URL:-}" ]]; then
  echo "No Slack webhook configured; skipping notification."
  exit 0
fi

COUNTS="$(python3 scripts/summarize_live_tests.py --counts tests/live-results.xml \
  || echo "results unavailable")"
RUN_URL="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${GITHUB_RUN_ID}"
export COUNTS RUN_URL

jq -n '{
  text: ([
    env.NOTIFICATION_TITLE,
    ("Repository: `" + env.GITHUB_REPOSITORY + "`"),
    ("Results: " + env.COUNTS),
    ("Run: " + env.RUN_URL)
  ] | join("\n"))
}' |
  curl --fail-with-body --silent --show-error \
    --header 'Content-Type: application/json' \
    --data-binary @- \
    "$SLACK_WEBHOOK_URL"
