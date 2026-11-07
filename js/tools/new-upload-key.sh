#!/bin/sh
# Installs a fresh upload key for the catalogue page's remote results on the results host
# (over ssh, so only a machine that can ssh there can mint one) and copies it to the
# clipboard. The previous two keys stay valid, so jobs still running with them keep working.
# Usage: tools/new-upload-key.sh [ssh-host]   (default dcp.network)
set -eu
HOST=${1:-dcp.network}
KEY=$(python3 -c 'import secrets; print(secrets.token_urlsafe(32))')
printf '%s\n' "$KEY" | ssh -o BatchMode=yes "$HOST" \
  'umask 077; f=~/3kidnas-results/token; { cat; [ -f "$f" ] && head -n 2 "$f"; } > "$f.new" && mv "$f.new" "$f"'
printf '%s' "$KEY" | pbcopy
echo "New upload key installed on $HOST and copied to the clipboard. Paste it into the page's Upload key field."
