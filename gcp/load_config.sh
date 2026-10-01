# Load gcp/config.env, if it exists, as *defaults*.
#
# Sourced by fjordsim-gcp and build_image.sh. It reads the file rather than `source`-ing it so a
# value already in the environment wins, which is what makes the documented one-off work:
#
#   BUCKET=gs://other-bucket gcp/fjordsim-gcp run --config oslofjorden --steps run_simulation
#
# Plain `source` would assign unconditionally and silently ignore that override.

_fjordsim_load_config() {
    local file="$1" line name value
    [[ -f "$file" ]] || return 0
    while IFS= read -r line || [[ -n "$line" ]]; do
        # Anything that is not NAME=... — comments and blank lines — is skipped.
        [[ "$line" =~ ^[[:space:]]*([A-Za-z_][A-Za-z0-9_]*)=(.*)$ ]] || continue
        name="${BASH_REMATCH[1]}"
        [[ -n "${!name:-}" ]] && continue
        value="${BASH_REMATCH[2]}"
        value="${value%\"}"; value="${value#\"}"
        value="${value%\'}"; value="${value#\'}"
        printf -v "$name" '%s' "$value"
        export "$name"
    done < "$file"
}
