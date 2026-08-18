#!/usr/bin/env bash
set -euo pipefail

# The installed command owns defaults, validation, resume markers and provenance.
exec nucleosuite mnase-suite "$@"
