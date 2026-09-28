#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<EOF
Usage: ${0##*/} -i FILE -m INT -M INT -n INT

Extract random, potentially overlapping strands (substrings) from a text file.
Newlines in the input file are ignored, treating the text as contiguous.

Options:
  -i, --input FILE    Path to the input text file
  -m, --min INT       Minimum length of a strand
  -M, --max INT       Maximum length of a strand
  -n, --count INT     Number of strands to produce
  -h, --help          Display this help message
EOF
}

validate_args() {
  local -r file="$1"
  local -r min="$2"
  local -r max="$3"
  local -r count="$4"

  if [[ -z "${file}" || -z "${min}" || -z "${max}" || -z "${count}" ]]; then
    echo "Error: Missing required arguments." >&2
    usage >&2
    exit 1
  fi

  if [[ ! -f "${file}" ]]; then
    echo "Error: Input file '${file}' not found." >&2
    exit 1
  fi

  if [[ ! -r "${file}" ]]; then
    echo "Error: Input file '${file}' is not readable." >&2
    exit 1
  fi

  if ! [[ "${min}" =~ ^[1-9][0-9]*$ ]] \
      || ! [[ "${max}" =~ ^[1-9][0-9]*$ ]] \
      || ! [[ "${count}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --min, --max, and --count must be positive integers." >&2
    exit 1
  fi

  if (( min > max )); then
    echo "Error: --min (${min}) must not exceed --max (${max})." >&2
    exit 1
  fi
}

generate_strands() {
  local -r file="$1"
  local -r min="$2"
  local -r max="$3"
  local -r count="$4"
  local -r seed="${SRANDOM:-$(date +%s%N)}"

  # Use AWK for unbounded random generation and fast text extraction.
  awk -v min="${min}" -v max="${max}" -v num="${count}" -v seed="${seed}" '
    BEGIN {
      srand(seed)
    }
    {
      # Strip newlines and carriage returns to keep text contiguous
      gsub(/\r/, "")
      text = text $0
    }
    END {
      len = length(text)
      if (len < min) {
        print "Error: Input text is shorter than minimum strand length." \
          > "/dev/stderr"
        exit 1
      }

      for (i = 1; i <= num; i++) {
        strand_len = min + int(rand() * (max - min + 1))

        max_offset = len - strand_len
        if (max_offset < 0) {
          strand_len = len
          max_offset = 0
        }

        # AWK strings are 1-indexed
        offset = int(rand() * (max_offset + 1)) + 1
        print substr(text, offset, strand_len)
      }
    }
  ' "${file}"
}

main() {
  local input_file=""
  local min_len=""
  local max_len=""
  local num_strands=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -i|--input)
        [[ $# -ge 2 ]] || {
          echo "Error: $1 requires an argument." >&2
          exit 1
        }
        input_file="$2"; shift 2 ;;
      -m|--min)
        [[ $# -ge 2 ]] || {
          echo "Error: $1 requires an argument." >&2
          exit 1
        }
        min_len="$2"; shift 2 ;;
      -M|--max)
        [[ $# -ge 2 ]] || {
          echo "Error: $1 requires an argument." >&2
          exit 1
        }
        max_len="$2"; shift 2 ;;
      -n|--count)
        [[ $# -ge 2 ]] || {
          echo "Error: $1 requires an argument." >&2
          exit 1
        }
        num_strands="$2"; shift 2 ;;
      -h|--help)
        usage; exit 0 ;;
      *)
        echo "Error: Unknown option $1" >&2
        usage >&2
        exit 1
        ;;
    esac
  done

  validate_args "${input_file}" "${min_len}" "${max_len}" "${num_strands}"
  generate_strands "${input_file}" "${min_len}" "${max_len}" "${num_strands}"
}

main "$@"

