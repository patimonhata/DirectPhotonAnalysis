#!/bin/bash
set -euo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
input_directory=${1:-$module_dir/input/jet5}
output_file=${2:-$input_directory/segments.list}

declare -a list_names=(
  dst_calo_cluster.list
  dst_mbd_epd.list
  dst_truth_jet.list
  g4hits.list
)
declare -a prefixes=(
  DST_CALO_CLUSTER_
  DST_MBD_EPD_
  DST_TRUTH_JET_
  G4Hits_
)

temporary_directory=$(mktemp -d)
trap 'rm -rf "$temporary_directory"' EXIT

for index in "${!list_names[@]}"; do
  input_file="$input_directory/${list_names[$index]}"
  prefix=${prefixes[$index]}
  normalized_file="$temporary_directory/${list_names[$index]}"
  if [[ ! -s "$input_file" ]]; then
    echo "Missing or empty input list: $input_file" >&2
    exit 2
  fi
  awk -v prefix="$prefix" '
    index($0, prefix) != 1 || $0 !~ /\.root$/ {
      printf "Malformed row %d in %s: %s\n", NR, FILENAME, $0 > "/dev/stderr"
      failed = 1
      next
    }
    { print substr($0, length(prefix) + 1) }
    END { if (failed) exit 1 }
  ' "$input_file" > "$normalized_file"
done

reference_file="$temporary_directory/${list_names[0]}"
for list_name in "${list_names[@]:1}"; do
  candidate_file="$temporary_directory/$list_name"
  if ! cmp -s "$reference_file" "$candidate_file"; then
    echo "Stream lists do not contain identical ordered suffixes:" >&2
    echo "  $input_directory/${list_names[0]}" >&2
    echo "  $input_directory/$list_name" >&2
    diff -u "$reference_file" "$candidate_file" | head -n 40 >&2 || true
    exit 3
  fi
done

if [[ $(sort "$reference_file" | uniq -d | head -n 1) ]]; then
  echo "Duplicate suffix found in input lists" >&2
  exit 4
fi

mkdir -p "$(dirname "$output_file")"
temporary_output=$(mktemp "$(dirname "$output_file")/.segments.XXXXXX")
cp "$reference_file" "$temporary_output"
chmod 0644 "$temporary_output"
mv "$temporary_output" "$output_file"
echo "Wrote $(wc -l < "$output_file") synchronized input suffixes to $output_file"
