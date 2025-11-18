#!/usr/bin/env bash
#
# Interactive input script with validation, confirmation, and full terminal line editing.
#

# Enable vi-mode for all inputs
set -o vi

# ------------------------------------------------------------
# Function: prompt_input
# Prompt user for input, validate, and confirm before accepting.
# Full line-editing (cursor movement, insertion, history) enabled.
# ------------------------------------------------------------
prompt_input() {
  local prompt="$1"
  local varname="$2"
  local is_array="$3" # "true" if input should be treated as an array
  local value confirm

  while true; do
    if [[ "$is_array" == "true" ]]; then
      # Use readline editing
      read -e -p "$prompt" -a value
      echo "$varname: ${value[*]}"
    else
      # Use readline editing for string input
      read -e -p "$prompt" value
      echo "$varname: $value"
    fi

    # Save to history so arrow keys recall previous inputs
    history -s "$value"

    # Confirm
    read -e -p "Proceed? (Y/n): " confirm
    if [[ "$confirm" =~ ^[Nn]$ ]]; then
      echo "Re-enter $varname..."
      continue
    fi

    # Validation
    if [[ "$is_array" == "true" ]]; then
      if declare -F validate_array_input >/dev/null; then
        validate_array_input value || {
          echo "Invalid input. Please try again."
          echo
          continue
        }
      fi
      eval "$varname=(\"\${value[@]}\")"
    else
      if declare -F validate_string_input >/dev/null; then
        validate_string_input "$value" || {
          echo "Invalid input. Please try again."
          echo
          continue
        }
      fi
      eval "$varname=\"$value\""
    fi

    break
  done
}

# ------------------------------------------------------------
# Function: validate_array_input
# Validate an array of target genes.
# ------------------------------------------------------------
validate_array_input() {
  local -n arr=$1
  if [[ ${#arr[@]} -eq 0 ]]; then
    return 1 # invalid: empty
  fi

  for gene in "${arr[@]}"; do
    if [[ ! $gene =~ ^[A-Za-z0-9*]+$ ]]; then
      echo "Invalid gene format: $gene"
      return 1
    fi
  done
  return 0
}

# ------------------------------------------------------------
# Function: validate_string_input
# Validate a string input (e.g., organism name).
# ------------------------------------------------------------
validate_string_input() {
  local input="$1"
  if [[ -z "$input" ]]; then
    echo "Input cannot be empty."
    return 1
  fi

  if [[ ! $input =~ ^[A-Za-z0-9\ ]+$ ]]; then
    echo "Invalid input: $input (only letters, numbers, spaces allowed)"
    return 1
  fi

  return 0
}
