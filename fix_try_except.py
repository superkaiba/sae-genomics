#!/usr/bin/env python3
"""Remove try/except blocks from local_client.py"""

import re
from pathlib import Path

file_path = Path("/workspace/research/projects/sae-genomics/src/sae_genomics/validation/local_client.py")

with open(file_path, 'r') as f:
    content = f.read()

# Pattern to match try/except blocks in the database loading methods
# These follow the pattern:
#   try:
#       code
#   except Exception as e:
#       self._print(f"[red]✗ Failed to load X: {e}[/red]")

# Find and remove try blocks (just the "try:" line and dedent the content)
lines = content.split('\n')
new_lines = []
in_try_block = False
try_indent = 0
skip_next_except = False

i = 0
while i < len(lines):
    line = lines[i]

    # Detect try block in a method
    if line.strip() == 'try:' and i > 0:
        # Check if previous line is a docstring or method def
        prev_lines = [lines[j] for j in range(max(0, i-5), i) if lines[j].strip()]
        if any('def _load_' in l for l in prev_lines[-3:]):
            in_try_block = True
            try_indent = len(line) - len(line.lstrip())
            skip_next_except = True
            i += 1
            continue

    # Detect except block to remove
    if skip_next_except and line.strip().startswith('except Exception'):
        # Skip this line and the next line (the print statement)
        i += 1
        if i < len(lines) and '_print' in lines[i] and 'Failed to load' in lines[i]:
            i += 1
        skip_next_except = False
        in_try_block = False
        continue

    # Dedent lines that were inside try block
    if in_try_block and line.startswith(' ' * (try_indent + 4)):
        # Remove 4 spaces of indentation
        new_lines.append(line[4:])
    else:
        new_lines.append(line)

    i += 1

# Write back
with open(file_path, 'w') as f:
    f.write('\n'.join(new_lines))

print("Removed try/except blocks from local_client.py")
