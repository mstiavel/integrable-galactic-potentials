#!/bin/sh
# Build docs/research_note.pdf from docs/research_note.md (needs pandoc and LuaLaTeX).
cd "$(dirname "$0")/.." || exit 1
pandoc -f markdown+lists_without_preceding_blankline docs/research_note.md -o docs/research_note.pdf --pdf-engine=lualatex \
  -V mainfont="STIX Two Text" -V 'mainfontfallback=STIX Two Math:mode=harf' \
  -V 'mainfontfallback=Apple Symbols:mode=harf' \
  -V monofont="Menlo" -V 'monofontfallback=STIX Two Math:mode=harf' \
  -V geometry:margin=2.5cm -V fontsize=11pt -V colorlinks --toc \
  -M title="Integrable galactic potentials beyond Stäckel: research note" \
  -M author="Massimo Stiavelli (Space Telescope Science Institute)" \
  -M date="September 2026"
