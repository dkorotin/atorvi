#!/bin/bash

# to produce pdf from .md file I'm using pandoc with eisvogel template https://github.com/Wandmalfarbe/pandoc-latex-template

pandoc atorvi_manual.md -f gfm -o atorvi_manual.pdf --template eisvogel --listings