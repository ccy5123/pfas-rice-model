"""UI layer for the PFAS-rice Streamlit dashboard.

Split out of the original monolithic app.py (HANDOFF P3-1) into:
  common.py  -- constants, cached model helpers, shared render building blocks
  sidebar.py -- the sidebar (Simple vs Expert); returns a config namespace
  simple.py  -- the general-audience (Korean) view
  expert.py  -- the full research (English) view
app.py is now just the assembler.
"""
import os
import sys

_SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
