"""Path setup for the validation test suite: make both the ``validation``
directory (for ``hermite_eval`` / ``heldout``) and the builder package
importable regardless of where pytest is invoked."""

from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_VALIDATION = os.path.dirname(_HERE)
_BUILDER = os.path.dirname(_VALIDATION)
for p in (_VALIDATION, _BUILDER):
    if p not in sys.path:
        sys.path.insert(0, p)
