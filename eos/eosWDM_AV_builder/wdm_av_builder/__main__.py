"""Entry point so ``python -m wdm_av_builder ...`` works."""

from __future__ import annotations

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
