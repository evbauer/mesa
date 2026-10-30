import os
import sys
# Make the trusted wdm_av_builder package importable under python3 -I without
# adding the cwd (which may be an untrusted data dir). Only the builder's own
# directory (two levels up from this script) is inserted.
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from wdm_av_builder.cli import main
sys.exit(main())
