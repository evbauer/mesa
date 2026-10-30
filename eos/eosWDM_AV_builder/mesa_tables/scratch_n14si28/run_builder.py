import sys
# Make the trusted wdm_av_builder package importable under python3 -I without
# adding the cwd (which may be an untrusted data dir). Only this known-good
# project path is inserted.
sys.path.insert(0, '/Users/bauer39/work/mesas/evbauer/mesa/eos/eosWDM_AV_builder')
from wdm_av_builder.cli import main
sys.exit(main())
