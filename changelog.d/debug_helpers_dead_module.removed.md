Removed the dead `utils/debug_helpers.py` module: its only symbol `pl_sol` had
zero call sites anywhere in the tree, and the module pulled in a top-level
`matplotlib.pyplot` import for nothing.
