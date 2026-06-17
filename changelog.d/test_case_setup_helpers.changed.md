Factored the two repeated, gotcha-prone ``UserData`` idioms in the regression
cases into ``tests/case_setup.py``: ``build_bdry`` (the per-instance,
object-dtype boundary-type triple — never a shared class attribute) and
``make_diag_state`` (centralises the ``Nx = inx - 1`` / ``Ny = iny - 1`` /
``steps = [stepmax - 1]`` offsets while forwarding case-specific keywords).
Behaviour-preserving: ``vars(UserData())`` is byte-identical for all 11 cases
and the regression gate stays bit-identical.
