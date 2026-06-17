Fixed a latent ``TypeError`` in ``EnsembleState.set_members``: the assertion
``len(self.set_members == members)`` compared a bound method to a list and then
called ``len`` on the resulting bool. It now checks ``len(members) == len(self.members)``
(the ensemble size is preserved when forecast members are replaced by analysis members).
