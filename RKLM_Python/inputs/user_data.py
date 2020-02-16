
class UserDataInit(object):
    """
    Loads user defined initial conditions. Specifically, all attributes of the class object defined in the initial condition is loaded.

    Attributes
    ----------
    **kwargs: class object

    Todo
    -----
    List the essential attributes necessary for the simulation to run.
    """
    def __init__(self,**kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)