import ipywidgets as widgets

test_case = widgets.Dropdown(
    options=[('travelling vortex', 1), ('acoustic wave high', 2), ('long internal wave', 3), ('rising bubble',4)],
    value=2,
    description='test case:',
)

labels = widgets.Dropdown(
    options=[('initial condition', 'ic'),\
            ('before first advection routine', 'before_advect'),\
            ('after first advection routine', 'after_advect'),\
            ('after first euler backward expl part', 'after_ebnaexp'),\
            ('after first euler backward impl part', 'after_ebnaimp'),\
            ('after the half step', 'after_half_step'),\
            ('after the explicit euler timestep', 'after_efna'),\
            ('after advection routine for full timestep','after_full_advect'),\
            ('after full euler backward expl part', 'after_full_ebnaexp'),\
            ('after a full timestep', 'after_full_step')
    ],
    value='after_full_step',
    description='time label:',
)

time_step = widgets.Text(
    value='000',
    description='time-step:',
    disabled=False,
)