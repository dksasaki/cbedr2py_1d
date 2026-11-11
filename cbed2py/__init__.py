"""
cbed2py -- transient CBED (v2 coefficients) at a single point, through rpy2.

    import cbed2py
    p   = cbed2py.default_params()
    ss  = cbed2py.steady_state(params)
    out = cbed2py.run_transient(times, forcing=..., params=..., yini=None)

See README.md of the project for inputs, parameters and output.
"""
from .wrapper import (
    default_params,
    steady_state,
    run_transient,
    set_model,
    get_model,
    MODELS,
    JPOC_MODE_PARS,
    JPOC_MODES,
    FORCING,
    STATES,
    MODEL_DIR,
)

__all__ = ["default_params", "steady_state", "run_transient",
           "set_model", "get_model", "MODELS",
           "JPOC_MODE_PARS", "JPOC_MODES", "FORCING", "STATES", "MODEL_DIR"]
__version__ = "0.1.0"
