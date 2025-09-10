"""
OXASL module for multi-TE ASL data

This module is designed to operate within the OXASL pipeline.
If installed, then it will be called by ``oxasl.oxford_asl.oxasl``
whenever multi-TE data is supplied.

The relevant processing function can also be called independently
on a ``Workspace`` object, however this will not include the
standard oxasl preprocessing or registration.
"""
from .api import run, model_multite, MultiTEOptions
from ._version import __version__

__all__ = ["__version__", "run", "model_multite", "MultiTEOptions"]
