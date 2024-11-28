# -*- coding: utf-8 -*-
"""Custom error types."""


class LowCoordinationNumber(KeyError):
    """Error for low coordination number."""


class HighCoordinationNumber(KeyError):
    """Error for high coordination number."""


class NoOpenDefined(KeyError):
    """Error in case the open check is not defined for this coordination number."""
