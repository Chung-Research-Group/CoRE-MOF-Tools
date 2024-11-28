# -*- coding: utf-8 -*-
"""Encode chemical elements numerically and decode numerical representations of elements."""
from loguru import logger

from .decode import decode, decode_many
from .encode import encode, encode_many

__all__ = ("encode", "encode_many", "decode", "decode_many")
logger.disable("element_coder")
