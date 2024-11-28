# -*- coding: utf-8 -*-
"""Base classes for checks."""
import abc
from typing import List

from backports.cached_property import cached_property


class AbstractCheck(abc.ABC):
    """Base class for checks."""

    @property
    @abc.abstractmethod
    def description(self) -> str:
        """Return a description of the check."""
        pass

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """Return the name of the check."""
        pass

    @cached_property
    def is_ok(self) -> bool:
        """Return whether the check passed."""
        return self._run_check()

    @abc.abstractmethod
    def _run_check(self):
        pass


class AbstractIndexCheck(abc.ABC):
    """Base class for checks that return indices."""

    @property
    @abc.abstractmethod
    def description(self):
        """Return a description of the check."""
        pass

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """Return the name of the check."""
        pass

    @cached_property
    def is_ok_and_indices(self):
        """Return whether the check passed and the indices that failed."""
        result, indices = self._run_check()
        return result, indices

    @cached_property
    def is_ok(self) -> bool:
        """Return whether the check passed."""
        result, _ = self.is_ok_and_indices
        return result

    @cached_property
    def flagged_indices(self) -> List[int]:
        """Return the indices that failed the check."""
        _, indices = self.is_ok_and_indices
        return indices

    @abc.abstractmethod
    def _run_check(self):
        pass


class AbstractMissingCheck(abc.ABC):
    """Base class for checks that return candidate positions."""

    @property
    @abc.abstractmethod
    def description(self) -> str:
        """Return a description of the check."""
        pass

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """Return the name of the check."""
        pass

    @cached_property
    def is_ok_indices_positions(self):
        """Return whether the check passed and the indices and positions that failed."""
        result, indices, positions = self._run_check()
        return result, indices, positions

    @property
    def is_ok(self) -> bool:
        """Return whether the check passed."""
        result, _, _ = self.is_ok_indices_positions
        return result

    @property
    def flagged_indices(self) -> List[int]:
        """Return the indices that failed the check."""
        _, indices, _ = self.is_ok_indices_positions
        return indices

    @property
    def candidate_positions(self):
        """Return the candidate positions for the missing sites."""
        _, _, positions = self.is_ok_indices_positions
        return positions

    @abc.abstractmethod
    def _run_check(self):
        pass
