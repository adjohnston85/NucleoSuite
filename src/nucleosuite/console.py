"""Shared console messages for NucleoSuite commands."""

from __future__ import annotations

import random
from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class CommandMessagePair:
    """A command-start message and its matching success message."""

    startup: str
    completion: str


ORC_COMPLETION_MESSAGE = "Work complete."
PEASANT_COMPLETION_MESSAGE = "Job's done."

ORC_STARTUP_MESSAGES: tuple[str, ...] = (
    "Zug zug.",
    "Dabu.",
    "Swobu.",
    "Lok'tar.",
)

PEASANT_STARTUP_MESSAGES: tuple[str, ...] = (
    "Okay.",
    "Right-o.",
    "Alright.",
    "Yes, milord.",
)

MESSAGE_PAIRS: tuple[CommandMessagePair, ...] = tuple(
    CommandMessagePair(message, ORC_COMPLETION_MESSAGE)
    for message in ORC_STARTUP_MESSAGES
) + tuple(
    CommandMessagePair(message, PEASANT_COMPLETION_MESSAGE)
    for message in PEASANT_STARTUP_MESSAGES
)

STARTUP_MESSAGES: tuple[str, ...] = tuple(pair.startup for pair in MESSAGE_PAIRS)


def message_pair(
    pairs: Sequence[CommandMessagePair] = MESSAGE_PAIRS,
) -> CommandMessagePair:
    """Return one randomly selected, internally matched message pair."""
    if not pairs:
        raise ValueError("At least one command message pair is required.")
    return random.choice(tuple(pairs))


def startup_message(messages: Sequence[str] = STARTUP_MESSAGES) -> str:
    """Return one randomly selected command-start message.

    Use :func:`message_pair` when the matching completion message is also required.
    """
    if not messages:
        raise ValueError("At least one startup message is required.")
    return random.choice(tuple(messages))


def print_startup(message: str | None = None) -> None:
    """Print a command-start message."""
    print(message if message is not None else startup_message(), flush=True)


def print_completion(message: str = ORC_COMPLETION_MESSAGE) -> None:
    """Print the successful command-completion message."""
    print(message, flush=True)
