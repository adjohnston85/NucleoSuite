"""Consistent command-line help formatting for NucleoSuite."""

from __future__ import annotations

import argparse
import re


_DEFAULT_WORD = re.compile(r"\bdefault\b", flags=re.IGNORECASE)


class NucleoSuiteHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """Show each option default once while retaining formatted descriptions.

    Older command parsers often stated a default explicitly in their help text
    and also used :class:`argparse.ArgumentDefaultsHelpFormatter`.  Argparse
    consequently appended a second ``(default: ...)`` annotation.  Explicit
    wording is retained; automatic wording is added only when the option help
    does not already mention its default.
    """

    def _get_help_string(self, action: argparse.Action) -> str:
        help_text = action.help or ""
        if _DEFAULT_WORD.search(help_text):
            return help_text
        return super()._get_help_string(action)

