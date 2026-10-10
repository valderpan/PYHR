#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author : Haoran Pan
# Date: 2026/06/30



"""
Terminal color constants and helper functions for PYHR command-line output.

This module mainly provides ANSI colors for argparse/help messages and
Rich-style color names for rich logging or Console output.
"""

import os
import sys


# =============================================================================
# Color control
# =============================================================================

def _support_color() -> bool:
    """
    Return True if terminal color should be enabled.

    Color will be disabled when:
    1. stdout is not a terminal, e.g. redirecting help to a file
    2. environment variable NO_COLOR is set
    """
    return sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


ENABLE_COLOR = _support_color()


# =============================================================================
# ANSI color codes
# Suitable for argparse help text
# =============================================================================

if ENABLE_COLOR:
    RESET = "\033[0m"
    BOLD = "\033[1m"

    BLACK = "\033[30m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"

    # 256-color / true-color style
    BURGUNDY = "\033[38;2;128;0;32m"
    ORANGE = "\033[38;2;255;140;0m"
    LIGHT_BLUE = "\033[38;2;80;160;255m"
    LIGHT_GREEN = "\033[38;2;80;200;120m"
else:
    RESET = ""
    BOLD = ""

    BLACK = ""
    RED = ""
    GREEN = ""
    YELLOW = ""
    BLUE = ""
    MAGENTA = ""
    CYAN = ""
    WHITE = ""

    BURGUNDY = ""
    ORANGE = ""
    LIGHT_BLUE = ""
    LIGHT_GREEN = ""


def color_text(text: str, color: str, bold: bool = False) -> str:
    """
    Wrap text with ANSI color.

    Parameters
    ----------
    text : str
        Text to colorize.
    color : str
        ANSI color variable, e.g. BURGUNDY, GREEN, CYAN.
    bold : bool
        Whether to make text bold.

    Returns
    -------
    str
        Colored text if terminal supports color, otherwise plain text.
    """
    if not ENABLE_COLOR:
        return text

    prefix = BOLD + color if bold else color
    return f"{prefix}{text}{RESET}"


# =============================================================================
# Rich style names
# Suitable for rich console/logging markup
# =============================================================================

RICH_BURGUNDY = "rgb(128,0,32)"
RICH_ORANGE = "rgb(255,140,0)"
RICH_LIGHT_BLUE = "rgb(80,160,255)"
RICH_LIGHT_GREEN = "rgb(80,200,120)"

RICH_SUCCESS = "bold green"
RICH_WARNING = "bold yellow"
RICH_ERROR = "bold red"
RICH_INFO = "bold cyan"


def color_command_line(description, command, color, bold=True):
    return description.replace(
        command,
        color_text(command, color, bold=bold)
    )

