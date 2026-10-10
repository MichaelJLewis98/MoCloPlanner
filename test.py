"""Compatibility entry point; use moclo_cli.py for the command-line app."""

from moclo_core import *  # noqa: F401,F403
from moclo_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
