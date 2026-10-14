"""The server's command line arguments, parsed once and read by the routers."""
from __future__ import annotations

import argparse
from typing import List, Optional

import argcomplete

from cors_handler import CorsHandler


class ArgsSettings:
    """The parsed command line, plus the CORS origins it asks for."""

    def __init__(self, args: argparse.Namespace, cors_handler: CorsHandler) -> None:
        self.args = args
        self.cors_handler = cors_handler

    @property
    def host(self) -> str:
        return self.args.host

    @property
    def port(self) -> int:
        return self.args.port

    @property
    def jupyter_port(self) -> int:
        return self.args.jupyter_port

    @property
    def debug(self) -> bool:
        return self.args.debug

    # Jupyter is only exposed beyond localhost when the user opted into CORS.
    @property
    def remote_allowed(self) -> bool:
        return self.args.cors is not None

    @property
    def custom_origins(self) -> Optional[List[str]]:
        return self.cors_handler.custom_origins

    def origins(self) -> List[str]:
        return self.cors_handler.get_origins(self.args)


def parse_args_settings(argv: Optional[List[str]] = None) -> ArgsSettings:
    cors_handler = CorsHandler()
    parser = argparse.ArgumentParser(description="Hera UI API server")
    cors_handler.add_argument(parser)
    parser.add_argument('--debug', action='store_true', help='Enable debugpy remote debugging on port 5678')
    parser.add_argument('-y', '--yes', action='store_true', help='Skip confirmation prompts')
    parser.add_argument('--host', default='0.0.0.0', help='Address to bind to (default: 0.0.0.0, all interfaces; use 127.0.0.1 for local-only, e.g. behind a reverse proxy).')
    parser.add_argument('--port', type=int, default=8000, help='Port for the API server')
    parser.add_argument('--jupyter-port', type=int, default=8888, help='Port for Jupyter server (0 to disable)')
    argcomplete.autocomplete(parser)
    return ArgsSettings(parser.parse_args(argv), cors_handler)


_settings: Optional[ArgsSettings] = None


def set_args_settings(settings: ArgsSettings) -> None:
    global _settings
    _settings = settings


def get_args_settings() -> ArgsSettings:
    """The arguments in force; plain defaults when the server was not started from the CLI (tests)."""
    global _settings
    if _settings is None:
        _settings = parse_args_settings([])
    return _settings
