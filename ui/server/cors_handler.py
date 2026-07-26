import argparse
import socket
import sys

DEFAULT_SERVER_PORT = 8000


class CorsHandler:
    def __init__(self) -> None:
        self.custom_origins: list[str] | None = None

    def _local_origins(self, server_port: int = DEFAULT_SERVER_PORT) -> list[str]:
        return [f'http://{h}:{p}' for h in ['localhost', '127.0.0.1', '0.0.0.0'] for p in [5173, server_port]]

    def _remote_origins(self, server_port: int = DEFAULT_SERVER_PORT) -> list[str]:
        """Origins for this machine's own hostname and LAN IP(s), so it can be reached by name from other machines."""
        hosts: set[str] = set()
        hostname = socket.gethostname()
        if hostname:
            hosts.add(hostname)
            try:
                hosts.add(socket.gethostbyname(hostname))
            except OSError:
                # hostname may not resolve via DNS/hosts file; skip the IP form.
                pass
        try:
            # Pick the IP of the interface used to reach the network (no packets are sent for UDP connect).
            with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as probe:
                probe.connect(('8.8.8.8', 80))
                hosts.add(probe.getsockname()[0])
        except OSError:
            # No route to the network (offline); the hostname entry above still applies.
            pass
        return [f'http://{h}:{p}' for h in sorted(hosts) for p in [5173, server_port]]

    def add_argument(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            '--cors',
            default=None,
            metavar='ORIGINS',
            help=(
                'Enable CORS for external origins. '
                'Use "all" to allow all origins. '
                'Use "remote" to also allow this machine\'s own hostname and LAN IP '
                '(so it can be opened by name from other machines). '
                'Pass a comma-separated list of IPs to allow specific ones '
                '(e.g. --cors 192.168.1.10,10.0.0.5). '
                'Each IP is prefixed with http:// and the server port automatically. '
                'Only "all" asks for confirmation; "remote" and explicit IPs just print a warning.'
            ),
        )

    def get_origins(self, args: argparse.Namespace) -> list[str]:
        server_port = getattr(args, 'port', DEFAULT_SERVER_PORT)

        if args.cors is None:
            return self._local_origins(server_port)

        if args.cors == 'all':
            self.custom_origins = ['*']
            warning = "WARNING: CORS is enabled for ALL origins. This is insecure in production."
        elif args.cors == 'remote':
            self.custom_origins = self._remote_origins(server_port)
            warning = f"WARNING: CORS is enabled for this machine's remote origins: {', '.join(self.custom_origins)}"
        else:
            self.custom_origins = [f'http://{ip}:{server_port}' for ip in args.cors.split(',')]
            warning = f"WARNING: CORS is enabled for custom origins: {', '.join(self.custom_origins)}"

        print(warning)
        # Only "all" (['*']) is a real security downgrade, so only it needs confirmation.
        # "remote" and explicit IP lists are limited, deliberate origins — just warn and proceed.
        if args.cors == 'all' and not args.yes:
            try:
                answer = input("CORS weakens browser security. Continue? (use -y to skip) [y/N] ").strip().lower()
            except EOFError:
                answer = ''
            if answer != 'y':
                print("Aborted.")
                sys.exit(0)

        if args.cors == 'all':
            return ['*']
        return self._local_origins(server_port) + self.custom_origins
