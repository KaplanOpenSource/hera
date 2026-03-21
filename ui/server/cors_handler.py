import argparse
import sys

LOCAL_ORIGINS = [f'http://{h}:{p}' for h in ['localhost', '127.0.0.1', '0.0.0.0'] for p in [5173, 8000]]


class CorsHandler:
    def __init__(self) -> None:
        self.custom_origins: list[str] | None = None

    def add_argument(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            '--cors',
            default=None,
            metavar='ORIGINS',
            help=(
                'Enable CORS for external origins. '
                'Use "all" to allow all origins. '
                'Pass a comma-separated list of IPs to allow specific ones '
                '(e.g. --cors 192.168.1.10,10.0.0.5). '
                'Each IP is prefixed with http:// and port 8000 automatically.'
            ),
        )

    def get_origins(self, args: argparse.Namespace) -> list[str]:
        if args.cors is None:
            return LOCAL_ORIGINS

        if args.cors == 'all':
            self.custom_origins = ['*']
            warning = "WARNING: CORS is enabled for ALL origins. This is insecure in production."
        else:
            self.custom_origins = [f'http://{ip}:8000' for ip in args.cors.split(',')]
            warning = f"WARNING: CORS is enabled for custom origins: {', '.join(self.custom_origins)}"

        print(warning)
        if not args.yes:
            try:
                answer = input("CORS weakens browser security. Continue? (use -y to skip) [y/N] ").strip().lower()
            except EOFError:
                answer = ''
            if answer != 'y':
                print("Aborted.")
                sys.exit(0)

        if args.cors == 'all':
            return ['*']
        return LOCAL_ORIGINS + self.custom_origins
