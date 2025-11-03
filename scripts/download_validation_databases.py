#!/usr/bin/env python3
"""
Download biological validation databases for local use.

Usage:
    python scripts/download_validation_databases.py --all
    python scripts/download_validation_databases.py --databases hpo gwas_catalog string
    python scripts/download_validation_databases.py --tier 1
    python scripts/download_validation_databases.py --update-only
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import List, Dict

from rich.console import Console
from rich.table import Table

from sae_genomics.validation.downloaders import DOWNLOADERS

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

console = Console()

# Database tiers (by priority/size)
TIERS = {
    1: ['gencc', 'hpo', 'string'],  # Small, essential (~1 GB total)
    2: ['gwas_catalog', 'orphanet', 'clinvar'],  # Medium size (~2 GB total)
    3: ['monarch', 'open_targets'],  # Large (~10-15 GB total)
    # API-only databases (cache setup, no large downloads)
    'api': ['disgenet', 'omim'],
}


def download_database(
    name: str,
    output_dir: Path,
    force: bool = False,
    verify: bool = True,
    **kwargs
) -> Dict:
    """Download a single database.

    Args:
        name: Database name
        output_dir: Output directory
        force: Force re-download even if exists
        verify: Verify integrity after download
        **kwargs: Additional arguments for specific downloaders

    Returns:
        Metadata dict with download info

    Raises:
        Exception: If download or verification fails
    """
    logger.info(f"Downloading {name}...")

    downloader_cls = DOWNLOADERS[name]
    downloader = downloader_cls(output_dir=output_dir / name, **kwargs)

    # Check if already downloaded
    if not force and downloader.is_downloaded():
        console.print(f"[yellow]✓ {name} already downloaded. Use --force to re-download.[/yellow]")
        return downloader.get_metadata()

    # Download
    try:
        metadata = downloader.download()

        # Verify integrity
        if verify:
            console.print(f"[blue]Verifying {name}...[/blue]")
            if not downloader.verify():
                raise ValueError(f"Verification failed for {name}")

        console.print(f"[green]✓ {name} downloaded successfully[/green]")
        return metadata

    except Exception as e:
        console.print(f"[red]✗ Failed to download {name}: {e}[/red]")
        raise


def print_download_summary(results: Dict[str, Dict]):
    """Print a summary table of download results.

    Args:
        results: Dict mapping database names to result dicts
    """
    table = Table(title="Download Summary")
    table.add_column("Database", style="cyan")
    table.add_column("Status", style="green")
    table.add_column("Version", style="yellow")
    table.add_column("Size (MB)", justify="right", style="magenta")
    table.add_column("Files", justify="right")

    for name, result in results.items():
        if result['status'] == 'success':
            metadata = result['metadata']
            total_size = sum(f['size_mb'] for f in metadata.get('files', []))
            n_files = len(metadata.get('files', []))

            table.add_row(
                name,
                "✓ Success",
                metadata.get('version', 'unknown'),
                f"{total_size:.1f}",
                str(n_files)
            )
        else:
            table.add_row(
                name,
                "✗ Failed",
                "-",
                "-",
                "-"
            )

    console.print(table)


def main():
    parser = argparse.ArgumentParser(
        description="Download validation databases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download all databases
  python scripts/download_validation_databases.py --all

  # Download specific databases
  python scripts/download_validation_databases.py --databases hpo gencc string

  # Download by tier
  python scripts/download_validation_databases.py --tier 1

  # Force re-download
  python scripts/download_validation_databases.py --databases hpo --force
        """
    )

    parser.add_argument(
        '--databases',
        nargs='+',
        choices=list(DOWNLOADERS.keys()),
        help='Specific databases to download'
    )
    parser.add_argument(
        '--tier',
        type=int,
        choices=[1, 2, 3],
        help='Download all databases in a tier (1=small, 2=medium, 3=large)'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Download all databases'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/validation_databases'),
        help='Output directory for databases (default: data/validation_databases)'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force re-download even if files exist'
    )
    parser.add_argument(
        '--no-verify',
        action='store_true',
        help='Skip integrity verification'
    )
    parser.add_argument(
        '--update-only',
        action='store_true',
        help='Only update databases that have new versions'
    )

    # API keys for restricted databases
    parser.add_argument('--disgenet-api-key', help='DisGeNET API key')
    parser.add_argument('--omim-api-key', help='OMIM API key')

    args = parser.parse_args()

    # Determine which databases to download
    if args.all:
        databases = list(DOWNLOADERS.keys())
        console.print("[bold]Downloading all databases...[/bold]")
    elif args.tier:
        databases = TIERS[args.tier]
        console.print(f"[bold]Downloading Tier {args.tier} databases: {', '.join(databases)}[/bold]")
    elif args.databases:
        databases = args.databases
        console.print(f"[bold]Downloading: {', '.join(databases)}[/bold]")
    else:
        parser.error("Must specify --all, --tier, or --databases")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    console.print(f"Output directory: [cyan]{args.output_dir}[/cyan]\n")

    # Load existing metadata
    log_file = args.output_dir / 'download_log.json'
    if log_file.exists():
        with open(log_file) as f:
            download_log = json.load(f)
    else:
        download_log = {'last_updated': None, 'databases': {}}

    # Download each database
    results = {}
    for name in databases:
        try:
            # Skip if update-only and no new version
            # TODO: Implement version checking logic

            metadata = download_database(
                name,
                args.output_dir,
                force=args.force,
                verify=not args.no_verify,
                disgenet_api_key=args.disgenet_api_key,
                omim_api_key=args.omim_api_key,
            )
            results[name] = {'status': 'success', 'metadata': metadata}

        except Exception as e:
            logger.error(f"Failed to download {name}: {e}")
            results[name] = {'status': 'failed', 'error': str(e)}
            # Continue with other databases

    # Update download log
    download_log['last_updated'] = datetime.now().isoformat()
    for name, result in results.items():
        if result['status'] == 'success':
            download_log['databases'][name] = result['metadata']

    with open(log_file, 'w') as f:
        json.dump(download_log, f, indent=2)

    # Print summary
    console.print()
    print_download_summary(results)

    # Overall summary
    success = sum(1 for r in results.values() if r['status'] == 'success')
    failed = sum(1 for r in results.values() if r['status'] == 'failed')

    console.print(f"\n[bold]Download complete: {success} successful, {failed} failed[/bold]")
    console.print(f"Download log saved to [cyan]{log_file}[/cyan]")

    if failed > 0:
        console.print("\n[red]Failed databases:[/red]")
        for name, result in results.items():
            if result['status'] == 'failed':
                console.print(f"  - [red]{name}[/red]: {result['error']}")
        return 1

    console.print("\n[green]✓ All downloads completed successfully![/green]")
    return 0


if __name__ == '__main__':
    sys.exit(main())
