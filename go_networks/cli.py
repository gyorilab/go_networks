from typing import Optional

from generate_v2 import main as gen_networks

import click


@click.group("go-networks")
def main():
    """Update the go_networks"""


@main.command()
@click.option(
    "--regenerate-props",
    is_flag=True,
    help="Regenerate the properties file. Use this option if the underlying "
    "data has changed. If there is no local props file cached, a new one "
    "will be generated regardless of this flag.",
)
@click.option(
    "--go-term",
    type=str,
    default="GO:2001239",
    help="The GO term to generate the network for. Default: GO:2001239",
)
@click.option(
    "--style-network",
    type=str,
    default="4c2006cd-9fef-11ec-b3be-0ac135e8bacf",
    help="Network ID of the style network",
)
@click.option(
    "--local-sif",
    type=str,
    help="Path to a local SIF file. Only needed if we are regenerating the "
    "properties file.",
)
def test(
    regenerate_props: bool,
    go_term: str,
    style_network: str,
    local_sif: Optional[str] = None,
):
    """Build a test network."""
    # Use the test network set:
    # d7acdc1d-a08f-11ec-b3be-0ac135e8bacf
    gen_networks(
        local_sif=local_sif,
        network_set="d7acdc1d-a08f-11ec-b3be-0ac135e8bacf",
        style_network=style_network,
        regenerate=regenerate_props,
        test_go_term=go_term,
    )


@main.command()
@click.option(
    "--regenerate-props",
    is_flag=True,
    help="Regenerate the properties file. Use this option if the underlying "
    "data has changed.",
)
@click.option(
    "--style-network",
    type=str,
    default="4c2006cd-9fef-11ec-b3be-0ac135e8bacf",
    help="Network ID of the style network",
)
@click.option(
    "--network-set",
    help="Network set ID to add the new networks to.",
    default="bdba6a7a-488a-11ec-b3be-0ac135e8bacf",
)
@click.option(
    "--local-sif",
    type=str,
    help="Path to a local SIF file to use instead of downloading from S3.",
)
def run(
    regenerate_props: bool,
    style_network: str,
    network_set: str,
    local_sif: Optional[str] = None,
):
    """Run the go network generation."""
    gen_networks(
        local_sif=local_sif,
        network_set=network_set,
        style_network=style_network,
        regenerate=regenerate_props,
    )


if __name__ == "__main__":
    main()
