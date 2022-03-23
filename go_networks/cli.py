from typing import Optional

from generate_v2 import main as gen_networks

import click


TEST_SET = "d7acdc1d-a08f-11ec-b3be-0ac135e8bacf"
STYLE_NX = "4c2006cd-9fef-11ec-b3be-0ac135e8bacf"
PROD_SET = "303190ca-aac0-11ec-b3be-0ac135e8bacf"


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
    # See style network at
    # https://www.ndexbio.org/viewer/networks/4c2006cd-9fef-11ec-b3be-0ac135e8bacf
    default=STYLE_NX,
    help=f"Network ID of the style network. Default: {STYLE_NX}",
)
@click.option(
    "--local-sif",
    type=str,
    help="Path to a local SIF file. Only needed if we are regenerating the "
    "properties file.",
)
@click.option(
    "--ndex-server-style",
    type=str,
    default="http://ndexbio.org",
    help="The URL of the ndex server to use for the style network. "
    "Default: http://ndexbio.org",
)
def test(
    regenerate_props: bool,
    go_term: str,
    style_network: str,
    local_sif: Optional[str] = None,
    ndex_server_style: str = "http://ndexbio.org",
):
    """Build a test network."""
    # Use the test network set:
    # d7acdc1d-a08f-11ec-b3be-0ac135e8bacf
    gen_networks(
        local_sif=local_sif,
        network_set=TEST_SET,
        style_network=style_network,
        regenerate=regenerate_props,
        test_go_term=go_term,
        ndex_server_style=ndex_server_style,
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
    default=STYLE_NX,
    help=f"Network ID of the style network. Default: {STYLE_NX}",
)
@click.option(
    "--network-set",
    help=f"Network set ID to add the new networks to. Default: {PROD_SET}",
    default=PROD_SET,
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
