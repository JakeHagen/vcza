import click
import zarr

from .write_block import write_block
from .ingest_block import ingest_block
from .dvcza import dvcza


@click.group()
def cli():
    return


@cli.command()
@click.option("-b", "--block", type=int, default=None)
@click.argument("vcz")
def write(vcz, block):
    write_block(vcz, block, None)

@cli.command()
@click.option("-i", "--infos", multiple=True)
@click.option("-b", "--block", type=int, default=None)
@click.option("-z", "--vcz", type=click.Path(exists=True))
@click.argument("vcf", default="-")
def ingest(vcf, vcz, block, infos):
    ingest_block(vcf, vcz, block, infos)

@cli.command(name="dvcza")
@click.option("-z", "--vcz", type=click.Path(exists=True))
@click.option("-c", "--cmd", type=str)
@click.option("-i", "--infos", multiple=True)
def dvcza_(vcz, cmd, infos):
    root = zarr.open(vcz, mode="r+")
    dvcza(root, cmd, infos)

cli()
