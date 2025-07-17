import click
import zarr

from .write_block import write_block
from .ingest_block import ingest_block, create_info_array, parse_infos
from .dvcza import dvcza


@click.group()
def cli():
    return


@cli.command()
@click.option("-z", "--vcz")
@click.option("-b", "--block", type=int, default=None)
def write(vcz, block):
    root = zarr.open(vcz, mode="r+")
    write_block(root, block, None)

@cli.command()
@click.option("-i", "--infos", multiple=True)
@click.option("-z", "--vcz", type=click.Path(exists=True))
def initinfos(vcz, infos):
    root = zarr.open(vcz, mode="r+")
    for field, dtype, fill in parse_infos(infos):
        create_info_array(root, field, dtype, fill)

@cli.command()
@click.option("-i", "--infos", multiple=True)
@click.option("-b", "--block", type=int, default=None)
@click.option("-z", "--vcz", type=click.Path(exists=True))
@click.argument("vcf", default="-")
def ingest(vcf, vcz, block, infos):
    root = zarr.open(vcz, mode="r+")
    ingest_block(vcf, root, block, infos)

@cli.command(name="dvcza")
@click.option("-z", "--vcz", type=click.Path(exists=True))
@click.option("-c", "--cmd", type=str)
@click.option("-i", "--infos", multiple=True)
def dvcza_(vcz, cmd, infos):
    root = zarr.open(vcz, mode="r+")
    dvcza(root, cmd, infos)

cli()
