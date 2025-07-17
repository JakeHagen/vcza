from dask_jobqueue.slurm import SLURMCluster
from dask.distributed import Client
import numpy as np
from cyvcf2 import VCF
import subprocess
from vcztools import _vcztools

from .write_block import write_block
from .ingest_block import ingest_block, create_info_array, parse_infos

def get_process(cmd):
    return subprocess.Popen(
        cmd.split(" "), stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
    )

#def get_func(root, cmd, info_fields):
#    def f(block):
#        proc = get_process(cmd)
#        write_block(root, block, proc.stdin)
#        proc.stdin.close()

#        ingest_block(proc.stdout, root, block, info_fields)

#    return f
#
#def f(block, root=None, cmd=None, infos=None):
#    proc = get_process(cmd)
#    write_block(root, block, proc.stdin)
#    proc.stdin.close()
#    ingest_block(proc.stdout, root, block, infos)

def dvcza(root, cmd, infos):
    cluster = SLURMCluster(
        cores=6,
        processes=1,
        memory="30GB",
        account="shenlab",
        walltime="01:00:00",
    )
    cluster.adapt(minimum=1, maximum=100)

    client = Client(cluster)

    for field, dtype, fill in parse_infos(infos):
        create_info_array(root, field, dtype, fill)

    nblocks = root["variant_position"].nchunks


    def f(block, cmd=None, root=None, infos=None):
        proc = subprocess.Popen(cmd.split(' '), stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        contigs = [
            ("chr1", "248956422"),
            ("chr2", "242193529"),
            ("chr3", "198295559"),
            ("chr4", "190214555"),
            ("chr5", "181538259"),
            ("chr6", "170805979"),
            ("chr7", "159345973"),
            ("chr8", "145138636"),
            ("chr9", "138394717"),
            ("chr10", "133797422"),
            ("chr11", "135086622"),
            ("chr12", "133275309"),
            ("chr13", "114364328"),
            ("chr14", "107043718"),
            ("chr15", "101991189"),
            ("chr16", "90338345"),
            ("chr17", "83257441"),
            ("chr18", "80373285"),
            ("chr19", "58617616"),
            ("chr20", "64444167"),
            ("chr21", "46709983"),
            ("chr22", "50818468"),
            ("chrX", "156040895"),
            ("chrY", "57227415"),
            ("chrM", "16569"),
        ]
        print("##fileformat=VCFv4.4", file=proc.stdin)
        for chrom, length in contigs:
            print(f"##contig=<ID={chrom},length={length}>", file=proc.stdin)
        print(
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            sep="\t",
            file=proc.stdin,
        )
        encoder = _vcztools.VcfEncoder(
            len(root["variant_position"].blocks[block]),
            0,
            id=root["variant_id"].blocks[block].astype("S").reshape((-1, 1)),
            qual=root["variant_quality"].blocks[block],
            filter_ids=np.array([b"PASS"]).astype("S"),
            filter=np.full((len(root["variant_position"].blocks[block]), 1), False),
            chrom=root["contig_id"][:].astype("S")[root["variant_contig"].blocks[block]],
            pos=root["variant_position"].blocks[block].astype(np.int32),
            alt=root["variant_allele"].blocks[block][:, 1:].astype("S"),
            ref=root["variant_allele"].blocks[block][:, 0].astype("S"),
        )
        for j in range(len(root["variant_position"].blocks[block])):
            buflen = 1024
            while True:
                try:
                    line = encoder.encode(j, buflen)
                    print(line, file=proc.stdin)
                    break
                except _vcztools.VczBufferTooSmall:
                    buflen *= 2
        proc.stdin.close()

        TYPE_FILL_DICT = {
            "Integer": ("<i4", -2),
            "Float": (
                "<f4",
                np.array([0x7F800001, 0x7F800002], dtype=np.int32).view(np.float32)[1],
            ),
            "Flag": (np.bool, False),
            "String": ("|O", ""),
        }

        d_ = {}
        for info in infos:
            field, type_ = info.split(":")
            field, dtype, fill = (field, TYPE_FILL_DICT[type_][0], TYPE_FILL_DICT[type_][1])

            d_[field] = np.full(
                root["variant_position"].blocks[block].shape, fill, dtype=dtype
            )

        vcf = VCF(proc.stdout)

        for idx, v in enumerate(vcf):
            for field in d_:
                d_[field][idx] = v.INFO[field]

        for field in d_:
            root[f"variant_{field}"].blocks[block] = d_[field]

        return


    #def f(block, root=None, cmd=None, infos=None):
    #    proc = get_process(cmd)
    #    write_block(root, block, proc.stdin)
    #    proc.stdin.close()
    #    ingest_block(proc.stdout, root, block, infos)

    futures = client.map(
        f, range(nblocks), cmd=cmd, root=root, infos=infos
    )


    print(futures)

    _ = client.gather(futures)
    client.close()
    cluster.close()
