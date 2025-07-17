from vcztools import _vcztools
import numpy as np


def write_hdr(file):
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
    print("##fileformat=VCFv4.4", file=file)
    for chrom, length in contigs:
        print(f"##contig=<ID={chrom},length={length}>", file=file)
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
        file=file,
    )


def write_block(root, block, output):
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

    write_hdr(output)

    for j in range(len(root["variant_position"].blocks[block])):
        buflen = 1024
        while True:
            try:
                line = encoder.encode(j, buflen)
                print(line, file=output)
                break
            except _vcztools.VczBufferTooSmall:
                buflen *= 2
