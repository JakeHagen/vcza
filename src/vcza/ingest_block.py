import fasteners
import numpy as np
import numcodecs
from cyvcf2 import VCF
import zarr

from collections.abc import Iterator

TYPE_FILL_DICT = {
    "Integer": ("<i4", -2),
    "Float": (
        "<f4",
        np.array([0x7F800001, 0x7F800002], dtype=np.int32).view(np.float32)[1],
    ),
    "Flag": (np.bool, False),
    "String": ("|O", ""),
}


def parse_infos(infos: tuple[str]) -> Iterator[tuple[str, str, str]]:
    for info in infos:
        field, type_ = info.split(":")
        yield (field, TYPE_FILL_DICT[type_][0], TYPE_FILL_DICT[type_][1])


def create_info_array(root, info_name, lock_dir, dtype, fill_value, desc):
    array_name = f"variant_{info_name}"

    # create lock so we only make one array
    lock = fasteners.InterProcessLock(f'{lock_dir}/{array_name}.lock_file')
    lock.acquire()

    # recheck if array exists, other process might have gotten past inital check but should hang on acquiring lock
    # once process acquires lock, array should exi
    if array_name in root:
        lock.release()
        return

    root.full(
        name=array_name,
        shape=root["variant_position"].shape,
        chunks=root["variant_position"].chunks,
        dtype=dtype,
        fill_value=fill_value,
        object_codec=numcodecs.VLenUTF8() if dtype == "|O" else None
    )

    root[array_name].attrs.update(
        {
            "description": desc,
            "_ARRAY_DIMENSIONS": root["variant_position"].attrs[
                "_ARRAY_DIMENSIONS"
            ]
        }
    )

    lock.release()
    return

def ingest_block(in_, vcz, block, infos):
    d_ = {}
    root = zarr.open(vcz, mode="r+")

    vcf = VCF(in_)
    hdr_infos = {}
    for h_ in vcf.header_iter():
       if h_['HeaderType'] == 'INFO':
           hdr_infos[h_['ID']] = h_

    for field in infos:
        array_name = f"variant_{field}"

        try:
            hdr_info = hdr_infos[field]
        except KeyError as e:
            print(f'{field} is not in vcfs header')
            raise e

        desc = hdr_info['Description']
        dtype = TYPE_FILL_DICT[hdr_info['Type']][0]
        fill = TYPE_FILL_DICT[hdr_info['Type']][1]

        d_[field] = np.full(
            shape=root["variant_position"].blocks[block].shape, fill_value=fill, dtype=dtype
        )

        if array_name not in root:
            create_info_array(root, field, vcz, dtype, fill, desc)



    for idx, v in enumerate(vcf):
        for field in d_:
            # if the variant doesnt have the info field, continue (keep array default value)
            try:
                field_val = v.INFO[field]
            except KeyError:
                continue
            # if field value is none ie "." (None from cyvcf2), continue (keep the array fill value)
            try:
                d_[field][idx] = field_val
            except TypeError:
                if not field_val:
                    continue

    for field in d_:
        root[f"variant_{field}"].blocks[block] = d_[field]

    return
