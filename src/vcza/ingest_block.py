import numpy as np
import numcodecs
from cyvcf2 import VCF

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


def create_info_array(root, info_name, dtype, fill_value, desc=None):
    array_name = f"variant_{info_name}"
    if array_name in root:
        raise Exception(
            f"{array_name} already exists in this vcz, delete it first if it needs to be re-done"
        )

    if dtype == "|O":
        object_codec = numcodecs.VLenUTF8()
    else:
        object_codec = None

    root.full(
        name=array_name,
        shape=root["variant_position"].shape,
        chunks=root["variant_position"].chunks,
        dtype=dtype,
        fill_value=fill_value,
        object_codec=object_codec,
    )

    if desc:
        root[array_name].attrs.update(
            {
                "description": desc,
                "_ARRAY_DIMENSIONS": root["variant_position"].attrs[
                    "_ARRAY_DIMENSIONS"
                ],
            }
        )
    else:
        root[array_name].attrs.update(
            {
                "description": "added by vcza",
                "_ARRAY_DIMENSIONS": root["variant_position"].attrs[
                    "_ARRAY_DIMENSIONS"
                ],
            }
        )


def ingest_block(in_, root, block, infos):
    d_ = {}
    for field, dtype, fill in parse_infos(infos):
        d_[field] = np.full(
            root["variant_position"].blocks[block].shape, fill, dtype=dtype
        )

    vcf = VCF(in_)

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
