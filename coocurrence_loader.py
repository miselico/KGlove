import pathlib
from struct import unpack
from typing import BinaryIO, List, Optional, Tuple, cast

import numpy as np
import scipy.sparse


def _read_little_endian_crec(file: BinaryIO
                             ) -> Optional[Tuple[int, int, float]]:
    le_int = file.read(16)
    # https://docs.python.org/3/library/struct.html#format-strings
    if len(le_int) == 0:
        return None
    crec = cast(Tuple[int, int, float], unpack('<iid', le_int))
    return crec


def load(cooccurrence_file_content: BinaryIO) -> scipy.sparse.coo_matrix:

    row: List[int] = []
    column: List[int] = []
    data: List[float] = []
    while (cooccurrence_file_content.readable()):
        crec = _read_little_endian_crec(cooccurrence_file_content)
        if crec is None:
            break
        row.append(crec[0])
        column.append(crec[1])
        data.append(crec[2])
    result = scipy.sparse.coo_matrix((data, (row, column)), dtype=np.float64)
    return result


if __name__ == "__main__":
    p = pathlib.Path("output/cooccurrence_file.bin")
    with open(p, 'rb') as file:
        m = load(file)
    print(m.tocsc())
