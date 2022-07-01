import pathlib
from struct import unpack
from typing import BinaryIO, List, Optional, Tuple, cast

import numpy as np
import scipy.sparse


def load(cooccurrence_file_content: BinaryIO) -> scipy.sparse.coo_matrix:
    content = file.read()
    row: List[int] = []
    column: List[int] = []
    data: List[float] = []
    for i in range(len(content/16)):
        # get the records
        relevant_bytes = content[i*16:(i+1)*16]
        crec = cast(Tuple[int, int, float], unpack('<iid', relevant_bytes))
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
