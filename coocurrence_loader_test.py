import coocurrence_loader
from io import BytesIO


class TestCooccurrenceLoader:

    def test_read_crec(self) -> None:
        bio = BytesIO(b'\x03\x00\x00\x00\x02\x00\x00\x00ffffff\xe6?')
        crec = coocurrence_loader._read_little_endian_crec(bio)
        assert crec is not None
        assert crec[0] == 3
        assert crec[1] == 2
        assert crec[2] == 0.7
