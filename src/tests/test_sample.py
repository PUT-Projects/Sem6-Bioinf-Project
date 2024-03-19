from .. import dnaParser

def Test_EmptyDNA():
    dna = dnaParser.DNA()
    key = dna.getKey()
    length = dna.getLength()
    start = dna.getStart()
    probes = dna.getProbes()
    assert key == ""
    assert length == ""
    assert start == ""
    assert probes == []


Test_EmptyDNA()
