
class HITRANTIPS:

    def __init__(self, molid: int, isoid: int):
        self._molid = molid
        self._isoid = isoid

    
    def Q(self, temperature):
        from .hapi import partitionSum
        return partitionSum(self._molid, self._isoid, temperature)