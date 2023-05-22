import re


def find_orfs(seq: str) -> []:
    start_indices = __find_start_indices(seq)
    end_indices = __find_end_indices(seq)
    if not __number_of_indices_is_equal(start_indices, end_indices):
        end_indices = __remove_non_matching_indices(start_indices, end_indices)
    return __merge(start_indices, end_indices)


def __find_start_indices(seq: str) -> []:
    return [m.start() + 1 for m in re.finditer('ATG', seq)]


def __find_end_indices(seq: str) -> []:
    return [m.start() + 1 for m in re.finditer('TAA', seq)] + \
        [m.start() + 1 for m in re.finditer('TAG', seq)] + \
        [m.start() + 1 for m in re.finditer('TGA', seq)]


def __number_of_indices_is_equal(a: [], b: []) -> bool:
    return len(a) == len(b)


def __remove_non_matching_indices(starts: [], ends: []) -> []:
    for start in starts:
        for end in ends:
            if not __indices_match(start, end):
                ends.remove(end)
        break
    return ends


def __indices_match(start: int, end: int) -> bool:
    return (end - start) % 3 == 0


def __merge(start: [], end: []) -> []:
    result = []
    for i, v in enumerate(start):
        result.append((start[i], end[i] + 2))
    return result


class DNAGenome:
    def __init__(self, seq):
        self.seq = seq
        self.trans_dict = self.__get_trans_dict()

    def get_cds(self) -> []:
        orfs = find_orfs(self.seq)
        cds = []
        for orf in orfs:
            stripped = self.__strip_sequence(orf)
            cds.append(self.__split_sequence(stripped))
        return cds

    def get_translations(self) -> []:
        cds = self.get_cds()
        translations = []
        for cd in cds:
            acids = ''
            for i in cd:
                acids = acids + self.__translate_cd_to_aminoacid(i)
            translations.append(acids)
        return translations

    @staticmethod
    def __get_trans_dict() -> {}:
        bases = "TCAG"
        codons = [a + b + c for a in bases for b in bases for c in bases]
        aminoacids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        return dict(zip(codons, aminoacids))

    def __strip_sequence(self, orf: []) -> str:
        return self.seq[orf[0] - 1:orf[1]]

    @staticmethod
    def __split_sequence(stripped: str) -> []:
        return [stripped[i:i + 3] for i in range(0, len(stripped), 3)]

    def __translate_cd_to_aminoacid(self, cd: str) -> str:
        return self.trans_dict.get(cd)
