class Sequence:
    def __init__(self, initial_oli_nucleotide):
        self.sequence = initial_oli_nucleotide
        self.nucleotides_order = [initial_oli_nucleotide]

    def get_nucleotide_value(self, nucleotide):
        n = len(self.sequence)
        k = len(nucleotide)
        for i in range(k):
            if nucleotide.find(self.sequence[n - k + i + 1:]) == 0:
                return i + 1

    def add_nucleotide(self, nucleotide):
        self.nucleotides_order.append(nucleotide)
        self.sequence += nucleotide[len(nucleotide) - self.get_nucleotide_value(nucleotide):]
        return self.sequence

    def get_len(self):
        return len(self.sequence)

    def get_nucleotides_count(self):
        return len(self.nucleotides_order)

    def rebuild_sequence(self, subsequence_a, subsequence_b):
        self.sequence = subsequence_a
        self.add_nucleotide(subsequence_b)

    def remove_nucleotide(self, nucleotide):
        last_occurrence = len(self.nucleotides_order) - self.nucleotides_order[::-1].index(nucleotide) - 1
        if last_occurrence == 0:
            return
        begin_index = self.sequence.rfind(nucleotide)
        end_index = begin_index + len(nucleotide)
        if last_occurrence == len(self.nucleotides_order) - 1:
            self.sequence = self.sequence[end_index:]
            self.sequence.add_nucleotide(self.nucleotides_order[last_occurrence - 1])
            return
        subsequence_a = Sequence(self.sequence[:begin_index])
        subsequence_a.add_nucleotide(self.nucleotides_order[last_occurrence - 1])
        subsequence_b = Sequence(self.nucleotides_order[last_occurrence + 1])
        subsequence_b.add_nucleotide(self.sequence[end_index:])
        self.rebuild_sequence(subsequence_a.sequence, subsequence_b.sequence)