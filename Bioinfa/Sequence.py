class Sequence:
    def __init__(self, initial_oli_nucleotide):
        self.sequence = initial_oli_nucleotide
        self.nucleotides_order = [[initial_oli_nucleotide, 0]]

    def get_nucleotide_value(self, nucleotide):
        n = len(self.sequence)
        k = len(nucleotide)
        for i in range(k):
            if nucleotide.find(self.sequence[n - k + i + 1:]) == 0:
                return i + 1

    def add_nucleotide(self, nucleotide):
        self.sequence += nucleotide[len(nucleotide) - self.get_nucleotide_value(nucleotide):]
        size = self.get_len()
        self.nucleotides_order.append([nucleotide, size - len(nucleotide)])
        return self.sequence

    def get_len(self):
        return len(self.sequence)

    def get_nucleotides_count(self):
        return len(self.nucleotides_order)

    def get_nucleotide_length(self):
        return len(self.nucleotides_order[0][0])

    def rebuild_sequence(self, subsequence_a, subsequence_b):
        self.sequence = subsequence_a
        value_b = self.get_nucleotide_value(subsequence_b[:self.get_nucleotide_length()])
        self.sequence += subsequence_b[self.get_nucleotide_length() - value_b:]
        offset = self.get_len() - len(subsequence_b)
        return offset

    @staticmethod
    def get_common_nucleotides(sequence_a, sequence_b):
        n = len(sequence_a)
        k = len(sequence_b)
        for i in range(k):
            if sequence_a.find(sequence_b[n - k + i + 1:]) == 0:
                return k - (i + 1)

    def get_remove_score(self, index):
        current_nucleotide = self.nucleotides_order[index][0]
        previous_nucleotide = self.nucleotides_order[index - 1][0]
        next_nucleotide = self.nucleotides_order[index + 1][0]
        return max(0, len(current_nucleotide) - self.get_common_nucleotides(previous_nucleotide, current_nucleotide) - self.get_common_nucleotides(current_nucleotide, next_nucleotide))

    # TODO: unify get_nucleotide_value and get_remove_score functions
    # already done we can say

    def remove_nucleotide(self, index):
        initial_length = self.get_len()
        previous_nucleotide = self.nucleotides_order[index - 1]
        next_nucleotide = self.nucleotides_order[index + 1]
        sub_sequence_a = self.sequence[:previous_nucleotide[1] + 4]
        sub_sequence_b = self.sequence[next_nucleotide[1]:]
        self.nucleotides_order.pop(index)
        offset = self.rebuild_sequence(sub_sequence_a, sub_sequence_b)
        for nucleotide in self.nucleotides_order:
            nucleotide[1] += offset
        return self.get_len() - initial_length