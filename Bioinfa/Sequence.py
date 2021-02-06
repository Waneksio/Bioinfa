import random


class Sequence:
    def __init__(self, initial_oli_nucleotide, desired_length, oligos):
        self.desired_length = desired_length
        self.sequence = initial_oli_nucleotide
        self.nucleotides_order = [[initial_oli_nucleotide, 0]]
        self.tabu_list = []
        self.all_oligos = oligos
        self.tabu_length = int(len(oligos) / 10)

    def get_nucleotide_value(self, nucleotide):
        n = len(self.sequence)
        k = len(nucleotide)
        for i in range(k):
            if nucleotide.find(self.sequence[n - k + i + 1:]) == 0:
                return i + 1
        return k

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
        length = self.get_len()
        self.sequence = subsequence_a
        value_b = self.get_nucleotide_value(subsequence_b[:self.get_nucleotide_length()])
        self.sequence += subsequence_b[self.get_nucleotide_length() - value_b:]
        offset = length - self.get_len()
        return offset

    @staticmethod
    def get_common_nucleotides(sequence_a, sequence_b):
        n = len(sequence_a)
        k = len(sequence_b)
        for i in range(k):
            if sequence_b.find(sequence_a[n - k + i + 1:]) == 0:
                return k - (i + 1)

    def get_current_condensation_score(self):
        return self.get_nucleotides_count() / self.get_len()

    def get_nucleotide_removal_score(self, index):
        if index == 0:
            return 0

        current_nucleotide = self.nucleotides_order[index][0]
        previous_nucleotide = self.nucleotides_order[index - 1][0]
        oligonucleotides_in_solution = self.get_nucleotides_count() - 1

        if index + 1 == len(self.nucleotides_order):
            reduction = max(0, len(current_nucleotide) - self.get_common_nucleotides(previous_nucleotide, current_nucleotide))
            nucleotides_in_solution = self.get_len() - reduction
            return oligonucleotides_in_solution / nucleotides_in_solution

        next_nucleotide = self.nucleotides_order[index + 1][0]
        reduction = max(0, len(current_nucleotide) - self.get_common_nucleotides(previous_nucleotide, current_nucleotide) - self.get_common_nucleotides(current_nucleotide, next_nucleotide))

        nucleotides_in_solution = self.get_len() - reduction
        return oligonucleotides_in_solution / nucleotides_in_solution

    def get_element_to_remove(self):
        result = max(range(len(self.nucleotides_order)), key=lambda x: self.get_nucleotide_removal_score(x))
        if self.get_nucleotide_removal_score(result) >= self.get_current_condensation_score():
            return result
        return None

    def get_least_frequent_element(self):
        result = random.choice(self.all_oligos)
        while result in self.tabu_list:
            result = random.choice(self.all_oligos)
        return result

    def extension(self):
        while self.get_len() < self.desired_length:
            least_frequent_element = self.get_least_frequent_element()
            self.add_nucleotide(least_frequent_element)
            self.tabu_list.append(least_frequent_element)
            if len(self.tabu_list) > self.tabu_length:
                self.tabu_list.pop(0)

    def condensation(self):
        element_to_remove = self.get_element_to_remove()
        while element_to_remove:
            self.tabu_list.append(self.nucleotides_order[element_to_remove][0])
            if len(self.tabu_list) > self.tabu_length:
                self.tabu_list.pop(0)
            self.remove_nucleotide(element_to_remove)
            element_to_remove = self.get_element_to_remove()

    def remove_nucleotide(self, index):
        previous_nucleotide = self.nucleotides_order[index - 1]

        if index + 1 == self.get_nucleotides_count():
            current_nucleotide = self.nucleotides_order[index]
            self.sequence = self.sequence[:-(len(current_nucleotide[0]) - self.get_common_nucleotides(previous_nucleotide[0], current_nucleotide[0]))]
            self.nucleotides_order.pop(index)
            return

        next_nucleotide = self.nucleotides_order[index + 1]
        sub_sequence_a = self.sequence[:previous_nucleotide[1] + 4]
        sub_sequence_b = self.sequence[next_nucleotide[1]:]
        self.nucleotides_order.pop(index)
        offset = self.rebuild_sequence(sub_sequence_a, sub_sequence_b)

        for i in range(self.get_nucleotides_count()):
            if i >= index:
                self.nucleotides_order[i][1] -= offset

    def cluster(self):
        clustered_sequence = []
        sequence_copy = self.nucleotides_order[:]
        for i in range(len(sequence_copy)):
            if i+1 >= len(sequence_copy):
                if i==len(sequence_copy):
                    return clustered_sequence
                clustered_sequence.append(sequence_copy[i])
                return clustered_sequence
            if sequence_copy[i+1][1] - sequence_copy[i][1] == 1:
                clustered_nucleotide = sequence_copy[i]
                clustered_nucleotide[0] = sequence_copy[i][0] + sequence_copy[i+1][0][-1]
                clustered_sequence.append(clustered_nucleotide)
                del sequence_copy[i+1]
            else:
                clustered_sequence.append(sequence_copy[i])

    def remove_oli_at_position(self, position):
        for i in self.nucleotides_order:
            if i[1] == position:
                del i
                return