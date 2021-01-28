class Nucleotide:
    def __init__(self, sequence, next_element=None, previous_element=None):
        self.sequence = sequence
        self.next = next_element
        self.previous = previous_element

    def append(self, nucleotide):
        last_nucleotide = self.get_last()
        last_nucleotide.next = nucleotide
        nucleotide.previous = last_nucleotide

    def get_last(self):
        last_nucleotide = self
        while last_nucleotide.next:
            last_nucleotide = last_nucleotide.next
        return last_nucleotide

    def remove(self):
        pass

    def print(self):
        nucleotide = self
        while nucleotide.next:
            print(nucleotide.sequence)
            nucleotide = nucleotide.next
        print(nucleotide.sequence)