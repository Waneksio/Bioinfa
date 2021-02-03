import random
from Nucleotide import Nucleotide
from Sequence import Sequence


def generate_naive_solution(oli_nucleotides, required_length, k):
    sequence = Sequence(oli_nucleotides[0])
    while sequence.get_len() < required_length:
        values = get_values_2(sequence.sequence, oli_nucleotides)
        next_oli_index = select_max(values)
        sequence.add_nucleotide(oli_nucleotides[next_oli_index])
    return sequence


def get_value(sequence, nucleotide):
    n = len(sequence)
    k = len(nucleotide)
    for i in range(k):
        if nucleotide.find(sequence[n - k + i + 1:]) == 0:
            return i + 1


def get_values_2(sequence, nucleotides):
    values = [get_value(sequence, nucleotide) for nucleotide in nucleotides]
    sequences = [add_nucleotide(sequence, nucleotides[i], values[i]) for i in range(len(values))]
    values_2 = [values[i] + min(get_values(nucleotides[i], nucleotides)) for i in range(len(sequences))]
    return values_2


def select_max(values):
    min_value = 10000000
    min_indices = []
    for i in range(len(values)):
        if values[i] < min_value:
            min_value = values[i]
            min_indices = [i]
        elif values[i] == min_value:
            min_indices.append(i)
    return random.choice(min_indices)


def get_values(sequence, nucleotides):
    values = [get_value(sequence, nucleotide) for nucleotide in nucleotides]
    return values


def add_nucleotide(sequence, nucleotide, value):
    sequence += nucleotide[len(nucleotide) - get_value(sequence, nucleotide):]
    result = sequence
    return result


def main():
    sequence = generate_naive_solution(["ATCG", "TCGG", "GATG", "GGGT", "AGTA"], 12, 4)
    print(sequence.sequence)
    print(sequence.get_remove_score(2))
    print(sequence.remove_nucleotide(2))
    print(sequence.sequence)
    return
    print(generate_naive_solution(["ATCG", "TCGG", "GATG", "GGGT", "AGTA"], 12, 4))
    return 0


if __name__ == '__main__':
    main()