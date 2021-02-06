import random
import numpy as np
from Nucleotide import Nucleotide
from Sequence import Sequence
import Generator as gen
import copy
import time


def check_length(sequences, required_length):
    for sequence in sequences:
        if sequence.get_len() < required_length:
            return True
    return False


def generate_naive_solution(oli_nucleotides, required_length, k, max_solutions):
    sequences = []
    sequence = Sequence(oli_nucleotides[0], required_length, oli_nucleotides)
    sequences.append(sequence)
    solutions_left = max_solutions

    while check_length(sequences, required_length):
        for sequence1 in sequences:
            if sequence1.get_len() >= required_length:
                continue
            values = get_values_2(sequence1.sequence, oli_nucleotides)
            potential_next = select_max(values)

            for i in range(1, len(potential_next)):
                if solutions_left > 0:
                    sequences.append(copy.deepcopy(sequence1))
                    solutions_left -= 1
                    sequences[max_solutions - solutions_left].add_nucleotide(oli_nucleotides[potential_next[i]])

            sequence1.add_nucleotide(oli_nucleotides[potential_next[0]])
    return sequences


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
    return min_indices


def get_values(sequence, nucleotides):
    values = [get_value(sequence, nucleotide) for nucleotide in nucleotides]
    return values


def add_nucleotide(sequence, nucleotide, value):
    sequence += nucleotide[len(nucleotide) - get_value(sequence, nucleotide):]
    result = sequence
    return result


def get_Lev_distance(real_sequence, solution_sequence):
    if real_sequence == '' or real_sequence == None:
        if solution_sequence == ''or solution_sequence == None:
            return
        return len(solution_sequence)
    if solution_sequence == ''or solution_sequence == None:
        return len(real_sequence)

    n = len(real_sequence)
    m = len(solution_sequence)
    distances = np.zeros((n + 1, m + 1))
    for t1 in range(n+1):
        distances[t1][0] = t1
    for t2 in range(m+1):
        distances[0][t2] = t2

    a = 0
    b = 0
    c = 0

    for t1 in range(1, n + 1):
        for t2 in range(1, m + 1):
            if (real_sequence[t1 - 1] == solution_sequence[t2 - 1]):
                distances[t1][t2] = distances[t1 - 1][t2 - 1]
            else:
                a = distances[t1][t2 - 1]
                b = distances[t1 - 1][t2]
                c = distances[t1 - 1][t2 - 1]

                if (a <= b and a <= c):
                    distances[t1][t2] = a + 1
                elif (b <= a and b <= c):
                    distances[t1][t2] = b + 1
                else:
                    distances[t1][t2] = c + 1
    return distances[n][m]


def main():
    k = 8
    m = 190
    n = 200
    radom = gen.generate_random_test_sequence(n)
    digested = gen.snip_snap(k, m, radom)

    start_time = time.time()
    sequences = generate_naive_solution(digested, n, k, 10)

    absolutely_best_sequence = sequences[0]
    for sequence in sequences:
        best_sequence = sequence
        for _ in range(10):
            if get_Lev_distance(best_sequence.sequence, radom) > get_Lev_distance(sequence.sequence, radom):
                best_sequence = sequence.sequence
            sequence.condensation()
            sequence.extension()
        if get_Lev_distance(best_sequence.sequence, radom) < get_Lev_distance(absolutely_best_sequence.sequence, radom):
            absolutely_best_sequence = best_sequence
    stop_time = time.time()
    print("time: " + str(stop_time - start_time))
    print(absolutely_best_sequence.sequence)
    print(radom)
    print(get_Lev_distance(absolutely_best_sequence.sequence, radom))
    return
    print(get_Lev_distance(sequence.sequence, radom))
    sequence.condensation()
#    print(sequence.sequence)

    return 0


if __name__ == '__main__':
    main()