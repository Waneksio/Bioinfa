import random
import numpy as np
from Nucleotide import Nucleotide
from Sequence import Sequence
import Generator as gen


def generate_naive_solution(oli_nucleotides, required_length, k):
    sequence = Sequence(oli_nucleotides[0])
    sequence.initialize_tabu_list(oli_nucleotides)
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

    for t1 in range(n + 1):
        for t2 in range(m + 1):
            print(int(distances[t1][t2]), end=" ")
        print()
    return distances[n][m]




def main():
    sequence = generate_naive_solution(["ATCG", "TCGG", "GATG", "GGGT", "AGTA"], 25, 4)
#    print(sequence.sequence)
#    print(sequence.nucleotides_order)
#    sequence.cluster()
    #sequence.condensation()
#    print(sequence.sequence)
#    random = gen.generate_random_test_sequence(50)
#    print(random)
#    digested = gen.snip_snap(4, 20,random)
    distance = get_Lev_distance('kelm','hello')
    print(int(distance))
#    for i in digested:
 #       print(i)
    return
    sequence.condensation()
#    print(sequence.sequence)

    return 0


if __name__ == '__main__':
    main()