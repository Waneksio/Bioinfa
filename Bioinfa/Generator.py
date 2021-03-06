import random


def generate_random_test_sequence(length):
    result = ''
    for i in range(length):
        nucleo = random.randint(0,3)
        if nucleo ==0:
            result = result + "A"
        elif nucleo == 1:
            result = result + 'C'
        elif nucleo == 2:
            result = result + 'T'
        else:
            result = result + 'G'
    return result


def snip_snap(oli_len, num_snips, sequence):
    digested = [sequence[:oli_len]]
    while len(digested) < num_snips:
        snip_seed = random.randint(0,len(sequence) - oli_len)
        snippet = sequence[snip_seed:snip_seed+oli_len]
        if snippet not in digested:
            digested.append(snippet)
    with open('./test_instance.txt', "w") as f:
        f.write(sequence)
        for snip in digested:
            f.write("\n" + snip)


