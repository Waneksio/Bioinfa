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
    digested = []
    for i in range(num_snips):
        snip_seed = random.randint(0,len(sequence) - oli_len)
        snippet = sequence[snip_seed:snip_seed+oli_len]
        digested.append(snippet)
    return digested
