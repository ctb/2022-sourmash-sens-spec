import random
import screed


nucl = ['a', 'c', 'g', 't']


def generate_mutated_reads(records, readlen, error_rate, nucl=tuple(nucl)):
    records = [ r for r in records if len(r.sequence) >= readlen ]
    len_genome = sum([ len(x.sequence) for x in records ])

    target = []
    for idx, r in enumerate(records):
        target += [ idx ] * len(r.sequence)

    while 1:
        contig_num = random.choice(target)
        contig = records[contig_num].sequence

        start = random.randint(0, len(contig) - readlen)
        read = contig[start:start + readlen].upper()

        # reverse complement?
        if random.choice([0, 1]) == 0:
            read = screed.rc(read)

        # error?
        read_mutations = 0
        for _ in range(readlen):
            if error_rate > 0:
                while random.randint(1, int(1.0/error_rate)) == 1:
                    pos = random.randint(1, readlen) - 1
                    orig = read[pos]
                    new_base = random.choice(nucl)
                    if orig.lower() == new_base:
                        continue

                    read = read[:pos] + random.choice(nucl) + read[pos+1:]
                    read_mutations += 1


        assert len(read) == readlen
        yield start, read, read_mutations
        
        
