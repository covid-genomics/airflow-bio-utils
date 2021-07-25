import random
import string
from datetime import date
from math import log

random.seed(42)


def qual_to_prob(Q: float) -> float:
    return 10 ** (Q / -10)


def prob_to_qual(p: float) -> float:
    return -10 * log(p, 10)


def get_random_seq_id() -> str:
    letters = string.ascii_lowercase
    return "/".join(
        [
            "".join(
                random.choice(letters)
                for j in range(int(15 - i * random.randrange(1, 10000) / 5000))
            )
            for i in range(random.randint(1, 4))
        ]
    )


def generate_fastq_file(output_file="generated.fastq"):
    nucs = ("A", "C", "T", "G")
    with open(output_file, "w") as fastq:
        for i in range(1000):
            start_dt = date.today().replace(day=1, month=1).toordinal()
            end_dt = date.today().toordinal()
            random_day = date.fromordinal(random.randint(start_dt, end_dt))
            seq = [
                random.choice(nucs) if random.random() < 0.99 else "N"
                for _ in range(100)
            ]
            seq_name = f"{get_random_seq_id()}|seq/fastq:{len(seq)}|test/sequence/{i}|{random_day}"
            fastq.write(f">{seq_name}\n")
            fastq.write("".join(seq) + "\n")
            fastq.write("+\n")
            qual = [
                chr(int(prob_to_qual(1.0 - random.random()) + 33))
                if n in nucs
                else "!"
                for n in seq
            ]
            fastq.write("".join(qual) + "\n")
