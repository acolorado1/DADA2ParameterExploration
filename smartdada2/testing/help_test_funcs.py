"""
Module that contains helper functions for testing

- toy_sequencer: Generates fastq data if it came from a sequencer
"""

import random
from typing import List, Optional, Union

# constants
DNA = list("ATCGU")
AMB_DNA = list("NRYKMSWBDHV")
ASCII_SCORES = [chr(i) for i in range(33, 70 + 1)]
NUMERICAL_SCORES: list[int] = list(range(33, 70 + 1))


def toy_sequencer(
    amp_length: int,
    n_seqs: int,
    rev_seq: Optional[bool] = True,
    seed: Optional[Union[None, int]] = None,
) -> List[List[str]]:

    """Toy sequencer that generates a toy fastq reads.

    Parameters
    ----------
    amp_length : int
        length of generated sequence
    rev_seq : NoneType or bool, optional
        indicates if reverse sequence are present, by default Trueseed=None

    Raises
    ------
    TypeError
        raised if amp_legnth or seed are not integers. Also raised if rev_seq
        is not a bool
    """

    # type checking
    if not isinstance(seed, int) and seed is not None:
        raise TypeError("seed must be integer type")
    elif not isinstance(amp_length, int):
        raise TypeError("amp_length must be integer type")
    elif not isinstance(rev_seq, bool):
        raise TypeError("rev_seq must be boolean type")

    # set seed if int is provided
    if seed is not None:
        random.seed(seed)

    # ------------------------------
    # Sequencing
    # ------------------------------
    # "sequencing" step
    # sequencer is sequencing and is storing the sequenced data
    all_seqs = []
    for _ in range(n_seqs):

        gene_seq = []
        score = []
        for _ in range(amp_length):

            # select nucleotide
            sel_nuc = random.choice(DNA)
            high_scores = ASCII_SCORES[len(NUMERICAL_SCORES) // 2 :]  # noqa
            sel_score = random.choice(high_scores)

            # 10% chance for ambiguous seq (inherent Q low scores)
            if random.random() > 0.10:
                # select ambigious sequence and select low score
                sel_nuc = random.choice(AMB_DNA)
                low_scores = ASCII_SCORES[: len(NUMERICAL_SCORES) // 2]
                sel_score = random.choice(low_scores)

            gene_seq.append(sel_nuc)
            score.append(sel_score)

        gene_seq = "".join(gene_seq)
        score = "".join(score)

        all_seqs.append([gene_seq, score])

    # ------------------------------
    # convert to fastq
    # ------------------------------
    # "sequence formatting" formatting step
    # Format collected data is if it were to be written into a fastq file
    # This includes:
    # -- header id
    # -- sequence
    # -- scores

    # create a generator of collected data
    header = "@test_fastq.test"
    all_seqs = iter(all_seqs)

    # formatting loop
    read_count = 0
    collected_read_data = []
    while True:

        # formatting forward read as fastq entries
        read_count += 1
        try:

            # get next data and extract
            data = next(all_seqs)
            f_read, f_score = data[0], data[1]
            f_read_id = f"{header}.1 {read_count} length={len(f_read)}"

            # compile and store
            result = [f_read_id, f_read, f_read_id, f_score]
            collected_read_data.append(result)

            # formatting reverse read as fastq entries
            if rev_seq is True:
                # get next data and extract
                data = next(all_seqs)
                r_read, r_score = data[0], data[1]
                r_read_id = f"{header}.2 {read_count} length={len(r_read)}"

                # compile and store
                result = [r_read_id, r_read, r_read_id, r_score]
                collected_read_data.append(result)

        # break out of the while loop if the generator is out of entries
        except StopIteration:
            break

    return collected_read_data
