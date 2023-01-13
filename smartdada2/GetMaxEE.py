import numpy as np
import pandas as pd

import smartdada2.GetTrimParameters as GTP
from smartdada2.reader import reader


def read_size_by_maxEE(FastqEntries, left: int, right: int):
    """Takes reads and calculates the sum of the expected
    error for each read before any trimming and after
    obvious trimming.

    Args:
        FastqEntries (FastqEntry): reads taken from fastq files
        left (int): trim value output from trim_ends_less_than_threshold
        right (int): trunc value output from trim_ends_less_than_threshold

    Returns:
        pd.dataframe: dataframe containing two columns
    """
    # raise errors for wronog data types
    if not isinstance(FastqEntries, reader.FastqReader):
        raise TypeError("must be of FastqEntry type")

    # get list of expected errors per read
    raw_EE = []

    # for each entry
    for entry in FastqEntries.unpack_entries():
        phred_scores = list(entry.scores)

        # convert scores to numeric format to integer 0-42
        scores = np.array(phred_scores).view(np.int32) - 33

        # converts phred scores to expected error
        EE_entry = [10 ** (-score / 10) for score in scores]

        # create list of lists of expected error
        raw_EE.append(EE_entry)

    # get sum of EE for no trimming and obvious trimming
    no_trimming = []
    obv_trimming = []

    # get the sum of EE
    for EE_list in raw_EE:
        obv_trim = EE_list[left:right]
        sum_no_trim = sum(EE_list)
        sum_obv_trim = sum(obv_trim)

        # append to lists
        no_trimming.append(sum_no_trim)
        obv_trimming.append(sum_obv_trim)

    # create two column dataframe
    sumEE = pd.DataFrame(
        list(zip(no_trimming, obv_trimming)), columns=["NoTrimming", "ObviousTrimming"]
    )

    sumEE = sumEE.reset_index(drop=True)

    return sumEE
