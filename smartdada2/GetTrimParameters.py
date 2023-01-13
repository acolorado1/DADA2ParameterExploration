import warnings

import pandas as pd

from smartdada2.reader import reader


def trim_ends_less_than_threshold(avg_qs_list, threshold=30, max_trim_per=0.1):
    """Perform obvious trimming: moving from the middle index to the left
    and then to the right, the first index with a value lower than the
    threshold will be the potential trim index

    :param avg_qs_list: list of floats
        list of average quality scores per position
    :param threshold: integer
        quality score between 0 and 42
    :param max_trim_per: float
        maximum percentage of trimming performed at either end
    :return: return tuple of two integers
        potential left and right trim values
    """
    if not isinstance(avg_qs_list, list):
        raise TypeError("average list of quality scores must be of type list")
    if not isinstance(threshold, int):
        raise TypeError("threshold must be of type integer")
    if not isinstance(max_trim_per, float):
        raise TypeError("max trim percentage must be of type float")
    if threshold < 0 or threshold > 42:
        raise ValueError("threshold must be between 0 and 42")
    for avg_score in avg_qs_list:
        if not isinstance(avg_score, float):
            raise TypeError("quality scores must be of type float")
        if avg_score < 0 or avg_score > 42:
            raise ValueError("average scores must be between 0 and 42")

    # find the middle index
    len_list = len(avg_qs_list)
    mid_index = len_list // 2

    trim_left_index = 0
    trim_right_index = 0

    # get left index before value is below threshold
    current_index = mid_index

    # get left index before value is below threshold
    while current_index >= 0:
        # if value at current index is below threshold
        if avg_qs_list[current_index] < threshold:
            # get the prior index
            trim_left_index = current_index + 1
            break
        else:
            current_index -= 1

    # get right index before value is below threshold
    current_index = mid_index
    while current_index <= len_list and current_index < len_list:
        # if value at current index is below threshold
        if avg_qs_list[current_index] < threshold:
            # get prior checked index
            trim_right_index = current_index
            break
        else:
            current_index += 1

    # if not obvious trimming was performed at all
    if trim_left_index == 0 and trim_right_index == 0:
        warnings.warn("No obvious trimming performed")

    # if no obvious right trimming was performed
    if trim_right_index == 0:
        # right index is last index
        trim_right_index = len_list

    # warn if trim values create a short read
    perc_10 = round(len_list * max_trim_per)
    if trim_left_index > perc_10:
        warnings.warn("trim left value might be too high")

    if trim_right_index < len_list - perc_10:
        warnings.warn("trim right value might be too low")

    return trim_left_index, trim_right_index


def get_trim_length_avgEE(avgEE_list, left, right):
    """This function takes a list of average expected error
    per position and indexing sites to return a list of trim
    sites, read length, and average EE

    Args:
        avgEE_list (list): average expected error per position
        left (int): left index
        right (int): right index

    Returns:
        list: string of indexes, length of read, average EE per position
    """
    if not isinstance(avgEE_list, list):
        raise TypeError("average EE must be of type list")
    if not isinstance(left, int):
        raise TypeError("left index must be of type int")
    if not isinstance(right, int):
        raise TypeError("right index must be of type index")
    for EE in avgEE_list:
        if not isinstance(EE, float):
            raise TypeError("all values in avgEE must be of type float")

    # to get average EE per position of read
    current_read_len = len(avgEE_list[left:right])
    current_sumEE = sum(avgEE_list[left:right])

    # get indexes, read length, and avgEE
    read_len = current_read_len
    return [left, right, read_len, current_sumEE / read_len]


def read_size_by_avg_EE(FastqEntries, left: int, right: int, max_trim_perc=0.20):
    """Creates a dataframe containing average expected error per position
    for each length of read calculated within certain ranges of positions
    (e.g., not taking out more than 20% off of either end)

    Args:
        FastqEntries (FastqEntry): reads taken from fastq files
        left (integer): left index
        right (integer): right index
        max_trim_perc (float, optional): max percentage of the read trimmed on
        either end. Defaults to 0.2.

    Returns:
        pandas dataframe: dataframe containing 3 columns: trim positions, read
        length, and average EE per position
    """
    if not isinstance(FastqEntries, reader.FastqReader):
        raise TypeError("must be of FastqEntry type")
    if not isinstance(left, int):
        raise TypeError("left must be of type int")
    if not isinstance(right, int):
        raise TypeError("right must be of type int")
    if not isinstance(max_trim_perc, float):
        raise TypeError("max trim percentage must be of type float")

    # get average expected error per position
    avg_EE_df = FastqEntries.get_expected_error()
    avg_EE_list = list(avg_EE_df["AverageExpectedError"])

    for avgEE in avg_EE_list:
        if not isinstance(avgEE, float):
            raise TypeError("avg EEs must be of type float")

    # get length of reads
    read_len = len(avg_EE_list)

    # create lists with trim sites, read length, avg EE per position
    left_trim = []
    right_trunc = []
    read_len_list = []
    avgEE_position = []

    trim_bound = round(read_len * max_trim_perc)

    for current_left in range(left, trim_bound):
        for current_right in range(read_len - trim_bound, right):
            trim_readLen_avgEE = get_trim_length_avgEE(
                avg_EE_list, current_left, current_right
            )
            left_trim.append(trim_readLen_avgEE[0])
            right_trunc.append(trim_readLen_avgEE[1])
            read_len_list.append(trim_readLen_avgEE[2])
            avgEE_position.append(trim_readLen_avgEE[3])

    TrimInfo = pd.DataFrame(
        list(zip(left_trim, right_trunc, read_len_list, avgEE_position)),
        columns=["LeftIndex", "RightIndex", "ReadLength", "AvgEEPerPosition"],
    )

    TrimInfo = TrimInfo.reset_index(drop=True)

    return TrimInfo
