import itertools
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, List, Optional, Union

import numpy as np
import pandas as pd

from smartdada2.common.errors import FastqFormatError

# constants
DNA = list("ATCGU")
AMB_DNA = list("NRYKMSWBDHV")
ALL_DNA = DNA + AMB_DNA
ASCII_SCORES = [chr(i) for i in range(33, 70 + 1)]
NUMERICAL_SCORES = list(range(33, 70 + 1))


# @dataclass(slots=True)
@dataclass
class FastqEntry:
    """Contains the contents of a single read as a FastqEntry"""

    header: str
    seq: str
    scores: str
    length: int

    rseq: Union[None, bool] = None

    def __post_init__(self):

        # checking header (1 == forward, 2 == rev)
        self.check_seq_dir()

    def check_seq_dir(self):
        """Checks if the entry is a forward or reverse sequence"""

        # parse the header
        header_id = self.header.split()[0]
        if not header_id.startswith("@"):
            raise FastqFormatError("Unable to find header id")

        # find the direction of the sequence
        direction = header_id.split(".")[-1]

        # convert string into integer type
        # if not direction.isdigit():
        #     raise FastqFormatError("Unable to find the sequence direction")
        # else:
        #     direction = int(direction)

        # setting up sequence direction into FastqEntry
        # if direction == 1:
        #     self.rseq = False
        # elif direction == 2:
        #     self.rseq = True
        # else:
        #     raise FastqFormatError("Unable to find the sequence direction")

        # making sequences to be capital letters
        self.seq = self.seq.upper()


class FastqReader:
    """
    FastqReader is an efficient Fastq file reader that optimizes for memory
    usage and speed.

    Every entry within the Fastq file is stored within a FastqEntry object,
    which contains header, sequence, sequence length, score, and direction of
    each read.

    FastqReader contains methods that are meant for gathering information
    within fastq files.
    """

    def __init__(
        self,
        fpath: str,
        reverse_seq: Optional[bool] = False,
        technology: Optional[str] = "illumina",
        scale: Optional[str | None] = None,
    ):

        # FastqReader accessible parameters
        self.fpath: Path = Path(fpath)
        if not self.fpath.is_file():
            raise FileNotFoundError(f"{self.fpath} does not exist")

        self.reversed: bool = reverse_seq
        self.technology: str = technology

        self.scale = None
        if scale is not None:
            print(f"scaling reads to {scale}")
            self.scale = scale

        # states
        self.__n_entries: Optional[int] = 0
        self.__counted: bool = False

    def get_quality_scores(self) -> pd.DataFrame:
        """Returns scores in a per sequence bases

        Returns
        -------
        pd.DataFrame
            DataFrame containing quality scores per sequence
        """

        # converting to np.array
        all_scores = []
        for entry in self.iter_reads():
            phred_scores = list(entry.scores)
            scores = np.array(phred_scores).view(np.int32) - 33
            all_scores.append(scores)

        all_scores = pd.DataFrame(data=np.array(all_scores))
        return all_scores

    def get_average_score(self) -> pd.Series:
        """Returns average score of all sequences. Returns a a pd.Series object

        Returns
        -------
        pd.Series
            average score of all sequences

        """

        # get all scores df and take the average per column basis
        scores_df = self.get_quality_scores()
        average_score = scores_df.mean().reset_index()
        average_score.columns = ["Position", "AverageQualityScore"]
        return average_score

    def get_expected_error(self) -> pd.DataFrame:
        """Calculates average expected error per nucleotide

        Returns
        -------
        pd.DataFrame
            position and expected error in a dataframe
        """
        # loading average quality scores
        scores_df = self.get_average_score()

        # convert quality score values into expected errors values
        scores_df["AverageExpectedError"] = scores_df["AverageQualityScore"].apply(
            lambda score: 10 ** (-score / 10)
        )
        return scores_df.drop(columns="AverageQualityScore")

    def get_seq_ee_errors(self) -> pd.DataFrame:
        """Calculates expected error in all sequences

        Returns
        -------
        pd.DataFrame
            A dataframe that each row contains the expected error per
            nucleotide in each sequence. This is different from
            `get_expected_error()` as it does not take the average of all
            scores in within each sequence.
        """
        # quality scores
        qual_scores_df = self.get_quality_scores()
        return qual_scores_df.apply(lambda row: 10 ** (-row / 10), axis=1)

    def get_max_seq_ee(self) -> pd.Series:
        """Creates a pandas dataframe object that contains the sum of expected
        error of all sequences. The index will represent the sequence. The two
        columns "direction" and "total_error"

        Returns
        -------
        pd.DataFrame
            A two column dataframe that contains "direction" and "total_error"
        """

        # creates a base dataframe for calculating expected errors
        expected_error_df = self.__base_max_ee_df()

        # get raw
        raw_scores = self.get_seq_ee_errors()
        expected_error_df["max_ee"] = raw_scores.apply(
            lambda row: np.round(np.sum(row), 2), axis=1
        )

        # rearranging columns
        expected_error_df = expected_error_df[
            # ["length", "direction", "max_ee"]
            ["length", "max_ee"]
        ]

        return expected_error_df

    def get_avg_seq_ee(self) -> pd.DataFrame:
        """Get average expected error per sequence. Generates a pandas
        dataframe that contains sequence length, direction, and average error
        score

        Returns
        -------
        pd.DataFrame
            average expected error per sequence dataframe
        """
        # create empty df
        expected_error_df = self.__base_max_ee_df()

        # get raw ee scores
        raw_scores = self.get_average_score()
        expected_error_df["avg_ee"] = raw_scores.apply(lambda row: np.mean(row), axis=1)

        # return expected_error_df[["length", "direction", "avg_ee"]]
        return expected_error_df[["length", "avg_ee"]]

    def sequence_df(self) -> pd.DataFrame:
        """
        Returns a Dataframe structure of sequence reads. Row represents a
        sequence and the columns represents the individual nucleotides
        """
        return pd.DataFrame(data=(list(entry.seq) for entry in self.iter_reads()))

    def ambiguous_nucleotide_counts(self) -> pd.DataFrame:
        """Counts all ambiguous nucleotides in all reads.

        Returns
        -------
        pd.DataFrame
            DataFrame that contains the nucleotide position and the number
            of ambiguous nucleotides
        """
        seq_df = self.sequence_df()
        ambi_nuc_counts = (
            seq_df.apply(search_ambiguous_nucleotide).to_frame().reset_index()
        )
        ambi_nuc_counts.columns = ["Position", "AmbiguousCounts"]
        return ambi_nuc_counts

    def total_reads(self) -> int:
        """Returns total number of reads. Changes the value of self.n_entries
        to prevent any recalculation

        Returns
        -------
        int
            total number of reads
        """
        # checks if the n_entries have been counted before
        # -- if not, count and return value
        # -- if yes, return stored count
        if self.__n_entries == 0 and self.__counted is False:
            warnings.warn("No reads have been counted. Counting reads now...")
            self.__n_entries = sum(1 for _ in self.iter_reads())
            self.__counted = True
        return self.__n_entries

    def get_forward_reads(self) -> Iterable[FastqEntry]:
        """Obtains all forward reads within fastq files

        Parameters
        ----------
        tolist : bool, optional
            Saves the forward reads in memory (list), by default False

        Returns
        -------
        Iterable[FastqEntry]
            returns a generator if tolist is False. if tolist is True, then the
            reads will be saved into memory and will return a python list
            object.
        """

        # checking if the user wants to load reads in memory
        for read in self.__loader():
            if read.rseq is False:
                yield read

    def get_reverse_reads(self, tolist=False) -> Iterable[FastqEntry]:
        """Obtains all reverse reads within fastq files

        Parameters
        ----------
        tolist : bool, optional
            Saves the reverse reads in memory (list), by default False

        Returns
        -------
        Iterable[FastqEntry]

        """
        for read in self.__loader():
            if read.rseq is True:
                yield read

    def reservoir_sampling(
        self, n_samples: int, seed: Optional[bool] = None
    ) -> list[FastqEntry]:
        """Generates a uniform random subset sample of FastqEntry.

        Parameters
        ----------
        n_samples : int
            size of sub sample
        seed : Optional[int], optional
            Value that retains the state randomness, by default None

        Returns
        -------
        list[FastqEntry]
            subset of FastqEntry objects
        """

        # type checking
        if seed is not None and not isinstance(seed, int):
            raise ValueError(
                f"seed value must be an integer type, provided: {type(seed)}"
            )

        # set a results array and have elements inside as place holders
        # -- seen is equals to n_samples because those elements make up the
        # -- results array
        seen = n_samples
        subset_reads = []

        # populating results array with FastqEntry objects as place holders

        # -- checking if number of entries have counted before
        # -- -- if yes, set slice of FastqEntries as placeholders
        if self.__counted is True:
            if self.__n_entries < n_samples:
                raise ValueError(
                    "requested sample size is larger than number of entries"
                )
            else:
                subset_reads = list(itertools.islice(self.iter_reads(), n_samples))

        else:
            try:
                entries = self.iter_reads()
                subset_reads.extend(next(entries) for _ in range(n_samples))
            except StopIteration as e:
                raise StopIteration(
                    "number of requested samples exceeded number of entries"
                ) from e

        # next is to update the results array
        # -- checking if seed has been added
        if seed is not None:
            np.random.seed(seed)

        # -- now replacing placeholders by randomly selecting values as ind
        for entry in self.iter_reads():

            # select random int
            seen += 1
            rand_int = np.random.randint(0, seen - 1)

            # replace placeholders via index reassignment
            if rand_int < n_samples:
                subset_reads[rand_int] = entry

        return subset_reads

    def unpack_entries(
        self, full: Optional[bool] = False
    ) -> List[Union[FastqEntry, List[Any]]]:
        """Takes the iterator object produced from FastqReader and saves them
        into memory as a python list.

        Parameters
        ----------
        full : bool, optional
            indicates if the user wants a full unpacking of the dataset. If
            False, the returning object will be python list of FastqEntries. If
            set to True, the return object is a list of list containing
            entry information. by default False

        Returns
        -------
        List[FastqEntry]
            list of FastqEntry object will be returned if full is set to False.
            If full is set to True, then the FastqEntry object will be fully
            unpacked and returns unstructured entry.

            The unstructured dataset is a python list that contains
            [header id, sequence, scores, sequence length and direction]
        """

        # checking if user want list of FastqEntries or Raw Python data types
        # -- returning list of FastqEntry objects
        if full is False:
            return list(self.iter_reads())

        elif full is True:

            all_unpacked_entries = []
            for entry in self.iter_reads():

                # unpack all reads
                # -- unpacking
                unpacked_entries = [
                    entry.header,
                    entry.seq,
                    entry.scores,
                    int(entry.length),
                ]

                # -- converting sequence direction to string types
                if entry.rseq is True:
                    unpacked_entries.append("reverse")
                else:
                    unpacked_entries.append("forward")

                # storing
                all_unpacked_entries.append(unpacked_entries)

            return all_unpacked_entries

    def iter_reads(self) -> Iterable[FastqEntry]:
        """Returns a python generator containing FastqEntries

        Returns
        -------
        Sequence[FastqEntry]
            Generator object containing FastqEntries
        """
        for entry in self.__loader():

            # scaling information
            if self.scale is not None:
                entry.seq = entry.seq[: self.scale]
                entry.scores = entry.scores[: self.scale]
                entry.length = self.scale

            yield entry

    def slice_reads(
        self, range_idx: Union[tuple[int, int], list[int, int]], to_list=False
    ) -> Union[list[FastqEntry], Iterable[FastqEntry]]:
        """Allows index slicing for reads. If to_list is True, the method will
        return a list else it will return a generator

        Parameters
        ----------
        range_idx : tuple
            start and ending position (exclusive), to which reads to select.
        to_list : bool, optional
            convert to list object, If False, returns a generator.
            by default False

        Returns
        -------
        Union[list[FastqEntry], Iterable[FastqEntry]]
        return a list or iterator of FastqEntry objects

        Yields
        ------
        Iterator[Union[list[FastqEntry], Iterable[FastqEntry]]]
            if to_list is false, returns a generator of reads, else it will
            return a list of reads

        Raises
        ------
        TypeError
            Raised when a tuple or list is not provided
        ValueError
            Raised if range_idx does not contains 2 values
        TypeError
            Raised if values in range_idx are not integers
        ValueError
            raised if starting position value is larger than the ending
            position
        """
        if not isinstance(range_idx, tuple) and not isinstance(range_idx, list):
            raise TypeError(
                "Please provide a tuple or lists with starting and ending idx"
            )
        elif len(range_idx) != 2:
            raise ValueError("'range_idx' only takes two value (start, end)")
        elif not all(isinstance(value, int) for value in range_idx):
            raise TypeError("Values must be integers. Not floats, strings or booleans")
        elif range_idx[0] > range_idx[1]:
            raise ValueError("starting position cannot be larger than ending position")

        if to_list is True:
            return list(self.__slice(range_idx))
        else:
            return self.__slice(range_idx)

    # ----------------------------------------
    # private functions: users do not interact with this
    # ----------------------------------------
    def __loader(self) -> Iterable[FastqEntry]:
        """Creates a generator object that streams sequence reads as
        FastqEntry objects.

        While iterating, a count is also being conducted.

        Returns
        -------
        Sequence[FastqEntry]
            Generator object with FastqEntries
        """

        if self.fpath.suffix.lower() != ".fastq":
            raise ValueError("FastqReader only takes files '.fastq' or '.FASTQ' files")
        elif self.fpath.stat().st_size == 0:
            raise FastqFormatError("Fastq file contains no contents")

        # iterate all row contents in fastq file and collect entries
        entry_count = 0
        with open(self.fpath, "r") as fastq_file:

            contents_chunk = []
            for row_entry in fastq_file:

                # cleaning entries
                if row_entry == "":
                    continue

                content = row_entry.rstrip("\n")
                contents_chunk.append(content)

                # checking if there are 4 elements in the list
                # -- 4 lines = 1 entry
                if len(contents_chunk) == 4:

                    # check for valid sequences
                    seq_check = set(contents_chunk[1]) - set(ALL_DNA)
                    if len(seq_check) > 0:
                        raise FastqFormatError(
                            "File contains invalid sequence characters"
                        )

                    # check for valid ascii scores
                    # score_check = set(contents_chunk[3]) - set(ASCII_SCORES)
                    # if len(score_check) > 0:
                    #     raise FastqFormatError(
                    #         "File contains invalid score characters"
                    #     )

                    # convert into FastqEntry
                    entry_count += 1

                    # scaling sequence length
                    seq_length = len(contents_chunk[1][: self.scale])
                    if self.scale is not None:
                        seq_length = len(contents_chunk[1])

                    fastq_entry = FastqEntry(
                        header=contents_chunk[0],
                        seq=contents_chunk[1][: self.scale],
                        scores=contents_chunk[3][: self.scale],
                        length=seq_length,
                    )

                    # clear list
                    contents_chunk = []

                    # yield entry
                    yield fastq_entry

        self.__n_entries = entry_count
        self.__counted = True

    def __slice(self, range_idx: tuple[int, int]) -> Iterable[FastqEntry]:
        """ "Creates a generator of py FastqEntry object by a given range

        Parameters
        ----------
        idx_range : tuple
            start, end range of selecting FastqEntries

        Returns
        -------
        Iterable[FastqEntry]
            Generator with FastqEntry objects

        Yields
        ------
        Iterator[Iterable[FastqEntry]]
            FastqEntry objects
        """

        # unpacking range
        start, end = range_idx

        # iterating reads
        for idx, read in enumerate(self.iter_reads()):

            # break the iteration if the current idx == ending idx
            if idx == end:
                break

            # only selecting reads between start and end positions
            elif idx >= start and idx < end:
                yield read

    def __base_max_ee_df(self) -> pd.DataFrame:
        """Creates a base dataframe for expected errors calculations"""
        # create empty df
        base_ee_df = pd.DataFrame()

        # extract read direction
        # _dir = []
        # for read in self.iter_reads():
        #     if read.rseq is True:
        #         _dir.append("forward")
        #     else:
        #         _dir.append("reverse")

        # add sequence direction
        # base_ee_df["direction"] = _dir

        # add sequence length informatoin
        base_ee_df["length"] = [entry.length for entry in self.iter_reads()]

        return base_ee_df


# NOTE: Should this be in the FastqReader class or a separate function?
def search_ambiguous_nucleotide(nucleotides: pd.Series) -> int:
    """
    Parameters
    ----------
    nucleotides : pd.Series
        Series of nucleotides characters.

    Returns
    -------
    int
        total number of all ambiguous nucleotides found.
    """

    # counting all nucleotides
    count_series = nucleotides.value_counts()

    # searching for  nucleotides
    found_ambiguous_nucleotide = [
        ambi_nuc for ambi_nuc in AMB_DNA if ambi_nuc in count_series.index.tolist()
    ]

    return count_series[found_ambiguous_nucleotide].sum()
