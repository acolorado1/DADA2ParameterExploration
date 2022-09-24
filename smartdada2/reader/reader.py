import sys
import warnings
from pathlib import Path
from dataclasses import dataclass
from typing import Iterable
from typing import Optional
from typing import Union

import pandas as pd
import numpy as np

from smartdada2.seq_utils.search import binary_search


@dataclass(slots=True)
class FastqEntry:
    """Contains the contents of a single read as a FastqEntry"""

    header: str
    seq: str
    scores: str
    length: int

    # set to none unless user declares there is reversed sequences
    rseq: Union[None, int] = None


class FastqReader:
    """Memory efficient FastqReader"""

    def __init__(
        self,
        fpath: str,
        reverse_seq: Optional[bool] = False,
        technology: Optional[str] = "illumina",
    ):

        # FastqReader accessible parameters
        self.fpath: Path = Path(fpath)
        self.reversed: bool = reverse_seq
        self.technology: str = technology

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
        for entry in self.__loader():
            phred_scores = list(entry.scores)
            scores = [ord(phred_score) - 33 for phred_score in phred_scores]
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
        average_score = scores_df.mean()
        return average_score

    def iter_reads(self) -> Iterable[FastqEntry]:
        """Returns a python generator containing FastqEntries

        Returns
        -------
        Sequence[FastqEntry]
            Generator object containing FastqEntries
        """
        return self.__loader()

    def to_list(self) -> list[FastqEntry]:
        """Saves all lists into memory. Warning, large file sizes will use more
        memory but increase iteration performance.

        Returns
        -------
        list[FastqEntry]
            List of FastqEntries
        """
        return [read for read in self.__loader()]

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
            self.__n_entries = sum(1 for _ in self.__loader())
            self.__counted = True
            return self.__n_entries
        else:
            return self.__n_entries

    # TODO: Should this be a class method?
    # Should this write out temporary file creating a sub_sample read
    def sample(
        self,
        frac: Optional[float] = 0.3,
        sample_type: Optional[str] = "random",
        prob: Optional[float] = 0.5,
        seed: Optional[int] = None,
    ) -> Iterable[FastqEntry]:
        """Sub samples

        Parameters
        ----------
        frac : Optional[float], optional
            subsample size, by default 0.3
        sample_type : Optional[str], optional
            Subsample algorithms: random or reservoir, by default "random".
        prob : Optional[float], optional
            If sample_type is  reservoir, then a probability rate is required to select entries, by default 0.5
        seed : Optional[int], optional
            Value that conducts controlled randomness, by default None

        Yields
        ------
        Iterator[Iterable[FastqEntry]]
            Generator object that contains FastqEntry 

        Raises
        ------
        ValueError
            Raised when either frac or prob is larger than 1.0
        """
        # type checks
        if not isinstance(frac, float):
            frac = float(frac)
        
        
        # selecting sampling type
        match sample_type:
            case "random":

                # checking if fraction subsample is larger than the whole dataset
                if frac > 1.0:
                    raise ValueError("frac cannot be larger than 1.0")

                # checking if reads have been counted
                # this will update the `self.__n_entries` variable
                if self.__counted is False:
                    self.total_reads()

                # initializing random seed if Seed is not None
                # -- setting the seed will have a controlled randomness
                # -- good for reproducibility (share seed value)
                if seed is not None:
                    np.random.seed(seed)

                # calculate subsample_suze (how many entries we want)
                subsample_size = int(np.round(self.__n_entries * frac))

                # selecting random indices without replacement
                # range selection is (0, total_n_reads)
                # -- done without replacement (no repeated indices selected)
                random_indices = np.sort(
                    np.random.choice(
                        np.arange(0, self.__n_entries), 
                        size=subsample_size, 
                        replace=False
                    )
                )
                
                # to stop further computation once all indices are found
                stop_idx = np.max(random_indices)

                # loading the selected reads
                for idx, read_entry in enumerate(self.__loader()):

                    # stop if next read idx is largest than stop idx
                    if idx > stop_idx:
                        break

                    # if index exists within randomly selected indices, yield read
                    try:
                        binary_search(idx, random_indices, sorted=True)
                        yield read_entry
                    
                    # if it does not exists, go to the next
                    except ValueError:
                        continue

            case "reservoir":
                if frac > 1.0 or prob > 1.0:
                    raise ValueError("frac or prob cannot be greater than 1.0")
                raise NotImplementedError("Reservoir sampling has not been developed yet")
            case _:
                raise ValueError("Invalid sampling choice")

    # ----------------------------------------
    # private functions: users do not interact with this
    # ----------------------------------------
    def __loader(self) -> Iterable[FastqEntry]:
        """Creates a generator object that contains sequence read data as
        FastqEntry object.

        While iterating, a count is also being conducted

        Returns
        -------
        Sequence[FastqEntry]
            Generator object with FastqEntries
        """

        # iterate all row contents in fastq file
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

                    # convert into FastqEntry
                    entry_count += 1
                    fastq_entry = FastqEntry(
                        header=contents_chunk[0],
                        seq=contents_chunk[1],
                        scores=contents_chunk[3],
                        length=len(contents_chunk[1]),
                    )

                    # clear list
                    contents_chunk = []

                    # yield entry
                    yield fastq_entry

        self.__n_entries = entry_count
        self.__counted = True
