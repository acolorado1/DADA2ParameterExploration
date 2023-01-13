"""
Testing module tests all random, positive, and negative cases that FastqReader
may encounter
"""
import os
import unittest

import pandas
from pandas import DataFrame

from smartdada2.common.errors import FastqFormatError
from smartdada2.reader.reader import FastqEntry, FastqReader

# smartdada2 imports
from smartdada2.testing.help_test_funcs import toy_sequencer


class TestFastqEntry(unittest.TestCase):
    """FastqEntry is a datatype that is used within FastqReader. Each entry
    that is read within a FastqFile is stored in a FastqEntry data type.

    This datatype is important for all the methods to work within the
    FastqReader
    """

    # ------------------------------
    # Setup class methods
    # -- setups up files
    # -- removes files
    # ------------------------------
    @classmethod
    def set_up(cls) -> None:
        """Sets up files positive, negative, and random"""

    @classmethod
    def tear_down(cls) -> None:
        """removes all files"""
        os.remove(cls.upper_fastq)
        os.remove(cls.lower_fastq)

    # ------------------------------
    # Unit Testing
    # ------------------------------
    def test_entry(self):
        """Tests a positive case of entry"""

        # set entries
        raw_entries = [
            [
                "@test_fastq.test.1 1 length=50",
                "KDUHDNKSBARKSMKVNRBKWYRWDVBDYDRCUDYDVKANAAKKGBSBKG",
                "@test_fastq.test.1 1 length=50",
                "(#A!.&+,2@\"$&#&2('-231%-!+!$,!,;6%)',(>#;>20:\"$%#E",
            ],
            [
                "@test_fastq.test.2 1 length=50",
                "KWTBBNNDAVHWMYRRSBCRDSNGMDYYNSKTKNBVRDMUSHNRRSWDHK",
                "@test_fastq.test.2 1 length=50",
                "&\"@*%3\"&D((3'*!2#2D%%'#<&!2.,2&A)0,')+&<.'1+*--!+1",
            ],
        ]

        # create entry classes
        raw_reads_1 = raw_entries[0]
        test_entry_1 = FastqEntry(
            header=raw_reads_1[0],
            seq=raw_reads_1[1],
            scores=raw_reads_1[3],
            length=len(raw_reads_1[1]),
        )

        raw_reads_2 = raw_entries[1]
        test_entry_2 = FastqEntry(
            header=raw_reads_2[0],
            seq=raw_reads_2[1],
            scores=raw_reads_2[3],
            length=len(raw_reads_2[1]),
        )

        # setting expected variables
        expected_header_1 = "@test_fastq.test.1 1 length=50"
        expected_seq_1 = "KDUHDNKSBARKSMKVNRBKWYRWDVBDYDRCUDYDVKANAAKKGBSBKG"
        expected_score_1 = "(#A!.&+,2@\"$&#&2('-231%-!+!$,!,;6%)',(>#;>20:\"$%#E"
        expected_length_1 = 50
        expected_r_seq_1 = False

        expected_header_2 = "@test_fastq.test.2 1 length=50"
        expected_seq_2 = "KWTBBNNDAVHWMYRRSBCRDSNGMDYYNSKTKNBVRDMUSHNRRSWDHK"
        expected_score_2 = "&\"@*%3\"&D((3'*!2#2D%%'#<&!2.,2&A)0,')+&<.'1+*--!+1"
        expected_length_2 = 50
        expected_r_seq_2 = True

        # setting tests
        # -- testing first entry
        self.assertEqual(expected_header_1, test_entry_1.header)
        self.assertEqual(expected_seq_1, test_entry_1.seq)
        self.assertEqual(expected_score_1, test_entry_1.scores)
        self.assertEqual(expected_length_1, test_entry_1.length)
        # self.assertEqual(expected_r_seq_1, test_entry_1.rseq)

        # -- testing second entry
        self.assertEqual(expected_header_2, test_entry_2.header)
        self.assertEqual(expected_seq_2, test_entry_2.seq)
        self.assertEqual(expected_score_2, test_entry_2.scores)
        self.assertEqual(expected_length_2, test_entry_2.length)
        # self.assertEqual(expected_r_seq_2, test_entry_2.rseq)

    def test_entry_lowercases(self):
        """Tests if fastq files contains lower case sequences"""

        # set entries with lower case sequences
        raw_entries = [
            [
                "@test_fastq.test.1 1 length=50",
                "kduhdnksbarksmkvnrbkwyrwdvbdydrcudydvkanaakkgbsbkg",
                "@test_fastq.test.1 1 length=50",
                "(#A!.&+,2@\"$&#&2('-231%-!+!$,!,;6%)',(>#;>20:\"$%#E",
            ],
            [
                "@test_fastq.test.2 1 length=50",
                "kwtbbnndavhwmyrrsbcrdsngmdyynsktknbvrdmushnrrswdhk",
                "@test_fastq.test.2 1 length=50",
                "&\"@*%3\"&D((3'*!2#2D%%'#<&!2.,2&A)0,')+&<.'1+*--!+1",
            ],
        ]

        # create entry classes
        raw_reads_1 = raw_entries[0]
        test_entry_1 = FastqEntry(
            header=raw_reads_1[0],
            seq=raw_reads_1[1],
            scores=raw_reads_1[3],
            length=len(raw_reads_1[1]),
        )

        raw_reads_2 = raw_entries[1]
        test_entry_2 = FastqEntry(
            header=raw_reads_2[0],
            seq=raw_reads_2[1],
            scores=raw_reads_2[3],
            length=len(raw_reads_2[1]),
        )

        # setting expected variables
        expected_header_1 = "@test_fastq.test.1 1 length=50"
        expected_seq_1 = "KDUHDNKSBARKSMKVNRBKWYRWDVBDYDRCUDYDVKANAAKKGBSBKG"
        expected_score_1 = "(#A!.&+,2@\"$&#&2('-231%-!+!$,!,;6%)',(>#;>20:\"$%#E"
        expected_length_1 = 50
        expected_r_seq_1 = False

        expected_header_2 = "@test_fastq.test.2 1 length=50"
        expected_seq_2 = "KWTBBNNDAVHWMYRRSBCRDSNGMDYYNSKTKNBVRDMUSHNRRSWDHK"
        expected_score_2 = "&\"@*%3\"&D((3'*!2#2D%%'#<&!2.,2&A)0,')+&<.'1+*--!+1"
        expected_length_2 = 50
        expected_r_seq_2 = True

        # setting tests
        # -- testing first entry
        self.assertEqual(expected_header_1, test_entry_1.header)
        self.assertEqual(expected_seq_1, test_entry_1.seq)
        self.assertEqual(expected_score_1, test_entry_1.scores)
        self.assertEqual(expected_length_1, test_entry_1.length)
        # self.assertEqual(expected_r_seq_1, test_entry_1.rseq)

        # -- testing second entry
        self.assertEqual(expected_header_2, test_entry_2.header)
        self.assertEqual(expected_seq_2, test_entry_2.seq)
        self.assertEqual(expected_score_2, test_entry_2.scores)
        self.assertEqual(expected_length_2, test_entry_2.length)
        # self.assertEqual(expected_r_seq_2, test_entry_2.rseq)

    def test_invalid_header(self):
        """Negative case test to check if invalid headers are captured
        header id's require `@` in front of it"""

        # test entry
        raw_entries = [
            "test_fastq.test.1 1 length=50",
            "KDUHDNKSBARKSMKVNRBKWYRWDVBDYDRCUDYDVKANAAKKGBSBKG",
            "test_fastq.test.1 1 length=50",
            "(#A!.&+,2@\"$&#&2('-231%-!+!$,!,;6%)',(>#;>20:\"$%#E",
        ]

        header = raw_entries[0]
        seq = raw_entries[1]
        scores = raw_entries[2]
        length = len(seq)

        self.assertRaises(FastqFormatError, FastqEntry, header, seq, scores, length)


class TestFastqReader(unittest.TestCase):
    """FastqReader is the main engine that attempts to extract as much
    information from fastq sequences."""

    # ------------------------------
    # Unit testing
    # ------------------------------

    # -- Testing reader instantiation
    def test_reader_instance_1(self) -> None:
        """Positive case of a successful entry"""
        test_reader = FastqReader("./upper_seq.fastq")
        self.assertIsInstance(test_reader, FastqReader)

    def test_reader_instance_2(self) -> None:
        """Tests if can create fastq files with lower case sequences"""
        test_reader = FastqReader("./lower_seq.fastq")
        self.assertIsInstance(test_reader, FastqReader)

    def test_reader_invalid_seqs(self) -> None:
        """Tests if Fastq files contains unwanted types. For example, if digits
        or booleans or unknown characters are captured within the sequence"""
        test_reader = FastqReader("./invalid_seq.fastq")
        reads_iterator = test_reader.iter_reads()
        try:
            next(reads_iterator)
        except FastqFormatError as e:
            self.assertIsInstance(e, FastqFormatError)
        except Exception as e:
            self.fail(f"{e.__class__.__name__} captured, not FastqFormatError")

    def test_reader_invalid_scores(self) -> None:
        """Tests FastqReader for invalid scores"""
        test_reader = FastqReader("./invalid_scores.fastq")
        reads_iterator = test_reader.iter_reads()
        try:
            next(reads_iterator)
        except FastqFormatError as e:
            self.assertIsInstance(e, FastqFormatError)
        except Exception as e:
            self.fail(f"{e.__class__.__name__} captured, not FastqFormatError")

    def test_reader_file_not_exist(self) -> None:
        """Tests if exception is raised if the file is not found"""
        self.assertRaises(FileNotFoundError, FastqReader, "./does_not_exist.fastq")

    def test_reader_via_extension(self) -> None:
        """Positive test of using .fastq or .FASTQ as extensions"""
        # creating reader object with ".fastq" , ".FASTQ" ext
        test_reader_1 = FastqReader("./capital_ext_seq.FASTQ")

        # testing
        self.assertIsInstance(test_reader_1, FastqReader)

    # -- Testing quality scores
    def test_quality_score_type(self) -> None:
        """Checks if the correct type is returned is a DataFrame"""

        # setting up
        test_reader = FastqReader("./small.fastq")
        test_quality_score_test = test_reader.get_quality_scores()

        # test
        self.assertIsInstance(test_quality_score_test, DataFrame)

    def test_quality_score_values(self) -> None:
        """Tests if all values are equal. Using control randomized fastq files
        ./small.fastq"""
        expected_scores = [
            [7, 2, 32, 0, 13, 5, 10, 11, 17, 31, 1, 3, 5, 2, 5],
            [17, 7, 6, 12, 17, 18, 16, 4, 12, 0, 10, 0, 3, 11, 0],
            [11, 26, 21, 4, 8, 6, 11, 7, 29, 2, 26, 29, 17, 15, 25],
            [1, 3, 4, 2, 36, 5, 1, 31, 9, 4, 18, 1, 5, 35, 7],
            [7, 18, 6, 9, 0, 17, 2, 17, 35, 4, 4, 6, 2, 27, 5],
        ]

        # instantiating reader
        test_reader = FastqReader("./small.fastq")
        test_quality_score_test = test_reader.get_quality_scores().values.tolist()

        # testing
        self.assertEqual(expected_scores, test_quality_score_test)

    # -- Testing average quality score method
    def test_average_quality_score_type(self) -> None:
        """Checks if correct type is returned"""
        # setting up reader
        test_reader = FastqReader("./small.fastq")
        test_avg_scores = test_reader.get_average_score()

        # Testing
        self.assertIsInstance(test_avg_scores, DataFrame)

    def test_average_quality_scores_cols(self) -> None:
        """Checks if the correct column names are generated"""
        # expected answers
        expected_list = ["Position", "AverageQualityScore"]

        # setting up reader
        test_reader = FastqReader("./small.fastq")
        test_avg_scores = test_reader.get_average_score()
        test_colnames = test_avg_scores.columns.tolist()

        # test
        self.assertEqual(expected_list, test_colnames)

    def test_average_quality_score_values(self) -> None:

        # expected avg scores
        expected_avg_scores = [
            [0.0, 8.6],
            [1.0, 11.2],
            [2.0, 13.8],
            [3.0, 5.4],
            [4.0, 14.8],
            [5.0, 10.2],
            [6.0, 8.0],
            [7.0, 14.0],
            [8.0, 20.4],
            [9.0, 8.2],
            [10.0, 11.8],
            [11.0, 7.8],
            [12.0, 6.4],
            [13.0, 18.0],
            [14.0, 8.4],
        ]

        # setting up test
        test_reader = FastqReader("./small.fastq")
        test_avg_scores = test_reader.get_average_score().values.tolist()

        self.assertEqual(expected_avg_scores, test_avg_scores)

    # -- Testing Max expected error function
    def test_max_ee_type(self) -> None:
        """checks values produced"""
        # setting up reader
        test_reader = FastqReader("./small.fastq")
        test_avg_scores = test_reader.get_max_seq_ee()

        # Testing
        self.assertIsInstance(test_avg_scores, DataFrame)

    def test_max_ee_colnames(self) -> None:
        """Checks if the correct column names are generated"""
        # expected answers
        expected_list = ["length", "max_ee"]

        # setting up reader
        test_reader = FastqReader("./small.fastq")
        test_max_scores = test_reader.get_max_seq_ee()
        test_colnames = test_max_scores.columns.tolist()

        # test
        self.assertEqual(expected_list, test_colnames)

    def test_max_ee_values(self) -> None:
        """checks values produced"""

        # expected
        expected_values = [
            [15.0, 4.96],
            [15.0, 4.74],
            [15.0, 1.87],
            [15.0, 5.29],
            [15.0, 4.26],
        ]

        # setup test
        test_reader = FastqReader("./small.fastq")
        test_max_ee = test_reader.get_max_seq_ee().values.tolist()

        # self
        self.assertEqual(expected_values, test_max_ee)

    # -- Testing Max average expected error
    def test_avg_ee_scores_type(self) -> None:
        """checking types"""
        # setting up reader
        test_reader = FastqReader("./small.fastq")
        test_avg_scores = test_reader.get_avg_seq_ee()

        # Testing
        self.assertIsInstance(test_avg_scores, DataFrame)

    def test_avg_ee_scores(self) -> None:
        """check values of average scores"""

        # expected values
        expected_values = [
            [15.0, 4.3],
            [15.0, 6.1],
            [15.0, 7.9],
            [15.0, 4.2],
            [15.0, 9.4],
        ]

        # testing
        test_reader = FastqReader("./small.fastq")
        test_max_ee = test_reader.get_avg_seq_ee().values.tolist()

        self.assertEqual(expected_values, test_max_ee)

    def test_avg_ee_type(self) -> None:
        """checks values produced"""
        # setting up reader
        test_reader = FastqReader("./small.fastq")
        test_avg_scores = test_reader.get_average_score()

        # Testing
        self.assertIsInstance(test_avg_scores, DataFrame)

    # -- testing unpacking of sequence entries
    def test_unpack_entries(self):
        """Takes in FastqReader object and extract. Checks if the generated
        list is an instance of FastqEntry
        """

        # init reader
        test_reader = FastqReader("./small.fastq")

        # unpack entries
        unpack_reads = test_reader.unpack_entries()

        # check
        for entry in unpack_reads:
            self.assertIsInstance(entry, FastqEntry)

    def test_unpack_full_entries(self):
        """Unpacks reads and returns a python object. This is considered as a
        full unpacking due to the destructuring of the FastqEntry and only
        returning python objects."""

        # expected output
        expected_arr = [
            [
                "@test_fastq.test.1 1 length=15",
                "KDUHDNKSBARKSMK",
                '(#A!.&+,2@"$&#&',
                15,
                "forward",
            ],
            [
                "@test_fastq.test.2 1 length=15",
                "VNRBKWYRWDVBDYD",
                "2('-231%-!+!$,!",
                15,
                "forward",
            ],
            [
                "@test_fastq.test.1 2 length=15",
                "RCUDYDVKANAAKKG",
                ",;6%)',(>#;>20:",
                15,
                "forward",
            ],
            [
                "@test_fastq.test.2 2 length=15",
                "BSBKGKWTBBNNDAV",
                '"$%#E&"@*%3"&D(',
                15,
                "forward",
            ],
            [
                "@test_fastq.test.1 3 length=15",
                "HWMYRRSBCRDSNGM",
                "(3'*!2#2D%%'#<&",
                15,
                "forward",
            ],
        ]
        # init reader
        test_reader = FastqReader("./small.fastq")

        # unpack entries
        fully_unpacked_reads = test_reader.unpack_entries(full=True)

        # testing if all the values inside array are the same
        self.assertEqual(expected_arr, fully_unpacked_reads)

    def test_sequence_df(self) -> None:
        """Builds a pandas dataframe containing a sequence"""

        # building expected output
        expected_seqs = []
        with open("./small.fastq", "r") as f:
            contents = f.readlines()
            for i in range(1, len(contents), 4):
                seq = list(contents[i].rstrip())
                expected_seqs.append(seq)

        # test reader
        test_reader = FastqReader("./small.fastq")

        # build sequence df
        test_seq_df = test_reader.sequence_df().values.tolist()

        self.assertEqual(test_seq_df, expected_seqs)

    # ------------------------------
    # Setup class methods
    # -- setups up files
    # -- removes files
    # ------------------------------
    @classmethod
    def setUp(cls) -> None:
        """Sets up files positive, negative, and random"""

        # setting vars
        cls.small_fastq = "small.fastq"
        cls.upper_fastq = "upper_seq.fastq"
        cls.lower_fastq = "lower_seq.fastq"
        cls.invalid_seq = "invalid_seq.fastq"
        cls.invalid_scores = "invalid_scores.fastq"
        cls.capital_ext = "capital_ext_seq.FASTQ"
        cls.invalid_ext = "invalid_ext_seq.fasta"
        cls.empty_file = "empty.fastq"

        # generating small fastq file
        with open(cls.small_fastq, "w") as f:
            reads = toy_sequencer(15, 5, rev_seq=True, seed=42)
            for read in reads:
                for read_data in read:
                    f.write(f"{read_data}\n")

        # generating regular fastq fille
        with open(cls.upper_fastq, "w") as f:
            reads = toy_sequencer(200, 30, rev_seq=True, seed=42)
            for read in reads:
                for read_data in read:
                    f.write(f"{read_data}\n")

        # generating sequence with lower case sequence
        with open(cls.lower_fastq, "w") as f:
            reads = toy_sequencer(200, 30, rev_seq=True, seed=42)
            for read in reads:
                read[1] = read[1].lower()
                for read_data in read:
                    f.write(f"{read_data}\n")

        # generating sequence with invalid seq values
        with open(cls.invalid_seq, "w") as f:
            read_entries = [
                "@test_fastq.test.3 1 length=50",
                "KWTBBNN1AVHWMYR23ZCRDSNGMDYYNSKTKNBV@DZUSHNRZSWDHK",
                "@test_fastq.test.2 1 length=50",
                "&\"@*%3\"&D((3'*!2#2D%%'#<&!2.,2&A)0,')+&<.'1+*--!+1",
            ]
            for read in read_entries:
                f.write(f"{read}\n")

        # generating invalid scores
        with open(cls.invalid_scores, "w") as f:
            read_entries = [
                "@test_fastq.test.3 1 length=50",
                "KWTBBNN1AVHWMYR23ZCRDSNGMDYYNSKTKNBV@DZUSHNRZSWDHK",
                "@test_fastq.test.2 1 length=50",
                "&\"@*%3\"&D((Z'*!2#2D%%'#<&Z2.,2ZA)0,')+&<.'1+*--!+1",
            ]
            for read in read_entries:
                f.write(f"{read}\n")

        # creating fastq file with capital extension
        with open(cls.capital_ext, "w") as f:
            reads = toy_sequencer(200, 30, rev_seq=True, seed=42)
            for read in reads:
                for read_data in read:
                    f.write(f"{read_data}\n")

        # creating blank file with invalid ext
        with open(cls.invalid_ext, "w") as f:
            f.write("invalid ext")

        # creating empty fastqfile
        with open(cls.empty_file, "w") as f:
            f.write("")

    @classmethod
    def tearDown(cls) -> None:
        """removes all files"""
        os.remove(cls.small_fastq)
        os.remove(cls.upper_fastq)
        os.remove(cls.lower_fastq)
        os.remove(cls.invalid_seq)
        os.remove(cls.invalid_scores)
        os.remove(cls.invalid_ext)
        os.remove(cls.capital_ext)
        os.remove(cls.empty_file)


if __name__ == "__main__":

    unittest.main()
