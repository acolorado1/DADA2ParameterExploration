import unittest
import warnings

import smartdada2.GetTrimParameters as GTP
from smartdada2.reader import reader


class MyTestCase(unittest.TestCase):
    def test_trim_ends_less_than_threshold(self):
        known_list = [
            15.0,
            20.0,
            22.0,
            30.0,
            30.0,
            40.0,
            30.0,
            30.0,
            25.0,
            20.0,
            15.0,
        ]
        result = GTP.trim_ends_less_than_threshold(known_list)
        self.assertEqual(result, (3, 8))

        # check raised errors
        string_list = ["30", "35", "36"]
        self.assertRaises(TypeError, GTP.trim_ends_less_than_threshold, string_list)

        impossible_values = [20.0, 50.0, 2.0]
        self.assertRaises(
            ValueError, GTP.trim_ends_less_than_threshold, impossible_values
        )

        # check return values for short return sizes
        trim_L_short = [
            15.0,
            18.0,
            18.0,
            30.0,
            30.0,
            40.0,
            30.0,
            30.0,
            25.0,
            20.0,
            15.0,
        ]
        trim_R_short = [
            15.0,
            20.0,
            22.0,
            30.0,
            30.0,
            40.0,
            30.0,
            30.0,
            18.0,
            18.0,
            15.0,
        ]

        with warnings.catch_warnings(record=True) as w:
            GTP.trim_ends_less_than_threshold(trim_L_short)
            GTP.trim_ends_less_than_threshold(trim_R_short)
            print(w[0].category)
            assert issubclass(w[0].category, UserWarning)
            assert issubclass(w[1].category, UserWarning)

        # check for no trimming
        perfect_list = [
            30.0,
            30.0,
            30.0,
            30.0,
            30.0,
            40.0,
            30.0,
            30.0,
            30.0,
            30.0,
            30.0,
        ]
        with warnings.catch_warnings(record=True) as w:
            result = GTP.trim_ends_less_than_threshold(perfect_list)
            self.assertEqual(result, (0, 11))
            assert issubclass(w[0].category, UserWarning)

    def test_get_trim_length_avgEE(self):
        # raise errors
        wrong_type = ["1", 2, 3]
        self.assertRaises(TypeError, GTP.get_trim_length_avgEE, wrong_type, 1, 2)

        # calculate average
        test_EE_list = [2.0, 3.0]
        manual_calc = (2 + 3) // 2
        results = GTP.get_trim_length_avgEE(test_EE_list, 0, 2)
        self.assertEqual(results[2], manual_calc)

    def test_read_size_by_avg_EE(self):
        # raise errors
        not_fastq = ["atcg"]
        self.assertRaises(TypeError, GTP.read_size_by_avg_EE, not_fastq, 0, 2)

        # test df output, this ensures reads are of equal length
        test_data_fp = reader.FastqReader("./test_data/SRR1591840_tunc.fastq")
        df = GTP.read_size_by_avg_EE(test_data_fp, 0, 150)
        self.assertEqual(len(df), 840)


if __name__ == "__main__":
    unittest.main()


# TODO
# test the no trim warning
