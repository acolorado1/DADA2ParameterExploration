import unittest

import numpy as np
import pandas as pd

import smartdada2.GetMaxEE as GME
import smartdada2.GetTrimParameters as GTP
from smartdada2.reader import reader


class MyTestCase(unittest.TestCase):
    def test_read_size_by_maxEE(self):

        fastqs = reader.FastqReader("./test_data/SRR1591840_tunc.fastq")

        # check input types
        self.assertRaises(
            TypeError,
            GME.read_size_by_maxEE,
            "./test_data/SRR1591840_tunc.fastq",
            1,
            150,
        )

        self.assertRaises(TypeError, GME.read_size_by_maxEE, fastqs, 1, "150")

        # make sure output is a pandas df
        df = GME.read_size_by_maxEE(fastqs, 1, 150)
        self.assertEqual(isinstance(df, pd.DataFrame), True)


if __name__ == "__main__":
    unittest.main()
