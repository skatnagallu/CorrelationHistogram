import unittest
import correlationHistogramAnalysis.dissociation_tracks as diss

class TestCorrHist(unittest.TestCase):
    def test_parse_formula(self):
        e,_= diss.parse_formula('Fe2O3')
        self.assertTrue(len(e)==1)
        self.assertEqual(list(e.keys())[0],'Fe')

if __name__ == "__main__":
    unittest.main()