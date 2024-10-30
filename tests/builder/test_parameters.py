from dtcc_core import builder

import unittest


class TestParameters(unittest.TestCase):
    def test_default_parameters(self):
        p = builder.parameters.default()
        self.assertIsInstance(p, dict)
        self.assertEqual(p["model_name"], "DTCC")

if __name__ == "__main__":
    unittest.main()