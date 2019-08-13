import unittest

from isodesmic import IncrementingDictionary


class TestConstraintMap(unittest.TestCase):
    """
    Test that ConstraintMap objects function properly
    """

    def test_constraint_map(self):
        """
        Test that ConstraintMap objects behave like dictionaries but add new elements to the dictionary whenever a
        value is requested for a key currently not in the dictionary
        """
        constraint_map = IncrementingDictionary()
        self.assertEqual(constraint_map['a'], 0)
        self.assertEqual(constraint_map['b'], 1)
        self.assertEqual(len(constraint_map), 2)

        # Try to reference the key 'a' again
        self.assertEqual(constraint_map['a'], 0)

        self.assertEqual(constraint_map['c'], 2)


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))