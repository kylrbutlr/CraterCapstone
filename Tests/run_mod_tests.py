import unittest

"""
Import test modules here
"""


"""
Add module to mod_tests list here
"""
mod_tests = [
    
]

"""
Code to run all module tests
"""
VERBOSITY = 2 # Test output detail level

def main():
    for mod_test in mod_tests:
        suite = unittest.TestLoader().loadTestsFromModule(mod_test)
        unittest.TextTestRunner(verbosity=VERBOSITY).run(suite)

if __name__ == "__main__":
    main()
