
# CORA CI

The CORA CI runs the test suites in a docker container (see `HowTo-CORA-docker.md`).
This file describes the structure of the CI pipeline.

The CORA CI has 4 stages:
- A short stage
- An extended stage
- A compatibility stage
- An examples stage

You can find more information about each stage below.

---

### Merging into Main Branches

If you want to merge your branch into one of the main branches of CORA (`devgeneral`, `public-bugfix`, (`PUBLIC`)), 
please create a merge request which triggers all necessary tests automatically.
You can also trigger the merge pipeline manually by appending your commit message with:

    git commit -m "<commit-message> --ci-run-merge"

which runs everything up to the example stage.

---

## CI Test Stages

### Short Test Stage

This test stage runs on every branch after each commit pushed to Bitbucket.
It runs the short test suite to test the basic functionality.

    runTestSuite('short')

This test suite runs all unit tests with the `test_` prefix.
Each test should finish within a few seconds.
Longer test should be considered to move to the `long` test suite
and only add a small test to the `short` test suite.


### Extended Test Stage

This stage runs extensive tests on your commits.
This test stage always runs on the git branches `PUBLIC`, `devgeneral`, and `public-bugfix`.

It runs the following test suites:

    runTestSuite('long')
    runTestSuite('nn')
    runTestSuite('flaky')

where the `flaky` tests can sometimes fail for various reason (e.g., bad rng), but never more than two should fail.

#### Neural Network Tests

For branches starting with the prefix `nn-*`, 
the neural network test suite is run automatically.

    runTestSuite('nn')

Note that some of these tests require the necessary toolboxes to be installed (see manual).


### Compatibility Test Stage

CORA tries to be compatible with all recent Matlab versions.
However, only the latest Matlab version when the last major version of CORA was released, is thoroughly tested.
Usually, thats Matlab R202xb as we aim to release major versions of CORA in fall.

It runs the following test suites for various matlab version:

    runTestSuite('short')

This stage can be triggered by ending the commit message with `--ci-run-compatibility`.


### Examples Test Stage

This stage tests all examples of CORA and can be triggered 
by ending the commit message with `--ci-run-examples`.

It runs the following test suites:

    runTestSuite('examples')
    runTestSuite('benchmarks') % long examples
    runTestSuite('header')     % examples in header

