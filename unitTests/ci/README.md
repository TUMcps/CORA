
# CORA CI

The CORA CI runs the test suites in a docker container (see `HowTo-CORA-docker.md`).
This file describes the structure of the CI pipeline.

The CORA CI has 5 stages:
- A short stage
- An extended stage
- An examples stage
- A benchmarks stage
- A status update stage

You can find more information about each stage below.

---

### Merging into `devgeneral`

If you want to merge your branch into `devgeneral`, 
please trigger the merge pipeline by appending your commit message with:

    git commit -m "<commit-message> --ci-run-merge"

which runs everything up to the example stage.

### Preparing a release

Before each (major) release, we run all tests by appending the commit message with:

    git commit -m "<commit-message> --ci-run-release"

which runs all test stages.

Note that the Bitbucket repository is mirrored to GitLab to run the CI locally on our servers.
The results are then pushed back to Bitbucket and are visible next to the respective commit.
A click on the icon tells you more about the result of the CI.
https://bitbucket.org/MatthiasAlthoff/cora/commits/

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
This stage can be triggered by appending your commit message with `--ci-run-ext`.

It runs the following test suites:

    runTestSuite('long')
    runTestSuite('nn')
    runTestSuite('flaky')

where the `flaky` tests can sometimes fail for various reason, but never more than two should fail.

#### Neural Network Tests

For branches starting with the prefix `nn-*`, 
the neural network test suite is run automatically.

    runTestSuite('nn')

Note that some of these tests require the necessary toolboxes to be installed (see manual).


### Examples Test Stage

This stage tests all examples of CORA and can be triggered 
by ending the commit message with `--ci-run-examples`.

It runs the following test suites:

    runTestSuite('examples')


### Benchmarks Stage

This stage tests all benchmarks of CORA and can be triggered 
by ending the commit message with `--ci-run-benchmarks`.

    runTestSuite('benchmarks')


### Status Update Stage

As the CI is run on GitLab, users are unable to see on Bitbucket 
if the pipeline already finished or if more tests are scheduled.
This stage provides this information to the user including whether all tests were successful.

---

## Technical Details

The mirror from Bitbucket to GitLab runs every 30 minutes automatically, 
and after each pushed commit another mirror is requested.
However, as GitLab does not mirror the Bitbucket repository within 5 minutes after the last mirror,
it might take a few minutes for the CI to get ready...


