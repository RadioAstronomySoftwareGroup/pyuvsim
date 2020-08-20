<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Motivation and Context
<!--- Why is this change required? What problem does it solve? -->
<!--- If it fixes an open issue, please link to the issue here. If this PR closes an issue, put the word 'closes' before the issue link to auto-close the issue when the PR is merged. -->

## Types of changes
<!--- What types of changes does your code introduce? Put an replace the space with an `x` in all the boxes that apply: -->
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to change)
- [ ] Reference simulation update or replacement
- [ ] Documentation change (documentation changes only)
- [ ] Version change
- [ ] Build or continuous integration change


## Checklist:
<!--- You may remove the checklists that don't apply to your change type(s) or just leave them empty -->
<!--- Go over all the following points, and replace the space with an `x` in all the boxes that apply. -->
<!--- If you're unsure about any of these, don't hesitate to ask. We're here to help! -->

Bug fix checklist:
- [ ] My fix includes a new test that breaks as a result of the bug (if possible).
- [ ] All new and existing tests pass.
- [ ] I have checked that I reproduce the reference simulations or if there are differences they are explained below (if appropriate). If there are changes that are correct, I will update the reference simulation files after this PR is merged.
- [ ] I have checked (e.g., using the benchmarking tools) that this change does not significantly increase typical runtimes. If it does, I have included a justification in the comments on this PR.
- [ ] I have updated the [CHANGELOG](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/blob/master/CHANGELOG.md).

New feature checklist:
- [ ] I have added or updated the docstrings associated with my feature using the [numpy docstring format](https://numpydoc.readthedocs.io/en/latest/format.html).
- [ ] I have updated the documentation to highlight my new feature (if appropriate).
- [ ] I have added tests to cover my new feature.
- [ ] All new and existing tests pass.
- [ ] I have checked that I reproduce the reference simulations or if there are differences they are explained below (if appropriate). If there are changes that are correct, I will update the reference simulation files after this PR is merged.
- [ ] I have checked (e.g., using the benchmarking tools) that this change does not significantly increase typical runtimes. If it does, I have included a justification in the comments on this PR.
- [ ] I have updated the [CHANGELOG](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/blob/master/CHANGELOG.md).

Breaking change checklist:
- [ ] I have updated the docstrings associated with my change using the [numpy docstring format](https://numpydoc.readthedocs.io/en/latest/format.html).
- [ ] I have updated the documentation to reflect my changes (if appropriate).
- [ ] My change includes backwards compatibility and deprecation warnings (if possible).
- [ ] I have added tests to cover my changes.
- [ ] All new and existing tests pass.
- [ ] I have checked that I reproduce the reference simulations or if there are differences they are explained below (if appropriate). If there are changes that are correct, I will update the reference simulation files after this PR is merged.
- [ ] I have checked (e.g., using the benchmarking tools) that this change does not significantly increase typical runtimes. If it does, I have included a justification in the comments on this PR.
- [ ] I have updated the [CHANGELOG](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/blob/master/CHANGELOG.md).

Reference simulation update or replacement:
- [ ] If this is a new reference simulation or if the settings files have changed, I have added all the relevant files to the reference_simulations folder.
- [ ] If this is a new reference simulation I have fully documented it in a memo in docs and updated any other docs as appropriate.
- [ ] I have updated the [CHANGELOG](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/blob/master/CHANGELOG.md).

Documentation change checklist:
- [ ] Any updated docstrings use the [numpy docstring format](https://numpydoc.readthedocs.io/en/latest/format.html).
- [ ] If this is a significant change to the readme or other docs, I have checked that they are rendered properly on ReadTheDocs. (you may need help to get this branch to build on RTD, just ask!)

Version change checklist:
- [ ] I have updated the [CHANGELOG](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/blob/master/CHANGELOG.md) to put all the unreleased changes under the new version (leaving the unreleased section empty).
- [ ] I have run the appropriate tests for this level of version change (described in the readme under versioning).

Build or continuous integration change checklist:
- [ ] If required or optional dependencies have changed (including version numbers), I have updated the readme to reflect this.
- [ ] If this is a new CI setup, I have added the associated badge to the readme and to references/make_index.py (if appropriate).
