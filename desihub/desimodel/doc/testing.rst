=================
Testing desimodel
=================

Introduction
------------

Tests using desimodel are a bit tricky since most of the code involves reading
data files that are not included with the git product.  In addition, the data
files can be rather large, both individually and as a package. This document
describes how to create lightweight test branches of the desimodel *data*
that can be used for rapid testing.

When to Create a Test Branch
----------------------------

**The word "branch(es)" below refers to svn branches not git branches,
unless otherwise noted.**

Test branches should track major changes to the code and data, but should
*not* change for minor bug fixes in the code.  For example, a branch named
'test-1.0' should work with all code that has a tag of the form 1.0.x.
Similarly a 'test-1.1' branch should work with all code that has a tag of the
form 1.1.x.

Test branches should be created immedately after the trunk has been tagged
to produce a new minor version.  For example, immedately after a svn tag
'1.0.0' is created, the 'test-1.0' branch should be created before there
are any additional changes to the trunkk.

How to Create a Test Branch
---------------------------

1. Create a branch in the standard way::

    base=https://desi.lbl.gov/svn/code/desimodel
    svn copy $base/trunk $base/branches/test-1.0

2. Check out the branch, if you did not create it in your own checkout::

    svn checkout $base/branches/test-1.0

3. Change to the branch directory::

    cd test-1.0

4. Run the trim code to create a separate datalite directory::

    python -c "from desimodel.trim import trim_data; trim_data('data', 'datalite')"

5. Add the datalite directory, remove data and commit::

    svn add datalite
    svn remove data
    svn commit -m "Adding lite files"

6. Rename the existing datalite directory::

    svn move datalite data
    svn commit -m "Rename datalite/ back to data/"

7. Now tests can get this lightweight branch on-the-fly with::

    svn export $base/branches/test-1.0

If you find any problems, just wipe out the branch and start again.
