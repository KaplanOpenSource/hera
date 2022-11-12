.. _projectPage:

Project
*********

Overview
========

A project is a class of hera that is used to simplify the access to data.
Since all the documents of a class are related to a project, useng the
datalayer would require using the project name in every call. To simplify
the access, the project class includes a projectName property that eliminates
the need to call it in every call.

That is the code that uses the datalayer to retrieves all the 'measurement' documents of
the project 'my-project'  is

..  code-block:: python

    from hera import datalayer

    projectName = "my-project"
    List1 = datalayer.Measurements.getDocuments(projectName=projectName)


With the project class it is:


..  code-block:: python

    from hera import project

    projectName = "my-project"
    p = Project(projectName=projectName)

    List1 = p.Measurements.getDocuments()



or

..  code-block:: python

    from hera import project

    projectName = "my-project"
    p = Project(projectName=projectName)

    List1 = p.getMeasurementsDocuments()



Therefore, for multiple access to documents, using the project class is
simpler.

Interface
==========

The complete interface of the Project class is:








