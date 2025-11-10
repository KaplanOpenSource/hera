#!/bin/bash
docker run --rm -it hera python -c "from hera.datalayer.project import getProjectList; getProjectList()"
