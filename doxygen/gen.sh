#!/bin/sh

set -e

PROJECT_NUMBER=$(git describe --tags --abbrev=0) doxygen susa.cfg