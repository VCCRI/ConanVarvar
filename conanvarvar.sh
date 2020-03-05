# This file is part of ConanVarvar, a tool for detection of copy number variants
#
# Copyright (C) 2020 Victor Chang Cardiac Research Institute
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is provided without warranty of any kind.
#
# <http://www.gnu.org/licenses/>


# Abort on any error
set -e

# Command-line mode
if [ "$1" = "false" ]
then
    shift
    if [ "$1" = "/bin/sh" ]
    then
        Rscript --vanilla ConanVarvar.R --help
    else
        Rscript --vanilla ConanVarvar.R "$@"
    fi
fi

# GUI mode
if [ "$1" = "true" ]
then
    Rscript -e "options('shiny.port' = 3838, shiny.host = '0.0.0.0'); shiny::runApp('app.R')"
fi
