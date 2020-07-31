#!/bin/bash
# 
# Install fdmnes
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh
source ${SCRIPT_ROOT}/funcs-python.sh


function fdmnes_build_dependencies()
{
    require_pip
    pip_install ConfigParser
}


function fdmnes_install_fromsource()
{
    if [[ ! -d fdmnes && ${ARG_SKIPLONG} == true ]]; then
        cprint "Skipping fdmnes installation"
        return
    fi

    local restorewd=$(pwd)

    mkdir -p fdmnes
    cd fdmnes

    if [ ! -d /sware/exp/fdmnes ]; then
        cprint "Download fdmnes ..."

        require_web_essentials
        local fdmneslink=$(wget -O - -q http://neel.cnrs.fr/spip.php?article3137 | grep  -o 'http://neel.cnrs.fr/IMG/zip/[^"]*')
        local fdmneszipname=$(basename ${fdmneslink})

        if [[ ! -f ${fdmneszipname} && $(dryrun) == false ]]; then
            rm -f *.zip
            curl -O ${fdmneslink}
        fi
    fi

    cprint "Install fdmnes ..."
    local prefix=$(project_opt)/fdmnes
    local prefixstr=$(project_optstr)/fdmnes
    if [[ $(dryrun) == false ]]; then
        fdmnes_build_dependencies

        mexec mkdir -p ${prefix}

        local prefix_src=${prefix}/src
        local prefixstr_src=${prefixstr}/src

        if [[ ! -f ${prefix}/fdmnes ]]; then
            # Binary source directory
            if [[ -d /sware/exp/fdmnes ]]; then
                # Link to ESRF sware source
                if [ ! -d ${prefix_src} ]; then
                    mexec ln -s /sware/exp/fdmnes ${prefix_src}
                fi
            else
                # Unzip
                unzip -o ${fdmneszipname} -d ${prefix}
                mexec mv ${prefix}/fdmnes ${prefix_src}
            fi

            # Link to the binary that will be used
            mexec ln -s src/fdmnes_linux64 ${prefix}/fdmnes
            mexec chmod 775 ${prefix}/fdmnes
        fi

        # Environment
        addProfile $(project_resource) "# Installed fdmnes: ${prefixstr}"
        addBinPath ${prefix}
        addBinPathProfile $(project_resource) ${prefixstr}
        addProfile $(project_resource) "# Source of fdmnes: ${prefixstr_src}"
        addBinPath ${prefix_src}
        addBinPathProfile $(project_resource) ${prefixstr_src}
    fi

    cprint "Download pyfdmnes ..."
    if [[ $(dryrun) == false && ! -d pyFDMNES ]]; then
        git clone https://github.com/woutdenolf/pyFDMNES.git pyFDMNES
    fi

    cprint "Install pyfdmnes ..."
    if [[ $(dryrun) == false ]]; then
        cd pyFDMNES

        echo "[global]" > setup.cfg
        echo "; enter here the path to the fdmnes executable:" >> setup.cfg
        #echo "fdmnes_path=$(readlink -m ${prefix}/fdmnes)" >> setup.cfg
        echo "fdmnes_path=${prefix}/fdmnes" >> setup.cfg

        pip_install .
    fi

    cd ${restorewd}
}


function require_fdmnes()
{
    cprintstart "Require fdmnes ${1}"

    # Requirements (for running)
    require_python

    # Check
    if [[ $(python_hasmodule fdmnes) == true ]]; then
        cprint "Python module \"fdmnes\" is installed"
        cprintend "Require fdmnes ${1}"
        return
    fi

    # Install from source
    fdmnes_install_fromsource

    # Check
    if [[ $(python_hasmodule fdmnes) == true ]]; then
        cprint "Python module \"fdmnes\" is installed"
    else
        cprint "Python module \"fdmnes\" is NOT installed"
    fi

    cprintend "Require fdmnes ${1}"
}


