#!/bin/bash
# 
# Build PyTMM
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh
source $SCRIPT_ROOT/funcs-python.sh


function pytmm_install_fromsource()
{
    local restorewd=$(pwd)

    cprint "Download PyTMM ..."
    mkdir -p PyTMM
    cd PyTMM

    if [[ $(dryrun) == false ]]; then
        requires_web_access
        if [[ ! -d db ]]; then
            git clone https://github.com/polyanskiy/refractiveindex.info-database db
        fi
        if [[ ! -d PyTMM ]]; then
            git clone https://github.com/kitchenknif/PyTMM PyTMM
        fi
    fi

    cprint "Install PyTMM ..."
    if [[ $(dryrun) == false ]]; then
        cd PyTMM

        $(python_bin) setup.py build -f
        $(python_bin) setup.py install
        
        cd ../db

        local pytmmdir=$(python_get "import PyTMM,os; print(os.path.dirname(PyTMM.__file__));")
        
        echo "Copy database to ${pytmmdir}/visirlib ..."
        cp -R database ${pytmmdir}/visirlib
        ls ${pytmmdir}/visirlib
    fi

    cd ${restorewd}
}


function require_pytmm()
{
    cprintstart
    cprint "Verify pytmm ..."

    # Check
    require_python

    if [[ $(python_hasmodule "PyTMM") == true ]]; then
        cprint "Python module \"PyTMM\" is installed"
        cprintend
        return
    fi

    # Install from source
    pytmm_install_fromsource

    if [[ $(python_hasmodule "PyTMM") == true ]]; then
        cprint "Python module \"PyTMM\" is installed"
    else
        cprint "Python module \"PyTMM\" is NOT installed"
    fi

    cprintend
}

