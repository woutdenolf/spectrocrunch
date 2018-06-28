#!/bin/bash
# 
# Functions related to Python.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs.sh


function python_get()
{
    if [[ $(cmdexists $(python_bin)) == true ]]; then
        $(python_bin) -c "$@"
    else
        echo ""
    fi
}


function python_version()
{
    local tmp="import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));print(t)"
    tmp=$(python_get "${tmp}") 
    if [[ -z ${tmp} ]];then
        tmp=-1
    fi
    echo $tmp
}


function python_full_version()
{
    local tmp="import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));print(t)"
    tmp=$(python_get "${tmp}") 
    if [[ -z ${tmp} ]];then
        tmp=-1
    fi
    echo $tmp
}


function python_include()
{
    python_get "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc());"
}


function python_lib()
{
    python_get "import distutils.sysconfig,os; print(os.path.join(distutils.sysconfig.get_config_var('LIBDIR'),distutils.sysconfig.get_config_var('LDLIBRARY')));"

}


function python_pkg()
{
    python_get "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib());"
}


function python_bin()
{
    if [[ -z ${PYTHONVREQUEST} ]]; then
        if [[ $(cmdexists python) == true ]]; then
            # Don't try python3 here because of virtual envs
            echo "python"
            return
        else
            PYTHONVREQUEST=3
        fi
    fi

    if [[ $(cmdexists "python${PYTHONVREQUEST:0:1}") == true ]]; then
        echo "python${PYTHONVREQUEST:0:1}"
    else
        echo "python"
    fi
}


function python_apt()
{
    strreplace $(python_bin) "python2" "python"
}


function python_full_bin()
{
    which $(python_bin)
}


function python_hasmodule()
{
    local restorewd=$(pwd)
    cd /
    local ret=$(python_get "import sys;sys.stderr=sys.stdout;import ${1}")
    cd ${restorewd}

    ret=$(echo ${ret}|wc -m)
    if [[ ${ret} == 1 ]]; then
        echo true
    else
        echo false
    fi
}


function python_info()
{
    cprint "Python version: $(python_full_version)"
    cprint "Python location: $(python_full_bin)"
    cprint "Python package directory: $(python_pkg)"
    cprint "Python include: $(python_include)"
    cprint "Python library: $(python_lib)"
}


function pip_version()
{
    $(pip_bin) --version | awk -v N=2 '{print $2}'
}


function pip_source()
{
    $(pip_bin) --version | awk '{$1= ""; print $0}'
}


function pip_info()
{
    cprint "Pip:$(pip_source)"
}


function pip_bin()
{
    if [[ $(python_hasmodule "pip") == true ]]; then
        echo "$(python_bin) -m pip"
        return
    fi

    if [[ -z ${PYTHONVREQUEST} ]]; then
        if [[ $(cmdexists pip) == true ]]; then
            echo "pip"
            return
        else
            PYTHONVREQUEST=3
        fi
    fi

    if [[ $(cmdexists "pip${PYTHONVREQUEST:0:1}") == true ]]; then
        echo "pip${PYTHONVREQUEST:0:1}"
    else
        echo "pip"
    fi
}


function python_virtualenv_active()
{
    if [[ "$VIRTUAL_ENV" != "" ]]; then
        echo true
    else
        echo false
    fi
}


function python_virtualenv_deactivate()
{
    if [[ $(python_virtualenv_active) == true ]]; then
        deactivate
    fi
}


function addProfilePythonUserBase()
{
    # Root under which packages/scripts are installed when no
    # permissions to the default package/script directory.
    # When virtualenv is active: ignores PYTHONUSERBASE
    local prefix=$(project_userbase)/python/$(python_full_version)
    local prefixstr=$(project_userbasestr)/python/$(python_full_version)
    addProfile $(project_resource) "# Python user base: ${prefixstr}"
    addVar "PYTHONUSERBASE" ${prefix}
    addVarProfile $(project_resource) "PYTHONUSERBASE" ${prefixstr}
    addBinPath ${prefix}/bin
    addBinPathProfile $(project_resource) ${prefixstr}/bin
    addLibPath ${prefix}/lib
    addLibPathProfile $(project_resource) ${prefixstr}/lib
}


function pip_install()
{
    addProfilePythonUserBase
    if [[ $(install_systemwide) == true || $(python_virtualenv_active) == true ]]; then
        $(pip_bin) install $@
    else
        $(pip_bin) install $@ --user
    fi
}


function pip_upgrade()
{
    addProfilePythonUserBase
    if [[ $(install_systemwide) == true || $(python_virtualenv_active) == true ]]; then
        $(pip_bin) install --upgrade $@
    else
        $(pip_bin) install --upgrade $@ --user
    fi
}


function pip_ensure()
{
    addProfilePythonUserBase
    if [[ $(python_hasmodule pip) == false ]]; then
        if [[ $(python_hasmodule ensurepip) == true ]]; then
            if [[ $(python_virtualenv_active) == true ]]; then
                if [[ $(install_systemwide) == true ]]; then
                    $(python_bin) -m ensurepip
                else
                    $(python_bin) -m ensurepip --user
                fi
            else
                if [[ $(install_systemwide) == true ]]; then
                    $(python_bin) -m ensurepip
                else
                    $(python_bin) -m ensurepip --user
                fi
            fi
        else
            mapt-get install $(python_apt)-pip
        fi
    fi
}



function pip_uninstall()
{
    $(pip_bin) uninstall --yes $@
}


function python_build_dependencies()
{
    require_build_essentials
    mapt-get install libsqlite3-dev libreadline-dev libncurses5-dev libssl-dev libgdbm-dev tk-dev libffi-dev
}


function python_url()
{
    echo "https://www.python.org/ftp/python/"
}


function _python_latest()
{
    curl -sL $(python_url) | grep -E "href=\"${PYTHONVREQUEST}" | while read f
    do
        version=$(echo ${f} | grep -o -E "[0-9\.]+[0-9]" | head -1)
        
        tmp=$(curl -sL $(python_url)/${version} | grep -o -E "Python-${version}.tgz" | wc -l)
        if [[ $tmp -ne 0 ]]; then
            echo ${version}
        fi
    done
}

function python_latest()
{
    _python_latest | tail -1
}


function python_download()
{
    curl -L $(python_url)/${1}/Python-${1}.tgz --output Python-${1}.tgz
}


function python_install_fromsource()
{
    local restorewd=$(pwd)

    cprint "Download python ..."
    mkdir -p python
    cd python

    local version=$(get_local_version_strict ${1})
    if [[ -z ${version} ]]; then
        require_web_essentials
        local version=$(python_latest)
    fi

    local sourcedir=Python-${version}
    if [[ $(dryrun) == false && ! -d ${sourcedir} ]]; then
        python_download ${version}
        tar -xzf ${sourcedir}.tgz
        rm -f ${sourcedir}.tgz
    fi
    cd ${sourcedir}

    local prefix=$(project_opt)/python/${version}
    local prefixstr=$(project_optstr)/python/${version}
    if [[ ! -d ./build ]]; then

        cprint "Configure python for ${prefix} ..."
        if [[ $(dryrun) == false ]]; then
            python_build_dependencies

            mexec mkdir -p ${prefix}
            ./configure --prefix=${prefix} \
                        --with-ensurepip=install \
                        --enable-shared \
                        LDFLAGS=-Wl,-rpath=${prefix}/lib
            # --enable-shared (produce shared libraries and headers): needed by xraylib
            # LDFLAGS=-Wl,-rpath=${prefix}/lib (linker argument "-rpath ${prefix}/lib"): put shared libraries in a customn location (avoid conflict with system python)
            # --enable-optimizations (profile guided optimizations): gives "profiling ... .gcda:Cannot open" warning on "import distutils"
        fi

        cprint "Build python ..."
        if [[ $(dryrun) == false ]]; then
            make -s -j2
        fi
    fi

    cprint "Install python in ${prefix} ..."
    if [[ $(dryrun) == false ]]; then
        if [[ ! -f ${prefix}/bin/python ]]; then
            mmakeinstall python-${version}
        fi

        addProfile $(project_resource) "# Installed python: ${prefixstr}"
        addBinPath ${prefix}/bin
        addBinPathProfile $(project_resource) "${prefixstr}/bin"
        addLibPath ${prefix}/lib
        addLibPathProfile $(project_resource) "${prefixstr}/lib"
    fi

    cd ${restorewd}
}


function require_python()
{
    if [[ ! -z ${1} ]]; then
        PYTHONVREQUEST=${1}
    fi
    if [[ -z ${PYTHONVREQUEST} ]]; then
        cprint "Verify python (no specific version requested) ..."
    else
        cprint "Verify python ${PYTHONVREQUEST} ..."
    fi

    # Try system installation
    if [[ $(cmdexists $(python_bin)) == false ]]; then
        mapt-get install $(python_apt)
    fi

    # Check version
    if [[ $(require_new_version_strict $(python_full_version) ${PYTHONVREQUEST}) == false ]]; then
        cprint "Python version $(python_full_version) is used"
        return
    fi

    # Install from source
    python_install_fromsource ${PYTHONVREQUEST}

    # Check version
    if [[ $(require_new_version_strict $(python_full_version) ${PYTHONVREQUEST}) == false ]]; then
        cprint "Python version $(python_full_version) is used"
    else
        if [[ $(cmdexists python) == false ]]; then
            cerror "Python is not installed"
        else
            cerror "Python version $(python_full_version) is used but ${PYTHONVREQUEST} is required"
        fi
        return 1
    fi
}


function require_pip()
{
    pip_ensure
    pip_upgrade pip
    pip_upgrade setuptools
    pip_upgrade wheel
}


function require_pythondev()
{
    mapt-get install $(python_apt)-dev
}


function require_pyqt4()
{
    mapt-get install $(python_apt)-pyqt4
}


function require_pyqt5()
{
    mapt-get install $(python_apt)-pyqt5
    pip_install pyqt5
}


function require_pyqt()
{
    require_pyqt4
    require_pyqt5
}



