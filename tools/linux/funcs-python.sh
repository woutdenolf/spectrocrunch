#!/bin/bash
# 
# Functions related to Python.
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${SCRIPT_ROOT}/funcs-make.sh

if [[ -z ${PYTHONVREQUEST} ]];then
    PYTHONVREQUEST=""
fi


function python_exists()
{
    cmdexists $(python_bin)
}


function python_get()
{
    if [[ $(python_exists) == true ]]; then
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


function python_major_version()
{
    local tmp="import sys;print(sys.version_info[0])"
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


function python_depdir()
{
    echo dep_python$(python_full_version)
}


function python_include()
{
    python_get "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc());"
}


function python_lib()
{
    local _filename=$(python_get "import distutils.sysconfig;import os.path as op;print(op.join(distutils.sysconfig.get_config_var('LIBDIR'),distutils.sysconfig.get_config_var('LDLIBRARY')));")
    if [[ ! -f ${_filename} ]];then
        _filename=$(python_get "from distutils import sysconfig;import os.path as op;v = sysconfig.get_config_vars();fpaths = [op.join(v[pv], v['LDLIBRARY']) for pv in ('LIBDIR', 'LIBPL')]; print(list(filter(op.exists, fpaths))[0])")
    fi
    echo ${_filename}
}

function python_pkg_all()
{
    local tmp=$'\n'
    if [[ $(python_get "import site; print(hasattr(site,'getsitepackages'))") == "False" ]];then
        python_get "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib());"
    else
        python_get "import site${tmp}for p in site.getsitepackages():${tmp} print(p)"
    fi
}

function python_pkg()
{
    for x in $(python_pkg_all);do
        if [[ -d ${x} ]];then
            if [[ ! -z "$(ls -A ${x})" ]]; then
               echo ${x}
            fi
        fi
    done
}


function python_system_pkg()
{
    if [[ $(python_virtualenv_active) == true ]];then
        (deactivate;python_pkg)
    else
        python_pkg
    fi
}


function python_bin()
{
    if [[ -z ${PYTHONVREQUEST} ]]; then
        if [[ $(cmdexists python) == true ]]; then
            local tmp="import sys;t='{v[0]}.{v[1]}.{v[2]}'.format(v=list(sys.version_info[:3]));print(t)"
            PYTHONVREQUEST=$(python -c "$tmp")
            echo "python"
            return
        else
            PYTHONVREQUEST="3"
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
    # e.g. python-tk and python3.tk
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
    cprint "Python virtual environment: $(python_virtualenv_active)"
    cprint "Python location: $(python_full_bin)"
    cprint "Python package directories: $(python_pkg)"
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
    cprint "Pip: $(pip_source)"
}


function pip_exists()
{
    cmdexists $(pip_bin)
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
    if [[ $(python_get "import sys;print(sys.prefix!=getattr(sys,'base_prefix',sys.prefix))") == 'True' ]];then
        echo true
        return 
    fi
    if [[ $(python_get "import sys;print(sys.prefix!=getattr(sys,'real_prefix',sys.prefix))") == 'True' ]];then
        echo true
        return 
    fi
    echo false
}


function python_virtualenv_deactivate()
{
    if [[ $(python_virtualenv_active) == true ]]; then
        deactivate
    fi
}


function python_virtualenv_system_link()
{
    if [[ $(python_virtualenv_active) == true ]]; then
        local _restore=$(pwd)
        local _target=""
        for _pkgdir in $(python_pkg); do
            cd ${_pkgdir}
            for _syspkgdir in $(python_system_pkg); do
                echo "Look in ${_syspkgdir} ..."
                for var in "$@";do
                    echo "Look for ${_syspkgdir}/${var} ..."
                    for _target in ${_syspkgdir}/${var}; do
                        if [[ -e ${_target} || -d ${_target} ]];then
                            echo "$(pwd): create soft link ${_target}"
                            ln -s ${_target}
                        fi
                    done
                done
            done
        done
        cd ${_restore}
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
        $(pip_bin) install --user $@ 
    fi
}


function pip_upgrade()
{
    addProfilePythonUserBase
    if [[ $(install_systemwide) == true || $(python_virtualenv_active) == true ]]; then
        $(pip_bin) install --upgrade $@
    else
        $(pip_bin) install --upgrade --user $@ 
    fi
}


function pip_ensure()
{
    addProfilePythonUserBase
    if [[ $(python_hasmodule pip) == false ]]; then
        if [[ $(python_hasmodule ensurepip) == true ]]; then
            if [[ $(python_virtualenv_active) == true ]]; then
                $(python_bin) -m ensurepip
            else
                if [[ $(install_systemwide) == true ]]; then
                    $(python_bin) -m ensurepip
                else
                    $(python_bin) -m ensurepip --user
                fi
                # ensurepip is disabled in Debian/Ubuntu for the system python
                if [[ $(python_hasmodule pip) == false ]]; then
                    mapt-get install $(python_apt)-pip
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
    require_openssl
    mapt-get install libbz2-dev
    mapt-get install zlib1g-dev
    mapt-get install libsqlite3-dev
    mapt-get install libreadline-dev
    mapt-get install libncurses5-dev
    mapt-get install libssl-dev
    mapt-get install libgdbm-dev
    mapt-get install tk-dev libffi-dev
}


function python_url()
{
    echo "https://www.python.org/ftp/python/"
}


function python_all_versions()
{
    versions_from_site $(python_url) "href=\"[0-9\.]+"
}


function python_latest()
{
    latest_version python_all_versions ${1}
}


function python_download()
{
    if [[ ! -f ${1}.tar.gz ]]; then
        local _version=$(echo ${1} | grep -E -o "[\.0-9]+")
        curl -L $(python_url)/${_version}/Python-${_version}.tgz --output python-${_version}.tar.gz
    fi
}


function python_system_install()
{
    mapt-get install $(python_apt)
}


function python_source_install()
{
    source_install python "${1}" \
        --with-ensurepip=install \
        --enable-shared \
        LDFLAGS=-Wl,-rpath='"${prefix}/lib"'
    # --enable-shared (produce shared libraries and headers):
    #   needed by some libraries like xraylib
    # LDFLAGS=-Wl,-rpath=${prefix}/lib (linker argument "-rpath ${prefix}/lib"):
    #   put shared libraries in a custom location (avoid conflict with system python)
    # --enable-optimizations (profile guided optimizations):
    #   gives "profiling ... .gcda:Cannot open" warning on "import distutils"
}


function require_python()
{
    if [[ ! -z ${1} ]]; then
        PYTHONVREQUEST=${1}
    fi
    require_software python ${1}
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
    cprintstart "Require PyQt4"
    if [[ $(python_hasmodule "PyQt4") == false ]]; then
        mapt-get install $(python_apt)-qt4
    fi
    if [[ $(python_hasmodule "PyQt4") == false ]]; then
        python_virtualenv_system_link PyQt4 sip*.so sipconfig*.py
    fi
    if [[ $(python_hasmodule "PyQt4") == true ]]; then
        cprint "Python module \"PyQt4\" is working"
        cprintend "Require PyQt4"
        return
    else
        cprint "Python module \"PyQt4\" is NOT working. Try PySide ..."
    fi
    pip_install PySide
    if [[ $(python_hasmodule "PySide") == true ]]; then
        cprint "Python module \"PySide\" is working"
    else
        cprint "Python module \"PySide\" is NOT working"
    fi
    cprintend "Require PyQt4"
}


function require_pyqt5()
{
    cprintstart "Require PyQt5"
    if [[ $(python_hasmodule "PyQt5") == false ]]; then
        pip_install pyqt5
    fi
    if [[ $(python_hasmodule "PyQt5") == false ]]; then
        python_virtualenv_system_link PyQt5 sip*.so sipconfig*.py
    fi
    if [[ $(python_hasmodule "PyQt5") == true ]]; then
        cprint "Python module \"PyQt5\" is working"
        cprintend "Require PyQt5"
        return
    else
        cprint "Python module \"PyQt5\" is NOT working. Try PySide2 ..."
    fi
    pip_install PySide2
    if [[ $(python_hasmodule "PySide2") == true ]]; then
        cprint "Python module \"PySide2\" is working"
    else
        cprint "Python module \"PySide2\" is NOT working"
    fi
    cprintend "Require PyQt5"
}


function require_pyqt()
{
    if [[ $(python_major_version) == 2 ]];then
        require_pyqt4
    else
        require_pyqt5
    fi
}


function require_venv()
{
    if [[ $(python_major_version) == 3 ]];then
        if [[ $(python_hasmodule venv) == false ]];then
            mapt-get install $(python_apt)-venv
        fi
        if [[ $(python_hasmodule venv) == false ]];then
            pip_install virtualenv
        fi
    else
        if [[ $(python_hasmodule virtualenv) == false ]];then
            mapt-get install $(python_apt)-virtualenv
        fi
        if [[ $(python_hasmodule virtualenv) == false ]];then
            pip_install virtualenv
        fi
    fi
}


function create_venv()
{
    python_virtualenv_deactivate
    if [[ $(python_major_version) == 3 ]];then
        $(python_bin) -m venv ${1}
    else
        $(python_bin) -m virtualenv ${1}
    fi
}
