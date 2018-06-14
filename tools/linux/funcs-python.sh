#!/bin/bash
# 
# This script will install Python
# 

SCRIPT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $SCRIPT_ROOT/funcs.sh


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


function python_full_bin()
{
    which $(python_bin)
}


function python_hasmodule()
{
    local ret=$(python_get "import sys;sys.stderr=sys.stdout;import ${1}")
    ret=$(echo ${ret}|wc -m)
    if [[ ${ret} == 1 ]]; then
        echo true
    else
        echo false
    fi
}


function python_info()
{
    cprint "Python version: $(python_version)"
    cprint "Python location: $(python_full_bin)"
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


function pip_install()
{
    if [[ $(install_systemwide) == true ]]; then
        $(pip_bin) install $@
    else
        $(pip_bin) install $@ --user
    fi
}


function pip_upgrade()
{
    if [[ $(install_systemwide) == true ]]; then
        $(pip_bin) install --upgrade $@
    else
        $(pip_bin) install --upgrade --user $@
    fi
}


function python_build_dependencies()
{
    require_build_essentials
    mapt-get "install libbz2-dev libsqlite3-dev libreadline-dev zlib1g-dev libncurses5-dev libssl-dev libgdbm-dev libssl-dev openssl tk-dev"
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

    require_web_essentials
    require_web_access

    cprint "Download python ..."
    mkdir -p python
    cd python

    local version=$(python_latest)
    local sourcedir=Python-${version}
    if [[ $(dryrun) == false && ! -d ${sourcedir} ]]; then
        python_download ${version}
        tar -xzf ${sourcedir}.tgz
    fi
    cd ${sourcedir}

    local prefix=$(project_opt)/python/${version}
    local prefixstr=$(project_optstr)/python/${version}
    if [[ ! -d ./build ]]; then

        cprint "Configure python for ${prefix} ..."
        if [[ $(dryrun) == false ]]; then
            python_build_dependencies

            mkdir -p ${prefix}
            ./configure --prefix=${prefix} \
                        --enable-shared \
                        --enable-optimizations \
                        --with-ensurepip=install
                        LDFLAGS=-Wl,-rpath=${prefix}/lib \
                        
        fi

        cprint "Build python ..."
        if [[ $(dryrun) == false ]]; then
            make -s -j2
        fi
    fi

    cprint "Install python in ${prefix} ..."
    if [[ $(dryrun) == false ]]; then
        mmakeinstall

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
    cprint "Verify python ${PYTHONVREQUEST} ..."

    # Try system installation
    if [[ $(cmdexists $(python_bin)) == false ]]; then
        mapt-get "install $(python_bin)"
    fi

    # Check version
    if [[ $(require_new_version_strict $(python_version) ${PYTHONVREQUEST}) == false ]]; then
        cprint "python version $(python_version) will be used"
        return
    fi

    # Install from source
    python_install_fromsource

    # Check version
    if [[ $(require_new_version_strict $(python_version) ${PYTHONVREQUEST}) == false ]]; then
        cprint "python version $(python_version) will be used"
    else
        if [[ $(cmdexists python) == false ]]; then
            cerror "python is not installed"
        else
            cerror "python version $(python_version) will be used but ${PYTHONVREQUEST} is required"
        fi
    fi
}


function require_pip()
{
    if [[ $(python_hasmodule "pip") == false ]]; then
        if [[ $(python_hasmodule "ensurepip") == true ]]; then
            if [[ $(install_systemwide) == true ]]; then
                $(python_bin) -m ensurepip
            else
                $(python_bin) -m ensurepip --user
            fi
        else
            if [[ $(python_bin) == "python3" ]];then
                mapt-get "install python3-pip"
            else
                mapt-get "install python-pip"
            fi
        fi
    fi

    pip_upgrade pip
    pip_upgrade setuptools
    pip_upgrade wheel
}


function require_pythondev()
{
    mapt-get "install $(python_bin)-dev"
}


function require_pyqt4()
{
    mapt-get "install $(python_bin)-qt4"
}


function require_pyqt5()
{
    mapt-get "install $(python_bin)-qt5"
}


function require_pyqt()
{
    require_pyqt4
    require_pyqt5
}



