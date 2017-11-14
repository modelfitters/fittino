if( INSTALL_PYTHON2 )

    set( miniconda "Miniconda2-4.3.11")

else()

    set( miniconda "Miniconda3-")

endif()

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64" )

    if( INSTALL_PYTHON2 )

        set( miniconda_url https://repo.continuum.io/miniconda/Miniconda2-4.3.11-MacOSX-x86_64.sh )
        set( miniconda_url_md5 d2edb7d4f3f55c35b9a25fd2d0ef6856 )

    else()

        set( miniconda_url https://repo.continuum.io/miniconda/Miniconda2-4.3.11-MacOSX-x86_64.sh )
        set( miniconda_url_md5 d2edb7d4f3f55c35b9a25fd2d0ef6856 )

    endif()

elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64" )

    if( INSTALL_PYTHON2 )

        set( miniconda_url https://repo.continuum.io/miniconda/Miniconda2-4.3.11-Linux-x86_64.sh )
        set( miniconda_url_md5 d573980fe3b5cdf80485add2466463f5 )

    else()

        set( miniconda_url https://repo.continuum.io/miniconda/Miniconda2-4.3.11-Linux-x86_64.sh )
        set( miniconda_url_md5 d573980fe3b5cdf80485add2466463f5 )

    endif()

elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86" )

    if( INSTALL_PYTHON2 )

        set( miniconda_url https://repo.continuum.io/miniconda/Miniconda2-4.3.11-Linux-x86.sh )
        set( miniconda_url_md5 09e9780288d29e46f3dd7eaa7ce280bc )

    else()

        set( miniconda_url https://repo.continuum.io/miniconda/Miniconda2-4.3.11-Linux-x86.sh )
        set( miniconda_url_md5 09e9780288d29e46f3dd7eaa7ce280bc )

    endif()

else()

    message( FATAL_ERROR "CMAKE_SYSTEM_NAME is ${CMAKE_SYSTEM_NAME} and CMAKE_HOST_SYSTEM_PROCESSOR is ${CMAKE_HOST_SYSTEM_PROCESSOR}." )

endif()

find_package( UnixCommands REQUIRED )

set( miniconda_install_command ${CMAKE_COMMAND} -E env HOME=<INSTALL_DIR> PYTHONNOUSERSITE=1 --unset=PYTHONPATH ${BASH} <DOWNLOADED_FILE> -b -f -p <INSTALL_DIR> )

externalproject_add(

    ${miniconda}
    URL ${miniconda_url}
    URL_MD5 ${miniconda_url_md5}
    CONFIGURE_COMMAND ""
    DOWNLOAD_NO_EXTRACT 1
    BUILD_COMMAND ""
    INSTALL_COMMAND ${miniconda_install_command}

)

externalproject_get_property( ${miniconda} install_dir )
set( miniconda_root_dir ${install_dir} )