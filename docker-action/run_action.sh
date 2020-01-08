#!/bin/bash -x

source ${ROOTSYS}/bin/thisroot.sh

mkdir /opt/fittino

mkdir build
pushd build

cmake -DCMAKE_INSTALL_PREFIX=/opt/fittino ..
make
make doc
make install

popd

cp -rv /opt/fittino/share/doc/fittino/html .
chmod -R a+rwX ./html
