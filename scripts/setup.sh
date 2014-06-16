if [ -n "${STUDYNOISE}" ] ; then
echo "already set"

else
export THISDIR=`pwd`
export STUDYNOISE=${THISDIR}/
export PATH=${THISDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${THISDIR}/lib:${LD_LIBRARY_PATH}
fi
